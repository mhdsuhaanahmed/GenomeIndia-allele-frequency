"""
sensitivity_analysis.py

Tests whether the conclusions depend on the assumed Indian allele number.

The main analysis assumes AN_india = 2 * 9768 * GENOTYPE_RATE with
GENOTYPE_RATE = 0.98, the consortium's stated minimum. That is a floor, not a
measurement: the true per-site genotyping rate is somewhere between 0.98 and
1.00 and is not reported per variant. This script re-runs every comparison at
the upper bound (GENOTYPE_RATE = 1.0) and reports what changes.

Nothing on the gnomAD side changes, so the only thing moving is the precision
of the Indian arm. Point estimates of the odds ratio are almost invariant;
what can move is the width of the CI and, at the margin, FDR significance.

Input : outputs/stats/statistical_comparison_v2.csv
Output: outputs/stats/sensitivity_genotype_rate.csv  (+ printed summary)
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE = r"D:\GENOMEINDIA"
STATS_DIR = os.path.join(BASE, "outputs", "stats")
IN_CSV = os.path.join(STATS_DIR, "statistical_comparison_v2.csv")
OUT_CSV = os.path.join(STATS_DIR, "sensitivity_genotype_rate.csv")

N_INDIVIDUALS = 9768
RATES = [0.98, 1.00]          # baseline, upper bound
BASELINE_RATE = 0.98
ALPHA = 0.05

GNOMAD_URL = "https://gnomad.broadinstitute.org/api"
POP_IDS = {
    "South Asian": "sas",
    "East Asian": "eas",
    "African": "afr",
    "European": "nfe",
    "Admixed American": "amr",
    "Middle Eastern": "mid",
}


# --------------------------------------------------------------------------
# population allele counts: reuse if already saved, otherwise re-fetch
# --------------------------------------------------------------------------

def fetch_pop_counts(variant_ids):
    """Fallback only. Re-queries gnomAD v4 for AC/AN by population."""
    import requests

    query = """
    query($vid: String!) {
      variant(variantId: $vid, dataset: gnomad_r4) {
        genome { populations { id ac an } }
        exome  { populations { id ac an } }
      }
    }
    """
    rows = []
    for vid in variant_ids:
        try:
            r = requests.post(GNOMAD_URL, json={"query": query,
                                                "variables": {"vid": vid}},
                              timeout=30)
            data = r.json().get("data", {}).get("variant") or {}
        except Exception as e:                      # network, not absence
            print(f"  ! {vid}: {e}", file=sys.stderr)
            continue

        pooled = {}
        for source in ("exome", "genome"):
            block = data.get(source) or {}
            for p in block.get("populations", []):
                pid = p["id"]
                ac, an = pooled.get(pid, (0, 0))
                pooled[pid] = (ac + p["ac"], an + p["an"])

        for label, pid in POP_IDS.items():
            if pid in pooled:
                ac, an = pooled[pid]
                rows.append({"variant_id": vid, "population": label,
                             "ac_pop": ac, "an_pop": an})
        time.sleep(0.3)
    return pd.DataFrame(rows)


def load_panel():
    if not os.path.exists(IN_CSV):
        sys.exit(f"Not found: {IN_CSV}\nRun statistical_analysis_v2.py first.")

    df = pd.read_csv(IN_CSV)

    if {"ac_pop", "an_pop"}.issubset(df.columns):
        print("Using AC/AN already stored in statistical_comparison_v2.csv")
        return df

    print("AC/AN not stored in the CSV; re-querying gnomAD v4 ...")
    if "variant_id" not in df.columns:
        sys.exit("No variant_id column either. Re-run the main script with "
                 "variant_id and ac_pop/an_pop written out.")
    counts = fetch_pop_counts(sorted(df["variant_id"].unique()))
    if counts.empty:
        sys.exit("gnomAD returned nothing. Check the connection and retry; "
                 "an empty result here is a timeout, not an absence.")
    merged = df.merge(counts, on=["variant_id", "population"], how="left")
    missing = merged["ac_pop"].isna().sum()
    if missing:
        print(f"  ! {missing} comparisons unmatched and dropped")
    return merged.dropna(subset=["ac_pop", "an_pop"])


# --------------------------------------------------------------------------
# recompute at a given genotyping rate
# --------------------------------------------------------------------------

def run_at_rate(df, rate):
    an_india = int(round(2 * N_INDIVIDUALS * rate))
    out = []

    for _, row in df.iterrows():
        ac_i = int(round(row["af_india"] * an_india))
        ac_p, an_p = int(row["ac_pop"]), int(row["an_pop"])

        a, b = ac_i, an_india - ac_i          # India: alt, ref
        c, d = ac_p, an_p - ac_p              # comparison: alt, ref

        table = np.array([[a, b], [c, d]])
        _, p = stats.fisher_exact(table)

        # Haldane-Anscombe correction only for the CI, so zero cells
        # do not produce an undefined interval
        a_, b_, c_, d_ = (a + 0.5, b + 0.5, c + 0.5, d + 0.5) if 0 in (a, b, c, d) \
            else (a, b, c, d)
        or_point = (a_ * d_) / (b_ * c_)
        se = np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
        lo = np.exp(np.log(or_point) - 1.96 * se)
        hi = np.exp(np.log(or_point) + 1.96 * se)

        out.append({
            "gene": row["gene"],
            "rsid": row["rsid"],
            "category": row.get("category", ""),
            "population": row["population"],
            "genotype_rate": rate,
            "an_india": an_india,
            "af_india": row["af_india"],
            "af_pop": row["af_pop"],
            "odds_ratio": or_point,
            "or_lower": lo,
            "or_upper": hi,
            "ci_width_log": np.log(hi) - np.log(lo),
            "p_value": p,
        })

    res = pd.DataFrame(out)
    rej, q, _, _ = multipletests(res["p_value"], alpha=ALPHA, method="fdr_bh")
    res["q_value"] = q
    res["significant_fdr"] = rej
    res["direction"] = np.where(res["af_india"] > res["af_pop"], "higher", "lower")
    return res


# --------------------------------------------------------------------------

def main():
    df = load_panel()
    print(f"\nComparisons loaded: {len(df)} "
          f"({df['gene'].nunique()} loci x {df['population'].nunique()} populations)\n")

    runs = {rate: run_at_rate(df, rate) for rate in RATES}
    for rate, res in runs.items():
        an = int(round(2 * N_INDIVIDUALS * rate))
        n_sig = int(res["significant_fdr"].sum())
        print(f"GENOTYPE_RATE = {rate:.2f}   AN_india = {an}   "
              f"significant after FDR: {n_sig}/{len(res)}")

    key = ["gene", "rsid", "population"]
    base = runs[BASELINE_RATE].set_index(key)
    alt_rate = [r for r in RATES if r != BASELINE_RATE][0]
    alt = runs[alt_rate].set_index(key)

    cmp = pd.DataFrame({
        "category": base["category"],
        "or_base": base["odds_ratio"],
        "or_alt": alt["odds_ratio"],
        "q_base": base["q_value"],
        "q_alt": alt["q_value"],
        "sig_base": base["significant_fdr"],
        "sig_alt": alt["significant_fdr"],
    })
    cmp["pct_change_or"] = 100 * (cmp["or_alt"] - cmp["or_base"]).abs() / cmp["or_base"]
    cmp["log_or_shift"] = (np.log(cmp["or_alt"]) - np.log(cmp["or_base"])).abs()
    cmp["status_changed"] = cmp["sig_base"] != cmp["sig_alt"]

    print("\n" + "=" * 70)
    print("SENSITIVITY TO THE INDIAN ALLELE-NUMBER ASSUMPTION")
    print("=" * 70)
    print(f"Max change in odds ratio      : {cmp['pct_change_or'].max():.3f} %")
    print(f"Median change in odds ratio   : {cmp['pct_change_or'].median():.3f} %")
    print(f"Max shift in log odds ratio   : {cmp['log_or_shift'].max():.5f}")
    print(f"Comparisons changing FDR status: {int(cmp['status_changed'].sum())}")

    changed = cmp[cmp["status_changed"]]
    if len(changed):
        print("\nThese comparisons changed significance status:")
        print(changed[["or_base", "or_alt", "q_base", "q_alt",
                       "sig_base", "sig_alt"]].to_string())
        print("\nEach of these must be named explicitly in Section IV-E.")
    else:
        print("\nNo comparison changes significance status between the bounds.")

    # does the central claim survive?
    print("\nPER-LOCUS ELEVATION (significantly higher in >= 4 of 6 populations)")
    for rate, res in runs.items():
        sig = res[res["significant_fdr"]]
        elevated = (sig.groupby(["category", "gene"])["direction"]
                      .apply(lambda s: (s == "higher").sum() >= 4))
        by_cat = elevated.groupby(level=0).agg(["sum", "count"])
        line = "   ".join(f"{c}: {int(r['sum'])}/{int(r['count'])}"
                          for c, r in by_cat.iterrows())
        print(f"  rate {rate:.2f}  ->  {line}")

    os.makedirs(STATS_DIR, exist_ok=True)
    cmp.reset_index().to_csv(OUT_CSV, index=False)
    print(f"\nSaved to {OUT_CSV}")


if __name__ == "__main__":
    main()