"""
sensitivity_analysis.py

Tests whether the conclusions depend on the assumed Indian allele number.

The main analysis assumes AN_india = 2 * 9768 * GENOTYPE_RATE with
GENOTYPE_RATE = 0.98, the consortium's stated minimum. That is a floor, not a
measurement: the true per-site genotyping rate lies somewhere between 0.98 and
1.00 and is not reported per variant. This script re-runs every comparison at
the upper bound and reports what changes.

Nothing on the gnomAD side changes, so the only thing moving is the precision
of the Indian arm. Point estimates are near-invariant; what can move is the
width of the CI and, at the margin, FDR significance.

Loci are keyed on (gene, rsid), never on gene alone, because CYP2C9*2 and
CYP2C9*3 share a gene name and their opposite deviations are a result of the
study. A gene-level groupby silently merges them.

Input : outputs/stats/statistical_comparison_v2.csv
Output: outputs/stats/sensitivity_genotype_rate.csv
        outputs/stats/sensitivity_genotype_rate_per_locus.csv
"""

import os
import sys
import time
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

BASE      = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATS_DIR = os.path.join(BASE, "outputs", "stats")
IN_CSV    = os.path.join(STATS_DIR, "statistical_comparison_v2.csv")
OUT_CSV   = os.path.join(STATS_DIR, "sensitivity_genotype_rate.csv")

N_INDIVIDUALS = 9768
RATES         = [0.98, 1.00]      # baseline, upper bound
BASELINE_RATE = 0.98
ALPHA         = 0.05
MAJORITY      = 4                 # elevated = significantly higher in >= 4 of 6

LOCUS_KEY = ["gene", "rsid"]

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
            r = requests.post(GNOMAD_URL,
                              json={"query": query, "variables": {"vid": vid}},
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

    for col in LOCUS_KEY + ["population", "af_india", "af_pop"]:
        if col not in df.columns:
            sys.exit(f"Required column missing from {IN_CSV}: {col}")

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
    missing = int(merged["ac_pop"].isna().sum())
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

        _, p = stats.fisher_exact(np.array([[a, b], [c, d]]))

        # Haldane-Anscombe correction for the interval only, so a zero cell
        # does not give an undefined CI. The Fisher p above is exact and
        # uncorrected.
        if 0 in (a, b, c, d):
            a_, b_, c_, d_ = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        else:
            a_, b_, c_, d_ = a, b, c, d

        or_point = (a_ * d_) / (b_ * c_)
        se = np.sqrt(1 / a_ + 1 / b_ + 1 / c_ + 1 / d_)
        lo = float(np.exp(np.log(or_point) - 1.96 * se))
        hi = float(np.exp(np.log(or_point) + 1.96 * se))

        out.append({
            "gene": row["gene"],
            "rsid": row["rsid"],
            "category": row.get("category", ""),
            "population": row["population"],
            "genotype_rate": rate,
            "an_india": an_india,
            "ac_india": ac_i,
            "af_india": row["af_india"],
            "af_pop": row["af_pop"],
            "odds_ratio": or_point,
            "or_lower": lo,
            "or_upper": hi,
            "ci_width_log": float(np.log(hi) - np.log(lo)),
            "p_value": p,
        })

    res = pd.DataFrame(out)
    rej, q, _, _ = multipletests(res["p_value"], alpha=ALPHA, method="fdr_bh")
    res["q_value"] = q
    res["significant_fdr"] = rej
    res["direction"] = np.where(res["af_india"] > res["af_pop"],
                                "higher", "lower")
    return res


def locus_elevation(res):
    """Per (gene, rsid): populations where India is significantly higher,
    out of populations tested, and whether that clears the majority rule."""
    rows = []
    for (gene, rsid), grp in res.groupby(LOCUS_KEY, sort=False):
        sig = grp[grp["significant_fdr"]]
        n_hi = int((sig["direction"] == "higher").sum())
        n_lo = int((sig["direction"] == "lower").sum())
        rows.append({
            "category": grp["category"].iloc[0],
            "gene": gene,
            "rsid": rsid,
            "n_higher": n_hi,
            "n_lower": n_lo,
            "n_tested": len(grp),
            "elevated": n_hi >= MAJORITY,
        })
    return pd.DataFrame(rows).sort_values(["category", "gene", "rsid"])


# --------------------------------------------------------------------------

def main():
    df = load_panel()
    n_loci = df.groupby(LOCUS_KEY).ngroups
    print(f"\nComparisons loaded: {len(df)} "
          f"({n_loci} loci x {df['population'].nunique()} populations)\n")

    runs = {rate: run_at_rate(df, rate) for rate in RATES}
    for rate, res in runs.items():
        an = int(round(2 * N_INDIVIDUALS * rate))
        print(f"GENOTYPE_RATE = {rate:.2f}   AN_india = {an}   "
              f"significant after FDR: "
              f"{int(res['significant_fdr'].sum())}/{len(res)}")

    key = LOCUS_KEY + ["population"]
    base = runs[BASELINE_RATE].set_index(key)
    alt_rate = [r for r in RATES if r != BASELINE_RATE][0]
    alt = runs[alt_rate].set_index(key).reindex(base.index)

    cmp = pd.DataFrame({
        "category": base["category"],
        "or_base": base["odds_ratio"],
        "or_alt": alt["odds_ratio"],
        "ci_width_base": base["ci_width_log"],
        "ci_width_alt": alt["ci_width_log"],
        "q_base": base["q_value"],
        "q_alt": alt["q_value"],
        "sig_base": base["significant_fdr"],
        "sig_alt": alt["significant_fdr"],
    })
    cmp["pct_change_or"] = (100 * (cmp["or_alt"] - cmp["or_base"]).abs()
                            / cmp["or_base"])
    cmp["log_or_shift"] = (np.log(cmp["or_alt"]) - np.log(cmp["or_base"])).abs()
    cmp["status_changed"] = cmp["sig_base"] != cmp["sig_alt"]

    print("\n" + "=" * 70)
    print("SENSITIVITY TO THE INDIAN ALLELE-NUMBER ASSUMPTION")
    print("=" * 70)
    print(f"Max change in odds ratio        : {cmp['pct_change_or'].max():.3f} %")
    print(f"Median change in odds ratio     : {cmp['pct_change_or'].median():.3f} %")
    print(f"Max shift in log odds ratio     : {cmp['log_or_shift'].max():.5f}")
    print(f"Comparisons changing FDR status : {int(cmp['status_changed'].sum())}")

    worst = cmp["pct_change_or"].idxmax()
    print(f"Largest shift at                : {worst[0]} {worst[1]} vs {worst[2]}")

    changed = cmp[cmp["status_changed"]]
    if len(changed):
        print("\nThese comparisons changed significance status:")
        print(changed[["or_base", "or_alt", "q_base", "q_alt",
                       "sig_base", "sig_alt"]].to_string())
        print("\nEach of these must be named explicitly in Section IV-G.")
    else:
        print("\nNo comparison changes significance status between the bounds.")

    # ---- per-locus, keyed on rsid so CYP2C9*2 and *3 stay separate ----
    print("\nPER-LOCUS ELEVATION "
          f"(significantly higher in >= {MAJORITY} of 6 populations)")
    elev = {}
    for rate, res in runs.items():
        e = locus_elevation(res)
        elev[rate] = e
        print(f"\n  GENOTYPE_RATE = {rate:.2f}")
        for _, r in e.iterrows():
            flag = "   ELEVATED" if r["elevated"] else ""
            print(f"    [{r['category']:<8}] {r['gene']:<9} {r['rsid']:<12} "
                  f"higher {r['n_higher']}/{r['n_tested']}  "
                  f"lower {r['n_lower']}/{r['n_tested']}{flag}")
        by_cat = e.groupby("category")["elevated"].agg(["sum", "count"])
        line = "   ".join(f"{c}: {int(v['sum'])}/{int(v['count'])}"
                          for c, v in by_cat.iterrows())
        print(f"    -> {line}")

    # does any locus cross the majority threshold between the bounds?
    a = elev[BASELINE_RATE].set_index(LOCUS_KEY)["elevated"]
    b = elev[alt_rate].set_index(LOCUS_KEY)["elevated"].reindex(a.index)
    flipped = a.index[a.values != b.values]
    print(f"\n  Loci crossing the majority threshold between bounds: {len(flipped)}")
    for gene, rsid in flipped:
        print(f"    {gene} {rsid}: {bool(a.loc[(gene, rsid)])} -> "
              f"{bool(b.loc[(gene, rsid)])}")

    os.makedirs(STATS_DIR, exist_ok=True)
    cmp.reset_index().to_csv(OUT_CSV, index=False)
    elev_path = OUT_CSV.replace(".csv", "_per_locus.csv")
    pd.concat([e.assign(genotype_rate=r) for r, e in elev.items()]) \
      .to_csv(elev_path, index=False)
    print(f"\nSaved to {OUT_CSV}")
    print(f"Saved to {elev_path}")


if __name__ == "__main__":
    main()