import requests
import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import time
import json
import os

PANEL_CSV = r"D:\GENOMEINDIA\outputs\multi_gene\full_panel.csv"
CACHE     = r"D:\GENOMEINDIA\outputs\multi_gene\gnomad_counts_cache_v2.json"
OUT_DIR   = r"D:\GENOMEINDIA\outputs\stats"
os.makedirs(OUT_DIR, exist_ok=True)

N_INDIVIDUALS = 9768
GENOTYPE_RATE = 0.98          # set to 1.0 for the sensitivity analysis

POPS = {"sas": "South Asian", "eas": "East Asian", "afr": "African",
        "nfe": "European", "amr": "Admixed American", "mid": "Middle Eastern"}


def fetch_counts(variant_id, max_retries=5):
    q = """query($vid: String!, $ds: DatasetId!) {
      variant(variantId: $vid, dataset: $ds) {
        genome { ac an populations { id ac an } } } }"""
    for _ in range(max_retries):
        try:
            r = requests.post("https://gnomad.broadinstitute.org/api",
                              json={"query": q, "variables": {"vid": variant_id, "ds": "gnomad_r4"}},
                              timeout=40)
            if r.status_code == 200:
                p = r.json()
                if p.get("errors"):
                    return None
                v = p.get("data", {}).get("variant")
                if v and v.get("genome"):
                    return v["genome"]
        except requests.exceptions.RequestException:
            pass
        time.sleep(3)
    return None


def get_all_counts(df):
    cache = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    for _, row in df.iterrows():
        vid = row["variant_id"]
        if vid in cache:
            continue
        print(f"  fetching {row['gene']} ({vid})")
        g = fetch_counts(vid)
        if g:
            cache[vid] = g
            json.dump(cache, open(CACHE, "w"), indent=2)
        else:
            print(f"    FAILED — {row['gene']} dropped from statistics")
        time.sleep(1)
    return cache


def wilson_ci(ac, an, z=1.96):
    if an == 0:
        return (np.nan, np.nan)
    p = ac / an
    d = 1 + z**2 / an
    c = (p + z**2 / (2*an)) / d
    h = (z/d) * np.sqrt(p*(1-p)/an + z**2/(4*an**2))
    return (max(0.0, c-h), min(1.0, c+h))


def odds_ratio_ci(a, b, c, d):
    if min(a, b, c, d) == 0:
        a, b, c, d = a+0.5, b+0.5, c+0.5, d+0.5
    orv = (a*d)/(b*c)
    se = np.sqrt(1/a + 1/b + 1/c + 1/d)
    return orv, np.exp(np.log(orv) - 1.96*se), np.exp(np.log(orv) + 1.96*se)


def hudson_fst(p1, n1, p2, n2):
    if n1 < 2 or n2 < 2:
        return np.nan, np.nan
    num = (p1-p2)**2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
    den = p1*(1-p2) + p2*(1-p1)
    return num, den


if __name__ == "__main__":
    df = pd.read_csv(PANEL_CSV)
    print(f"Panel: {len(df)} loci — {df['gene'].tolist()}\n")

    print("Fetching allele counts...")
    cache = get_all_counts(df)

    an_india = int(round(2 * N_INDIVIDUALS * GENOTYPE_RATE))
    print(f"\nIndian allele number: {an_india} "
          f"(2 x {N_INDIVIDUALS} x {GENOTYPE_RATE:.0%})\n")

    rows, fst_rows = [], []

    for _, r in df.iterrows():
        vid = r["variant_id"]
        if vid not in cache:
            continue
        g = cache[vid]
        pop_counts = {p["id"]: (p["ac"], p["an"]) for p in g.get("populations", [])}

        ac_india = int(round(r["india_af"] * an_india))
        lo_i, hi_i = wilson_ci(ac_india, an_india)

        for pid, pname in POPS.items():
            if pid not in pop_counts:
                continue
            ac_pop, an_pop = pop_counts[pid]
            if an_pop < 100:
                continue

            # Fisher exact rather than chi-square: several loci in this panel
            # have expected cell counts far below the chi-square approximation's
            # validity threshold (ALDH2 is one Indian alternate copy; DPYD is
            # zero against East Asian), and those are exactly the comparisons
            # the validations rest on.
            tbl = [[ac_india, an_india - ac_india], [ac_pop, an_pop - ac_pop]]
            _, p_val = stats.fisher_exact(tbl)

            orv, or_lo, or_hi = odds_ratio_ci(ac_india, an_india - ac_india,
                                              ac_pop, an_pop - ac_pop)

            af_i, af_p = ac_india/an_india, ac_pop/an_pop
            lo_p, hi_p = wilson_ci(ac_pop, an_pop)
            num, den = hudson_fst(af_i, an_india, af_p, an_pop)

            fst_rows.append({"gene": r["gene"], "category": r["category"],
                             "population": pname, "num": num, "den": den})

            rows.append({
                "gene": r["gene"], "rsid": r["rsid"], "phenotype": r["phenotype"],
                "category": r["category"], "population": pname,
                "af_india": round(af_i, 6), "af_india_ci": f"{lo_i:.4f}-{hi_i:.4f}",
                "af_pop": round(af_p, 6), "af_pop_ci": f"{lo_p:.4f}-{hi_p:.4f}",
                "ac_india": ac_india, "an_india": an_india,
                "ac_pop": ac_pop, "an_pop": an_pop,
                "odds_ratio": round(orv, 4), "or_95ci": f"{or_lo:.3f}-{or_hi:.3f}",
                "or_lower": or_lo, "or_upper": or_hi,
                "p_value": p_val,
                "fst_locus": round(num/den, 5) if den else np.nan,
            })

    res = pd.DataFrame(rows)
    reject, q, _, _ = multipletests(res["p_value"], method="fdr_bh")
    res["q_value"] = q
    res["significant_fdr"] = reject
    res["direction"] = np.where(res["odds_ratio"] > 1, "higher", "lower")
    res = res.sort_values(["category", "gene", "population"])
    res.to_csv(os.path.join(OUT_DIR, "statistical_comparison_v2.csv"), index=False)

    print("="*74)
    print(f"TESTS: {len(res)}   significant after FDR: {res['significant_fdr'].sum()}")
    print("="*74)

    # --- the key analysis: does direction differ by category? ---
    print("\nDIRECTION BY CATEGORY (FDR-significant comparisons only)")
    sig = res[res["significant_fdr"]]
    for cat, grp in sig.groupby("category"):
        hi = (grp["direction"] == "higher").sum()
        lo = (grp["direction"] == "lower").sum()
        print(f"  {cat:<10} higher {hi:>3}   lower {lo:>3}   ({len(grp)} tests)")

    # locus-level: is the locus elevated vs the majority of populations?
    print("\nPER-LOCUS: populations where India is significantly higher")
    for (cat, gene, rsid), grp in res.groupby(["category", "gene", "rsid"]):
        s = grp[grp["significant_fdr"]]
        hi = (s["direction"] == "higher").sum()
        print(f"  [{cat:<8}] {gene:<9} {rsid:<12} higher in {hi}/{len(grp)} populations")

    # Fisher test on the category split, counting loci not comparisons
    locus = (res[res["significant_fdr"]]
             .groupby(["category", "gene", "rsid"])["direction"]
             .apply(lambda s: (s == "higher").sum() >= 4)   # majority of 6 populations
             .reset_index(name="mostly_higher"))
    dis = locus[locus["category"] == "disease"]
    pgx = locus[locus["category"] == "pgx"]
    if len(dis) and len(pgx):
        table = [[pgx["mostly_higher"].sum(), (~pgx["mostly_higher"]).sum()],
                 [dis["mostly_higher"].sum(), (~dis["mostly_higher"]).sum()]]
        odds, pf = stats.fisher_exact(table)
        print(f"\nFisher exact, loci elevated vs majority of populations:")
        print(f"  pgx     {table[0][0]}/{sum(table[0])}")
        print(f"  disease {table[1][0]}/{sum(table[1])}")
        print(f"  p = {pf:.4f}   (n is small; this is descriptive, not confirmatory)")

    print("\nStrongest effects:")
    top = res.reindex(res["odds_ratio"].apply(lambda x: abs(np.log(x))).sort_values(ascending=False).index)
    print(top.head(15)[["gene", "category", "population", "af_india", "af_pop",
                        "odds_ratio", "or_95ci", "q_value"]].to_string(index=False))

    fst = pd.DataFrame(fst_rows).dropna()
    print("\nPanel-wide Hudson FST (ratio of averages):")
    for pop, grp in fst.groupby("population"):
        d = grp["den"].sum()
        if d > 0:
            print(f"  {pop:<20} {grp['num'].sum()/d:.5f}  ({len(grp)} loci)")

    print(f"\nSaved to {OUT_DIR}\\statistical_comparison_v2.csv")