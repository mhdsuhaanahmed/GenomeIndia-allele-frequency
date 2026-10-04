"""
multi_gene_comparison.py

Step one of the pipeline. Resolves each panel locus to a GRCh38 position,
validates allele identity between GenomeIndia and gnomAD v4, and writes the
validated panel that every downstream script reads.

Output: outputs/multi_gene/full_panel.csv

Two distinctions this script exists to enforce:

1. A network failure and a genuine absence are different events. A run with any
   network failure writes to full_panel_INCOMPLETE.csv rather than overwriting
   a good panel.

2. The allele named in the clinical literature and the alternate allele in the
   genomic representation are not always the same. On minus-strand genes they
   are complements — MTHFR C677T is G>A genomically, DPYD c.1905+1G>A is C>T.
   risk_allele records the literature label; expected_alt records what must
   actually be matched in gnomAD. Conflating them produces a false mismatch.
"""

import os
import sys
import time
import requests
import pandas as pd

BASE     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(BASE, "data", "9768GI_SummaryStats")
OUT_DIR  = os.path.join(BASE, "outputs", "multi_gene")
os.makedirs(OUT_DIR, exist_ok=True)

REQUEST_PAUSE = 1.5        # seconds between loci, applied win or lose
BACKOFF_BASE  = 4          # seconds; 4, 8, 16, 32 — gnomAD's rate-limit window
                           # is longer than a 1-2-4 ramp survives

GENE_PANEL = [
    # positive control
    {"rsid": "rs671",      "gene": "ALDH2",   "phenotype": "Alcohol flush / ESCC",
     "category": "control", "risk_allele": "A", "expected_alt": "A"},

    # disease susceptibility
    {"rsid": "rs1801133",  "gene": "MTHFR",   "phenotype": "Hyperhomocysteinemia",
     "category": "disease", "risk_allele": "T", "expected_alt": "A"},   # minus strand
    {"rsid": "rs9939609",  "gene": "FTO",     "phenotype": "Obesity / BMI",
     "category": "disease", "risk_allele": "A", "expected_alt": "A"},
    {"rsid": "rs429358",   "gene": "APOE",    "phenotype": "Alzheimer's disease",
     "category": "disease", "risk_allele": "C", "expected_alt": "C"},
    {"rsid": "rs10455872", "gene": "LPA",     "phenotype": "Coronary artery disease",
     "category": "disease", "risk_allele": "G", "expected_alt": "G"},
    {"rsid": "rs3764261",  "gene": "CETP",    "phenotype": "Low HDL cholesterol",
     "category": "disease", "risk_allele": "A", "expected_alt": "A"},
    {"rsid": "rs7903146",  "gene": "TCF7L2",  "phenotype": "Type 2 diabetes",
     "category": "disease", "risk_allele": "T", "expected_alt": "T"},

    # pharmacogenomic
    {"rsid": "rs4149056",  "gene": "SLCO1B1", "phenotype": "Statin myopathy (*5)",
     "category": "pgx", "risk_allele": "C", "expected_alt": "C"},
    {"rsid": "rs1799853",  "gene": "CYP2C9",  "phenotype": "Warfarin dose (*2)",
     "category": "pgx", "risk_allele": "T", "expected_alt": "T"},
    {"rsid": "rs1142345",  "gene": "TPMT",    "phenotype": "Thiopurine toxicity (*3C)",
     "category": "pgx", "risk_allele": "C", "expected_alt": "C"},
    {"rsid": "rs3918290",  "gene": "DPYD",    "phenotype": "Fluoropyrimidine toxicity",
     "category": "pgx", "risk_allele": "A", "expected_alt": "T"},       # minus strand
    {"rsid": "rs1057910",  "gene": "CYP2C9",  "phenotype": "Warfarin dose (*3)",
     "category": "pgx", "risk_allele": "C", "expected_alt": "C"},
    {"rsid": "rs4244285",  "gene": "CYP2C19", "phenotype": "Clopidogrel (*2)",
     "category": "pgx", "risk_allele": "A", "expected_alt": "A"},
]

KEY_POPS = ["sas", "eas", "afr", "nfe", "amr", "mid"]

# reasons that mean "try again", as opposed to a real finding about the variant
NETWORK_FAILURE = "LOOKUP FAILED (network — retry needed)"


def backoff(attempt, response=None, label=""):
    """
    Wait before retrying. If the server sent Retry-After, obey it — that is the
    server stating how long its rate-limit window has left to run. Otherwise
    ramp 4, 8, 16, 32 seconds.
    """
    wait = BACKOFF_BASE * (2 ** attempt)
    if response is not None:
        hdr = response.headers.get("Retry-After")
        if hdr:
            try:
                wait = max(int(hdr), 1)
            except ValueError:
                pass
    print(f"  {label}waiting {wait}s before retry")
    time.sleep(wait)


def get_position(rsid, max_retries=4):
    """
    Ensembl rate-limits at roughly 15 requests per second and returns 429 when
    exceeded; it also returns intermittent 500s under load. Neither is evidence
    that the variant does not exist.

    Returns (chrom, pos, reason). reason is "ok" on success and describes the
    failure otherwise, so the caller can tell a network fault from an absence.
    """
    url = f"https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"

    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=30)

            if r.status_code == 200:
                for m in r.json().get("mappings", []):
                    if m.get("assembly_name", "").startswith("GRCh38"):
                        return m.get("seq_region_name"), m.get("start"), "ok"
                print(f"  no GRCh38 mapping returned for {rsid}")
                return None, None, "no GRCh38 mapping"

            if r.status_code == 404:
                print(f"  {rsid} not found in Ensembl")
                return None, None, "not in Ensembl"

            print(f"  Ensembl returned HTTP {r.status_code} (attempt {attempt + 1})")
            backoff(attempt, r, "Ensembl: ")
            continue

        except requests.exceptions.RequestException as e:
            print(f"  network error on attempt {attempt + 1}: {type(e).__name__}")
            backoff(attempt, None, "Ensembl: ")

    print(f"  GAVE UP after {max_retries} attempts. This is a network failure, "
          f"not an absence — do not treat {rsid} as missing.")
    return None, None, NETWORK_FAILURE


def gnomad_variant_id(rsid, expected_alt=None, max_retries=4):
    """
    An rsID identifies a position, not a substitution. At a multi-allelic site
    the search endpoint returns several variant IDs, and taking the first one
    silently selects whichever gnomAD happens to order first. For CYP2C19
    rs4244285 that is G>T, not the G>A that defines *2; for DPYD rs3918290 the
    first hit happens to be correct. The first hit is not reliably wrong, it is
    reliably unverified, which is why the expected ALT must be declared.
    """
    q = """
    query($rsid: String!, $ds: DatasetId!) {
      rsid_search: variant_search(query: $rsid, dataset: $ds) { variant_id }
    }"""

    hits = None
    for attempt in range(max_retries):
        try:
            r = requests.post("https://gnomad.broadinstitute.org/api",
                              json={"query": q, "variables": {"rsid": rsid, "ds": "gnomad_r4"}},
                              timeout=30)
            if r.status_code == 200:
                hits = [h["variant_id"]
                        for h in r.json().get("data", {}).get("rsid_search", [])]
                break
            print(f"  gnomAD search returned HTTP {r.status_code} (attempt {attempt + 1})")
            backoff(attempt, r, "gnomAD: ")
        except requests.exceptions.RequestException as e:
            print(f"  gnomAD search network error on attempt {attempt + 1}: "
                  f"{type(e).__name__}")
            backoff(attempt, None, "gnomAD: ")

    if hits is None:
        print(f"  GAVE UP searching gnomAD for {rsid} — network failure, not an absence.")
        return None, NETWORK_FAILURE

    if not hits:
        return None, "not indexed in gnomAD"

    if len(hits) > 1:
        print(f"  MULTI-ALLELIC: {len(hits)} variant IDs at this rsID: {hits}")

    if expected_alt:
        matches = [v for v in hits if v.split("-")[3] == expected_alt]
        if len(matches) == 1:
            if len(hits) > 1:
                print(f"  selected {matches[0]} on expected ALT '{expected_alt}'")
            return matches[0], "ok"
        if len(matches) > 1:
            print(f"  AMBIGUOUS: several hits carry ALT '{expected_alt}': {matches}")
            return None, "ambiguous ALT match"
        print(f"  NO HIT carries the expected ALT '{expected_alt}'. If this gene is "
              f"on the minus strand, expected_alt may need to be the complement of "
              f"the literature allele. Confirm against Ensembl VEP or dbSNP.")
        return None, f"no hit with expected ALT '{expected_alt}'"

    if len(hits) > 1:
        print("  No expected ALT declared, so the first hit cannot be trusted. "
              "Add expected_alt to GENE_PANEL and re-run.")
        return None, "multi-allelic, no expected ALT declared"

    return hits[0], "ok"


def gnomad_frequencies(variant_id, max_retries=4):
    q = """
    query($vid: String!, $ds: DatasetId!) {
      variant(variantId: $vid, dataset: $ds) {
        genome { af populations { id ac an } }
      }
    }"""
    for attempt in range(max_retries):
        try:
            r = requests.post("https://gnomad.broadinstitute.org/api",
                              json={"query": q, "variables": {"vid": variant_id, "ds": "gnomad_r4"}},
                              timeout=30)
            if r.status_code == 200:
                v = r.json().get("data", {}).get("variant")
                if v and v.get("genome"):
                    return v["genome"], "ok"
                return None, "no genome record in gnomAD"
            print(f"  gnomAD frequency query returned HTTP {r.status_code} "
                  f"(attempt {attempt + 1})")
            backoff(attempt, r, "gnomAD: ")
        except requests.exceptions.RequestException as e:
            print(f"  gnomAD frequency network error on attempt {attempt + 1}: "
                  f"{type(e).__name__}")
            backoff(attempt, None, "gnomAD: ")

    print(f"  GAVE UP fetching frequencies for {variant_id} — network failure.")
    return None, NETWORK_FAILURE


def india_record(chrom, pos):
    path = os.path.join(DATA_DIR, f"GI_9768_CBR-NIBMG_JointCall_AF_chr{chrom}.tsv")
    if not os.path.exists(path):
        print(f"  No file for chr{chrom} at {path}")
        print("   (chrX and chrY are not in this release; otherwise check DATA_DIR)")
        return None
    cols = ["CHROM", "POS", "ID", "REF", "ALT", "ALT_AF"]
    for chunk in pd.read_csv(path, sep="\t", header=None, names=cols, chunksize=500_000):
        hit = chunk[chunk["POS"] == pos]
        if len(hit):
            return hit.iloc[0].to_dict()
    return None


def check_alleles(india, gnomad_vid):
    """gnomAD IDs look like 10-112998590-C-T. Compare REF/ALT against GenomeIndia."""
    parts = gnomad_vid.split("-")
    g_ref, g_alt = parts[2], parts[3]
    i_ref, i_alt = india["REF"], india["ALT"]
    if (i_ref, i_alt) == (g_ref, g_alt):
        return "match", g_ref, g_alt
    if (i_ref, i_alt) == (g_alt, g_ref):
        return "flipped", g_ref, g_alt
    return "mismatch", g_ref, g_alt


def main():
    if not os.path.isdir(DATA_DIR):
        sys.exit(f"GenomeIndia data not found at {DATA_DIR}\n"
                 f"Download the summary statistics from "
                 f"https://ibdc.dbt.gov.in/genomeindia and extract the TSVs there.")

    rows, dropped = [], []

    for e in GENE_PANEL:
        print(f"\n=== {e['gene']} ({e['rsid']}) — {e['phenotype']} ===")

        try:
            chrom, pos, why = get_position(e["rsid"])
            if not chrom:
                dropped.append((e["gene"], e["rsid"], why))
                continue
            print(f"  chr{chrom}:{pos}")

            vid, why = gnomad_variant_id(e["rsid"], e.get("expected_alt"))
            if not vid:
                dropped.append((e["gene"], e["rsid"], why))
                continue
            print(f"  gnomAD ID: {vid}")

            india = india_record(chrom, pos)
            if india is None:
                dropped.append((e["gene"], e["rsid"], "absent from GenomeIndia release"))
                continue

            status, g_ref, g_alt = check_alleles(india, vid)
            print(f"  Alleles — India {india['REF']}>{india['ALT']}, "
                  f"gnomAD {g_ref}>{g_alt} => {status.upper()}")

            if status == "mismatch":
                print("  DISCARDING: different substitution at the same position; "
                      "frequencies are not comparable.")
                dropped.append((e["gene"], e["rsid"], "allele mismatch"))
                continue

            india_af = india["ALT_AF"]
            if status == "flipped":
                india_af = 1 - india_af
                print(f"  Orientation flipped, India AF corrected to {india_af:.4f}")

            g, why = gnomad_frequencies(vid)
            if not g:
                dropped.append((e["gene"], e["rsid"], why))
                continue

            pops = {}
            for p in g.get("populations", []):
                if p["id"] in KEY_POPS and p["an"] > 0:
                    pops[p["id"]] = round(p["ac"] / p["an"], 6)

            rows.append({
                "gene": e["gene"], "rsid": e["rsid"], "phenotype": e["phenotype"],
                "category": e["category"], "risk_allele": e["risk_allele"],
                "expected_alt": e["expected_alt"],
                "variant_id": vid, "chrom": chrom, "pos": pos,
                "allele_check": status,
                "india_ref": india["REF"], "india_alt": india["ALT"],
                "india_af": round(india_af, 6),
                "gnomad_global_af": round(g["af"], 6),
                "deviation": round(india_af - g["af"], 6),
                **{f"gnomad_{k}_af": v for k, v in pops.items()},
            })
            print(f"  India {india_af:.4f} | global {g['af']:.4f} "
                  f"| deviation {india_af - g['af']:+.4f}")

        finally:
            # applied whether the locus succeeded or failed, so a failure does
            # not accelerate the next request into a rate limit
            time.sleep(REQUEST_PAUSE)

    df = pd.DataFrame(rows)
    network_failures = [d for d in dropped if NETWORK_FAILURE in d[2]]

    print(f"\n{'='*70}")
    print(f"Validated {len(df)} of {len(GENE_PANEL)} loci")

    if dropped:
        print("\nDropped:")
        for gene, rsid, why in dropped:
            print(f"  {gene:<9} {rsid:<12} {why}")

    if network_failures:
        out = os.path.join(OUT_DIR, "full_panel_INCOMPLETE.csv")
        df.to_csv(out, index=False)
        print(f"\n{len(network_failures)} locus/loci failed on the network rather "
              f"than on the data.")
        print(f"Writing to {out} so the existing full_panel.csv is not overwritten.")
        print("Re-run when the connection is stable.")
    else:
        out = os.path.join(OUT_DIR, "full_panel.csv")
        df.to_csv(out, index=False)
        print(f"\nSaved to {out}")

    print(f"{'='*70}")
    if len(df):
        print(df[["gene", "rsid", "category", "india_af",
                  "gnomad_global_af", "deviation"]].to_string(index=False))


if __name__ == "__main__":
    main()