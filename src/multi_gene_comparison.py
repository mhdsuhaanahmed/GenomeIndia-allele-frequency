import requests
import time
import pandas as pd
import os

OUT_DIR = r"D:\GENOMEINDIA\outputs\multi_gene"
os.makedirs(OUT_DIR, exist_ok=True)

GENE_PANEL = [
    # already validated
    {"rsid": "rs7903146",  "gene": "TCF7L2",  "disease": "Type 2 diabetes",          "risk_allele": "T"},
    {"rsid": "rs429358",   "gene": "APOE",    "disease": "Alzheimer's disease",      "risk_allele": "C"},
    {"rsid": "rs10455872", "gene": "LPA",     "disease": "Coronary artery disease",  "risk_allele": "G"},
    # new
    {"rsid": "rs1801133",  "gene": "MTHFR",   "disease": "Hyperhomocysteinemia",     "risk_allele": "T"},
    {"rsid": "rs9939609",  "gene": "FTO",     "disease": "Obesity / BMI",            "risk_allele": "A"},
    {"rsid": "rs3764261",  "gene": "CETP",    "disease": "Low HDL cholesterol",      "risk_allele": "A"},
    {"rsid": "rs4244285",  "gene": "CYP2C19", "disease": "Clopidogrel non-response", "risk_allele": "A"},
    {"rsid": "rs671",      "gene": "ALDH2",   "disease": "Alcohol flush / ESCC",     "risk_allele": "A"},
]

KEY_POPS = ["sas", "eas", "afr", "nfe", "amr", "mid"]


def get_position(rsid):
    url = f"https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"
    try:
        r = requests.get(url, timeout=20)
        if r.status_code == 200:
            for m in r.json().get("mappings", []):
                if m.get("assembly_name", "").startswith("GRCh38"):
                    return m.get("seq_region_name"), m.get("start")
    except Exception as e:
        print(f"  Ensembl lookup failed: {e}")
    return None, None


def gnomad_variant_id(rsid):
    q = """
    query($rsid: String!, $ds: DatasetId!) {
      rsid_search: variant_search(query: $rsid, dataset: $ds) { variant_id }
    }"""
    try:
        r = requests.post("https://gnomad.broadinstitute.org/api",
                          json={"query": q, "variables": {"rsid": rsid, "ds": "gnomad_r4"}},
                          timeout=25)
        if r.status_code == 200:
            hits = r.json().get("data", {}).get("rsid_search", [])
            if hits:
                return hits[0]["variant_id"]
    except Exception as e:
        print(f"  gnomAD search failed: {e}")
    return None


def gnomad_frequencies(variant_id):
    q = """
    query($vid: String!, $ds: DatasetId!) {
      variant(variantId: $vid, dataset: $ds) {
        genome { af populations { id ac an } }
      }
    }"""
    try:
        r = requests.post("https://gnomad.broadinstitute.org/api",
                          json={"query": q, "variables": {"vid": variant_id, "ds": "gnomad_r4"}},
                          timeout=25)
        if r.status_code == 200:
            v = r.json().get("data", {}).get("variant")
            if v and v.get("genome"):
                return v["genome"]
    except Exception as e:
        print(f"  gnomAD frequency fetch failed: {e}")
    return None


def india_record(chrom, pos):
    path = rf"D:\GENOMEINDIA\data\9768GI_SummaryStats\GI_9768_CBR-NIBMG_JointCall_AF_chr{chrom}.tsv"
    if not os.path.exists(path):
        print(f"  No file for chr{chrom} (chrX/chrY are not in this release)")
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


if __name__ == "__main__":
    rows = []

    for e in GENE_PANEL:
        print(f"\n=== {e['gene']} ({e['rsid']}) — {e['disease']} ===")

        chrom, pos = get_position(e["rsid"])
        if not chrom:
            print("  Could not resolve position. Skipping.")
            continue
        print(f"  chr{chrom}:{pos}")

        vid = gnomad_variant_id(e["rsid"])
        if not vid:
            print("  Not indexed in gnomAD. Skipping.")
            continue
        print(f"  gnomAD ID: {vid}")

        india = india_record(chrom, pos)
        if india is None:
            print("  Not found in GenomeIndia. Skipping.")
            continue

        status, g_ref, g_alt = check_alleles(india, vid)
        print(f"  Alleles — India {india['REF']}>{india['ALT']}, gnomAD {g_ref}>{g_alt} => {status.upper()}")

        if status == "mismatch":
            print("  DISCARDING: different substitution at same position, frequencies are not comparable.")
            continue

        india_af = india["ALT_AF"]
        if status == "flipped":
            india_af = 1 - india_af
            print(f"  Orientation flipped, India AF corrected to {india_af:.4f}")

        g = gnomad_frequencies(vid)
        pops = {}
        if g:
            for p in g.get("populations", []):
                if p["id"] in KEY_POPS and p["an"] > 0:
                    pops[p["id"]] = round(p["ac"] / p["an"], 4)

        rows.append({
            "gene": e["gene"], "rsid": e["rsid"], "disease": e["disease"],
            "risk_allele": e["risk_allele"], "chrom": chrom, "pos": pos,
            "allele_check": status,
            "india_ref": india["REF"], "india_alt": india["ALT"],
            "india_af": round(india_af, 6),
            "gnomad_global_af": round(g["af"], 6) if g else None,
            **{f"gnomad_{k}_af": v for k, v in pops.items()},
        })
        print(f"  India {india_af:.4f} | global {g['af']:.4f}" if g else "  gnomAD frequencies unavailable")

        time.sleep(1)

    df = pd.DataFrame(rows)
    out = os.path.join(OUT_DIR, "multi_gene_comparison_extended.csv")
    df.to_csv(out, index=False)
    print(f"\n\nSaved {len(df)} validated loci to {out}\n")
    print(df.to_string())