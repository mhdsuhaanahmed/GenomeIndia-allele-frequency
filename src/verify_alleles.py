import requests
import time

# rsID -> what we expect to find, and which allele our panel script matched
TARGETS = [
    {
        "rsid": "rs1057910",
        "gene": "CYP2C9",
        "star": "*3",
        "expect": "p.Ile359Leu  (amino_acids I/L at 359)",
        "matched_allele": "C",
    },
    {
        "rsid": "rs1799853",
        "gene": "CYP2C9",
        "star": "*2",
        "expect": "p.Arg144Cys  (amino_acids R/C at 144)",
        "matched_allele": "T",
    },
    {
        "rsid": "rs3918290",
        "gene": "DPYD",
        "star": "*2A",
        "expect": "splice donor variant c.1905+1G>A (no amino acid change)",
        "matched_allele": "T",
    },
    {
        "rsid": "rs4149056",
        "gene": "SLCO1B1",
        "star": "*5",
        "expect": "p.Val174Ala  (amino_acids V/A at 174)",
        "matched_allele": "C",
    },
    {
        "rsid": "rs1142345",
        "gene": "TPMT",
        "star": "*3C",
        "expect": "p.Tyr240Cys  (amino_acids Y/C at 240)",
        "matched_allele": "C",
    },
    {
        "rsid": "rs4244285",
        "gene": "CYP2C19",
        "star": "*2",
        "expect": "splice defect c.681G>A (synonymous / splice consequence)",
        "matched_allele": "A",
    },
]


def vep(rsid, retries=3):
    url = f"https://rest.ensembl.org/vep/human/id/{rsid}?content-type=application/json"
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=40)
            if r.status_code == 200:
                return r.json(), None
            return None, f"HTTP {r.status_code}"
        except requests.exceptions.RequestException as e:
            if attempt < retries - 1:
                time.sleep(3)
                continue
            return None, type(e).__name__
    return None, "failed"


def summarise(entry, gene):
    """
    Collapse the many transcript rows into one line per (allele, consequence).
    Ensembl reports every annotated transcript, which is why the raw output is long.
    """
    by_allele = {}
    for tc in entry.get("transcript_consequences", []):
        if gene and tc.get("gene_symbol") != gene:
            continue
        allele = tc.get("variant_allele")
        aa = tc.get("amino_acids")
        pos = tc.get("protein_start")
        terms = tuple(sorted(tc.get("consequence_terms", [])))
        key = (allele, aa, pos, terms)
        by_allele[key] = by_allele.get(key, 0) + 1
    return by_allele


if __name__ == "__main__":
    for t in TARGETS:
        print(f"\n{'='*70}")
        print(f"{t['gene']} {t['star']}  ({t['rsid']})")
        print(f"  expecting      : {t['expect']}")
        print(f"  panel matched  : alt allele '{t['matched_allele']}'")
        print(f"{'-'*70}")

        data, err = vep(t["rsid"])
        if err:
            print(f"  VEP FAILED: {err}")
            continue

        for entry in data:
            print(f"  alleles at site: {entry.get('allele_string')}")
            rows = summarise(entry, t["gene"])
            if not rows:
                print("  no transcript consequences returned for this gene")
                continue

            print(f"  {'allele':<8}{'aa':<10}{'pos':<8}{'consequence':<45}{'n'}")
            for (allele, aa, pos, terms), n in sorted(
                    rows.items(), key=lambda x: (str(x[0][0]), str(x[0][2]))):
                marker = "  <-- matched" if allele == t["matched_allele"] else ""
                print(f"  {str(allele):<8}{str(aa or '-'):<10}{str(pos or '-'):<8}"
                      f"{','.join(terms)[:44]:<45}{n}{marker}")

        time.sleep(1)

    print(f"\n{'='*70}")
    print("Check each 'matched' row against the expected consequence above it.")
    print("A mismatch means the panel picked the wrong alternate allele at that site.")