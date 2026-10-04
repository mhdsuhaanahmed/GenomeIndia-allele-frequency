"""
verify_alleles.py

Independent confirmation that the alternate allele selected for each panel
locus is the one the clinical literature names. This is the check behind the
paper's methodological claim, and it is run after the panel is built rather
than as part of it, so that the selection and its verification are separate.

Two classes of locus need different treatment:

  protein-coding — the amino acid change identifies the allele unambiguously.
      CYP2C9 *3 is p.Ile359Leu; only one substitution at that position produces
      it, so a VEP lookup settles the matter.

  splice-site — every substitution at a canonical splice donor produces the
      same consequence term and no amino acid change, so VEP cannot distinguish
      them. DPYD *2A and CYP2C19 *2 fall here and are resolved instead through
      HGVS notation, which names the base change explicitly.

Exit status is 0 only if every locus verifies.
"""

import sys
import time
import requests

ENSEMBL = "https://rest.ensembl.org"

TARGETS = [
    {
        "rsid": "rs1057910", "gene": "CYP2C9", "star": "*3",
        "method": "protein",
        "expect_aa": "I/L", "expect_pos": 359,
        "expect_hgvs": None,
        "matched_allele": "C",
        "describe": "p.Ile359Leu",
    },
    {
        "rsid": "rs1799853", "gene": "CYP2C9", "star": "*2",
        "method": "protein",
        "expect_aa": "R/C", "expect_pos": 144,
        "expect_hgvs": None,
        "matched_allele": "T",
        "describe": "p.Arg144Cys",
    },
    {
        "rsid": "rs4149056", "gene": "SLCO1B1", "star": "*5",
        "method": "protein",
        "expect_aa": "V/A", "expect_pos": 174,
        "expect_hgvs": None,
        "matched_allele": "C",
        "describe": "p.Val174Ala",
    },
    {
        "rsid": "rs1142345", "gene": "TPMT", "star": "*3C",
        "method": "protein",
        "expect_aa": "Y/C", "expect_pos": 240,
        "expect_hgvs": None,
        "matched_allele": "C",
        "describe": "p.Tyr240Cys",
    },
    {
        "rsid": "rs3918290", "gene": "DPYD", "star": "*2A",
        "method": "hgvs",
        "expect_aa": None, "expect_pos": None,
        "expect_hgvs": "c.1905+1G>A",
        "matched_allele": "T",
        "describe": "splice donor c.1905+1G>A (minus strand: genomic C>T)",
    },
    {
        "rsid": "rs4244285", "gene": "CYP2C19", "star": "*2",
        "method": "hgvs",
        "expect_aa": None, "expect_pos": None,
        "expect_hgvs": "c.681G>A",
        "matched_allele": "A",
        "describe": "splice defect c.681G>A",
    },
]


def get(url, retries=4):
    """Ensembl returns 429 under rate limiting and intermittent 500s under load.
    Neither is a result; both are retried."""
    for attempt in range(retries):
        try:
            r = requests.get(url, timeout=40)
            if r.status_code == 200:
                return r.json(), None
            if r.status_code == 404:
                return None, "not found in Ensembl"
            wait = r.headers.get("Retry-After")
            wait = int(wait) if wait and wait.isdigit() else 4 * (2 ** attempt)
            print(f"    HTTP {r.status_code} (attempt {attempt + 1}); waiting {wait}s")
            time.sleep(wait)
        except requests.exceptions.RequestException as e:
            wait = 4 * (2 ** attempt)
            print(f"    network error (attempt {attempt + 1}): {type(e).__name__}; "
                  f"waiting {wait}s")
            time.sleep(wait)
    return None, "lookup failed after retries"


def vep(rsid):
    return get(f"{ENSEMBL}/vep/human/id/{rsid}"
               f"?content-type=application/json&hgvs=1")


def consequences_by_allele(entry, gene):
    """
    Collapse transcript rows to one record per alternate allele. Ensembl reports
    every annotated transcript; the transcript count is not informative, the
    distinct consequence per allele is.
    """
    out = {}
    for tc in entry.get("transcript_consequences", []):
        if gene and tc.get("gene_symbol") != gene:
            continue
        allele = tc.get("variant_allele")
        if allele is None:
            continue
        rec = out.setdefault(allele, {"aa": set(), "pos": set(),
                                      "terms": set(), "hgvsc": set()})
        if tc.get("amino_acids"):
            rec["aa"].add(tc["amino_acids"])
        if tc.get("protein_start"):
            rec["pos"].add(tc["protein_start"])
        rec["terms"].update(tc.get("consequence_terms", []))
        if tc.get("hgvsc"):
            rec["hgvsc"].add(tc["hgvsc"].split(":")[-1])
    return out


def verify_protein(t, by_allele):
    """The matched allele must be the one producing the expected amino acid
    change at the expected position, and no other allele may produce it."""
    rec = by_allele.get(t["matched_allele"])
    if rec is None:
        return False, (f"allele '{t['matched_allele']}' not annotated on "
                       f"{t['gene']} transcripts")

    aa_ok = t["expect_aa"] in rec["aa"]
    pos_ok = t["expect_pos"] in rec["pos"]

    if not aa_ok:
        return False, (f"expected amino acids {t['expect_aa']}, "
                       f"VEP reports {sorted(rec['aa']) or 'none'}")
    if not pos_ok:
        return False, (f"expected protein position {t['expect_pos']}, "
                       f"VEP reports {sorted(rec['pos']) or 'none'}")

    others = [a for a, r in by_allele.items()
              if a != t["matched_allele"] and t["expect_aa"] in r["aa"]]
    if others:
        return False, (f"allele(s) {others} also produce {t['expect_aa']} — "
                       f"the consequence does not identify the allele uniquely")

    return True, f"{t['expect_aa']} at {t['expect_pos']}, unique to this allele"


def verify_hgvs(t, by_allele):
    """
    Splice-site variants: every substitution at a canonical donor gives the same
    consequence term, so the protein check cannot discriminate. HGVS coding
    notation names the base change explicitly and does discriminate.
    """
    rec = by_allele.get(t["matched_allele"])
    if rec is None:
        return False, (f"allele '{t['matched_allele']}' not annotated on "
                       f"{t['gene']} transcripts")

    found = {h for h in rec["hgvsc"] if t["expect_hgvs"] in h}
    if not found:
        return False, (f"expected HGVS {t['expect_hgvs']}, VEP reports "
                       f"{sorted(rec['hgvsc'])[:4] or 'none'}")

    others = [a for a, r in by_allele.items()
              if a != t["matched_allele"]
              and any(t["expect_hgvs"] in h for h in r["hgvsc"])]
    if others:
        return False, f"allele(s) {others} also carry {t['expect_hgvs']}"

    return True, f"HGVS {t['expect_hgvs']} confirmed, unique to this allele"


def main():
    results = []

    for t in TARGETS:
        print(f"\n{'='*72}")
        print(f"{t['gene']} {t['star']}  ({t['rsid']})")
        print(f"  expecting     : {t['describe']}")
        print(f"  panel matched : alt allele '{t['matched_allele']}'")
        print(f"  method        : {t['method']}")
        print(f"{'-'*72}")

        data, err = vep(t["rsid"])
        if err:
            print(f"  LOOKUP FAILED: {err}")
            print("  This is a failure to verify, not a failed verification.")
            results.append((t, None, err))
            time.sleep(1)
            continue

        verdict, reason = None, "no transcript consequences returned"
        for entry in data:
            print(f"  alleles at site: {entry.get('allele_string')}")
            by_allele = consequences_by_allele(entry, t["gene"])
            if not by_allele:
                continue

            print(f"  {'allele':<9}{'aa':<9}{'pos':<7}{'consequence':<40}")
            for allele in sorted(by_allele):
                r = by_allele[allele]
                aa = "/".join(sorted(r["aa"])) or "-"
                pos = ",".join(str(p) for p in sorted(r["pos"])) or "-"
                terms = ",".join(sorted(r["terms"]))[:39]
                mark = "  <-- matched" if allele == t["matched_allele"] else ""
                print(f"  {allele:<9}{aa:<12}{pos[:28]:<30}{terms:<40}{mark}")
            check = verify_protein if t["method"] == "protein" else verify_hgvs
            verdict, reason = check(t, by_allele)
            break

        print(f"  {'VERIFIED' if verdict else 'FAILED'}: {reason}")
        results.append((t, verdict, reason))
        time.sleep(1)

    print(f"\n{'='*72}")
    print("SUMMARY")
    print(f"{'='*72}")

    passed = [r for r in results if r[1] is True]
    failed = [r for r in results if r[1] is False]
    errored = [r for r in results if r[1] is None]

    for t, verdict, reason in results:
        status = {True: "VERIFIED", False: "FAILED", None: "NOT CHECKED"}[verdict]
        print(f"  {status:<12} {t['gene']:<9} {t['star']:<5} {t['rsid']:<12} {reason}")

    print(f"\n{len(passed)} verified, {len(failed)} failed, "
          f"{len(errored)} not checked")

    if failed:
        print("\nA FAILED locus means the panel selected an allele that does not "
              "produce the expected consequence. Do not use the panel until it "
              "is resolved.")
    if errored:
        print("\nNOT CHECKED means Ensembl could not be reached. Re-run; this is "
              "not evidence either way.")

    sys.exit(0 if (not failed and not errored) else 1)


if __name__ == "__main__":
    main()