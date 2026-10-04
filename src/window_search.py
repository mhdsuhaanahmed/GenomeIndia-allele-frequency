import requests, json

# VEP endpoint returns protein consequences; the plain variation endpoint does not
r = requests.get(
    "https://rest.ensembl.org/vep/human/id/rs1057910?content-type=application/json",
    timeout=30)

if r.status_code != 200:
    print("Status:", r.status_code, r.text[:400])
else:
    data = r.json()
    for entry in data:
        print(f"\nInput: {entry.get('input')}  alleles: {entry.get('allele_string')}")
        for tc in entry.get("transcript_consequences", []):
            aa = tc.get("amino_acids")
            pos = tc.get("protein_start")
            if aa:  # only rows with a protein change
                print(f"  allele={tc.get('variant_allele')}  "
                      f"{aa} at {pos}  gene={tc.get('gene_symbol')}  "
                      f"{tc.get('consequence_terms')}")