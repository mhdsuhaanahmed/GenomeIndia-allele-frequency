# Allele-frequency comparison of thirteen clinically annotated loci between GenomeIndia and global reference populations

Analysis code and derived data for a preprint comparing allele frequencies at
thirteen disease-susceptibility and pharmacogenomic loci between the GenomeIndia
pooled release (9,768 individuals) and six gnomAD v4 reference populations.

## Summary of results

Thirteen loci against six reference populations gives 78 comparisons, of which 60
are significant after Benjamini-Hochberg FDR correction.

No disease-susceptibility locus in the panel is elevated in India against a
majority of reference populations. Two of six pharmacogenomic loci are:
CYP2C19 rs4244285 against all six, and CYP2C9\*3 rs1057910 against four of six.
A Fisher exact test of the category hypothesis gives p = 0.4545, so the pattern
is not supported as a category effect at this panel size and is reported as a
locus-level observation.

Three independent validations support the comparison:

- ALDH2 rs671, included as a positive control, is near-absent in India
  (0.00005) against 0.225 in East Asians
- CYP2C9\*2 rs1799853 reproduces its documented near-absence in East Asian
  populations, odds ratio 106.8 (95% CI 26.7 to 428.0)
- Panel-wide Hudson FST between GenomeIndia and gnomAD South Asian is 0.00182,
  against 0.032 to 0.076 for all other reference populations

## Data sources

Neither source dataset is redistributed here. Both are publicly available.

**GenomeIndia pooled allele frequencies.** Indian Biological Data Centre,
https://ibdc.dbt.gov.in/genomeindia. Download the per-chromosome summary
statistics archive and extract the TSVs into `data/`. No registration or
data use agreement is required.

**gnomAD v4.** Fetched at runtime through the public GraphQL API at
https://gnomad.broadinstitute.org/api. No key required. The API rate-limits
under load; the scripts back off and retry, and cache responses so reruns are
fast.

## Requirements

```bash
pip install -r requirements.txt
```

Python 3.11 or later. The scripts use pandas, numpy, scipy, statsmodels,
matplotlib and requests.

## Running

```bash
python src/multi_gene_comparison.py     # locate and validate each locus
python src/statistical_analysis_v2.py   # Fisher exact, FDR, odds ratios, FST
python src/sensitivity_analysis.py      # allele-number sensitivity check
python src/plot_panel_v2.py             # the five paper figures
```

`multi_gene_comparison.py` must run first; the others read its output. Paths in
all scripts are relative to the repository root.

Three further scripts are not part of the pipeline. `verify_alleles.py` confirms
each selected allele against Ensembl VEP independently of the code that selected
it, and exits non-zero if any locus fails. `window_search.py` is the VKORC1
absence check. `query_gnomad_variant.py` dumps raw VEP consequences for a single
rsID, which is useful when a new locus misbehaves.

## Outputs

| File | Contents |
|---|---|
| `outputs/multi_gene/full_panel.csv` | 13 validated loci, frequencies across 7 groups |
| `outputs/stats/statistical_comparison_v2.csv` | 78 comparisons: odds ratios, CIs, q-values, per-locus FST |
| `outputs/stats/sensitivity_genotype_rate.csv` | baseline against upper-bound allele number |
| `outputs/stats/sensitivity_genotype_rate_per_locus.csv` | per-locus elevation counts at both bounds |
| `figures/panel_bars_light.png` | per-locus frequencies across seven groups |
| `figures/panel_heatmap_light.png` | all loci by population, India column outlined |
| `figures/panel_deviation_light.png` | deviation from global pooled frequency, by category |
| `figures/panel_forest_light.png` | odds ratios with 95% CIs, all 78 comparisons |
| `figures/panel_elevation_light.png` | populations each locus is significantly elevated against |

## A note on allele identity

Allele identity at every locus was confirmed against Ensembl VEP or dbSNP before
any frequency was used. This is not optional and it is the main methodological
point of the study.

An rsID identifies a genomic position, not a substitution. Four loci in this
panel sit at multi-allelic sites, where the gnomAD rsID search endpoint returns
several variant identifiers. The first result is not reliably wrong — it is
reliably unverified. For CYP2C19 rs4244285 the first hit is G>T rather than the
G>A that defines \*2; for DPYD rs3918290 the first hit happens to be correct.
Nothing distinguishes the two cases without checking. SLC30A8 rs13266634 was
discarded entirely after the two databases were found to describe different
substitutions at the same coordinate.

A second and separate problem is strand. The allele named in the clinical
literature and the alternate allele in the genomic representation are
complements on minus-strand genes: MTHFR C677T is G>A genomically, and DPYD
c.1905+1G>A is C>T. Two of the thirteen loci are affected, and careful rsID
handling does not catch it — the gap is between clinical naming and genomic
representation. The panel therefore records the literature allele and the
expected genomic ALT as separate fields.

Resolution method depends on the locus. Missense variants are settled by protein
consequence through Ensembl VEP, since only one substitution at a given position
produces a given amino acid change. Splice variants cannot be settled that way:
every substitution at a canonical splice donor returns the same consequence term
and no amino acid change. DPYD \*2A and CYP2C19 \*2 are resolved instead through
HGVS coding notation, which names the base change explicitly.

## Known limitation in the source data

VKORC1 rs9923231, one half of the standard warfarin dosing pair, is absent from
the public GenomeIndia release. This is not a coverage problem: 48 variants
appear within 500 bp of the position, including one 6 bp away. The chromosome 16
frequency distribution is unremarkable, with 188,956 variants above 5 percent,
so it is not a symptom of aggressive filtering. The exclusion appears isolated
and the reason is not determinable from the summary release. Anyone planning
pharmacogenomic work on this resource should check locus by locus rather than
assuming coverage.

## Citation

Preprint DOI to follow.

## License

MIT. See LICENSE.
