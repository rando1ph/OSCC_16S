# OSCC 16S Multi-cohort Analysis Log

## Cohorts selected

### PRJNA666746
- 100 runs
- 50 OSCC tumor samples
- 50 matched normal samples
- 50 complete patient pairs
- Illumina MiSeq
- V3-V4
- Paired-end 251 bp

### PRJNA822685
- 162 runs
- 81 tumor samples
- 81 adjacent normal samples
- 81 complete patient pairs
- Illumina MiSeq
- V3-V4
- Paired-end

### PRJNA813034
Original dataset:
- 60 tumor
- 20 adjacent normal
- 15 precancer

Main analysis subset:
- 20 tumor
- 20 matched adjacent normal
- 20 complete patient pairs
- T1: 10 pairs
- T4: 10 pairs
- Illumina HiSeq 2500
- Paired-end

Unpaired tumors and precancer samples are excluded from the main
Tumor-vs-Matched_Normal analysis.

## PRJNA666746 primer orientation investigation

Raw paired reads have mixed library orientation.

Pilot sample:
SRR12750669

R1 primer orientation:
- Forward: 42.86%
- Reverse: 43.01%
- Unclassified: 14.13%

R2 primer orientation:
- Forward: 42.12%
- Reverse: 44.91%
- Unclassified: 12.97%

Pair orientation:
- FR: 38.47%
- RF: 38.28%

Cutadapt pilot:
- FR retained: 359,628 pairs (38.8%)
- RF retained: 359,262 pairs (38.7%)

Preprocessing strategy:
1. Identify FR and RF orientations using the expected V3-V4 primers.
2. Require both paired reads to contain the expected primer orientation.
3. Trim primers.
4. Reorient RF pairs so all processed reads use the same biological direction.
5. Discard pairs whose primer orientation cannot be confidently identified.

Primers:
- Forward: CCTACGGGNGGCWGCAG
- Reverse: GACTACHVGGGTATCTAATCC

Cutadapt pilot parameters:
- error rate: 0.15
- no indels
- discard untrimmed pairs

## Current analysis design

Each cohort is processed independently.

Raw reads and ASVs will not be merged directly across studies.

Downstream cross-cohort comparison will use a shared taxonomic level,
primarily genus-level abundance, followed by assessment of consistent
Tumor-vs-Matched_Normal effects across cohorts.

## Final validated state (2026-09-02)

All three cohort-level genus pipelines are complete and validated:

- PRJNA666746: DADA2_ASV, SILVA 138.2, 50 complete pairs, 480 genera.
- PRJNA822685: DADA2_ASV, SILVA 138.2, 56 complete pairs, 239 genera.
- PRJNA813034: OTU_VSEARCH_97 sensitivity cohort, Greengenes 13_8,
  20 complete pairs, 283 genera (0 genera q<0.1; 60 nominal p<0.05).

Genus-level statistics were re-derived from underlying genus count tables,
metadata and depth status and matched the on-disk files to file precision
(see scripts/overnight/validate_cohorts.py, PASS for all three cohorts).

Three-cohort genus-level integration, harmonization, cross-cohort
consistency, candidate ranking, tables, plots and the updated report are in
combined/three_cohort/. The preserved two-cohort results remain in
combined/two_cohort/. The integration treats PRJNA813034 as a Greengenes
OTU sensitivity cohort; counts are never pooled across cohorts and
taxonomy-database naming differences are not treated as biological absence.

Key generation scripts:
- scripts/overnight/three_cohort_integration.py
- scripts/overnight/validate_cohorts.py
