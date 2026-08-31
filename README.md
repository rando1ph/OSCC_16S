# OSCC_16S

Multi-cohort analysis of tissue-associated microbiota in oral squamous
cell carcinoma (OSCC).

## Main cohorts

- PRJNA666746: 50 paired tumor/normal patients
- PRJNA822685: 81 paired tumor/adjacent-normal patients
- PRJNA813034: 20 strictly paired tumor/adjacent-normal patients used in the main analysis

Total main-analysis design:
- 151 paired patients
- 302 biological samples

## Repository contents

- `scripts/`: analysis and preprocessing scripts
- `docs/`: analysis decisions and workflow notes
- `*/metadata/`: cohort metadata, manifests and accession information
- `*/results/`: derived analysis outputs
- `combined/`: cross-cohort results

Raw and processed FASTQ files are intentionally excluded from Git.
