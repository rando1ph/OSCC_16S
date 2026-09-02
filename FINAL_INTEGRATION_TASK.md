Continue the OSCC 16S project from the CURRENT VALIDATED STATE.

CRITICAL: all three cohort-level pipelines are complete. Do NOT rerun raw preprocessing, downloading, DADA2, read merging, dereplication, UCHIME, OTU clustering, taxonomy classification, or genus collapsing unless a validation check proves an existing result is corrupt.

Validated cohort state:

1. PRJNA666746
- DADA2 cohort
- taxonomy completed
- genus-level abundance and paired tumor/adjacent statistics completed

2. PRJNA822685
- DADA2 cohort
- taxonomy completed
- genus-level abundance and paired tumor/adjacent statistics completed

3. PRJNA813034
- sensitivity cohort
- source-compatible 97% OTU workflow
- singleton OTUs removed
- representative sequences filtered to retained OTUs
- taxonomy assigned using Greengenes 13_8 97% reference with consensus-vsearch
- genus-level collapse completed
- 283 genera
- 20 complete tumor/control pairs
- 0 excluded samples
- paired_genus_stats.tsv has 284 lines including header
- columns are:
  genus, complete_pairs, median_tumor_relabund,
  median_normal_relabund, median_paired_clr_effect,
  wilcoxon_p, bh_fdr_q
- depth_status.tsv sample orientation bug has been repaired
- taxonomy_processing.py now supports Greengenes k__Bacteria taxonomy
- otu_sensitivity.sh has been repaired and is reproducible

Existing successful results, especially combined/two_cohort/, MUST be preserved.

Your task now is ONLY downstream three-cohort analysis.

Required work:

1. Validate the three existing genus-level cohort outputs and statistics.
2. Harmonize genus labels across cohorts.
   - PRJNA666746 and PRJNA822685 use the existing DADA2 taxonomy workflow.
   - PRJNA813034 uses Greengenes 13_8 and must be explicitly treated as a sensitivity cohort.
   - Do not interpret taxonomy-database naming differences as biological absence.
   - Document unmatched or database-specific genus labels.
3. Perform three-cohort genus-level integration.
   - Keep cohort-specific effect estimates.
   - Compare effect direction and magnitude across cohorts.
   - Prefer shared/harmonizable genera for primary cross-cohort conclusions.
   - Do not pool raw counts across cohorts.
4. Produce a cross-cohort consistency table.
5. Produce a defensible candidate ranking for OSCC intratumoral microbes.
   Ranking should prioritize:
   - repeated direction across cohorts
   - effect magnitude
   - statistical evidence
   - prevalence / usable abundance where available
   - robustness across the two DADA2 cohorts and the PRJNA813034 sensitivity cohort
6. Generate useful non-empty summary tables and plots.
7. Update the preliminary report.
   Clearly document:
   - two DADA2 cohorts versus one OTU sensitivity cohort
   - Greengenes 13_8 taxonomy for PRJNA813034
   - genus-level integration rationale
   - taxonomy database heterogeneity as a limitation
   - that this is candidate prioritization, not causal proof
8. Do NOT spend time on optional machine learning yet.
9. Validate every final output for non-empty content, sensible sample counts, and coherent directions before completion.

Do not declare completion simply because a script was launched.

Only after the three-cohort integration, candidate ranking, plots/tables, report update, and final validation are actually complete, print:

RECOVERY_COMPLETE

followed by a concise summary of the final results and output locations.
