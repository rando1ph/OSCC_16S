Act as the persistent autonomous supervisor for this OSCC 16S project.

Do not exit after launching another script. Remain resident until the analysis is complete or a genuinely unrecoverable scientific issue requires human judgment.

Current state:
- PRJNA666746: download, MD5, preprocessing, DADA2 and export complete. Taxonomy previously exited 1 after a long VSEARCH run.
- PRJNA822685: download, MD5, preprocessing, DADA2 and export complete. Its taxonomy is currently running in tmux session twotax. Do not interrupt a healthy VSEARCH process.
- SILVA 138.2 dereplicated reference reads and taxonomy artifacts are complete.
- posttax watcher exists.
- PRJNA813034 has 40 selected paired samples and strict FR/RF processed reads. DADA2 must not be retried because its discrete quality scores prevent error-model estimation. It must use an OTU-based sensitivity-analysis route later.

Your responsibilities:

1. Stay running and supervise the project.
2. While long healthy jobs are active, check status/logs/processes about every 2 minutes and wait.
3. Do not make the human manually poll routine progress.
4. On mechanical/software failures, inspect the exact error, patch deterministic scripts, validate outputs and resume from the nearest valid checkpoint.
5. Never rerun an expensive successful stage merely because a downstream stage failed.
6. For PRJNA666746 taxonomy, first check whether taxonomy.qza and taxonomy-search-results.qza already exist and validate. If they do, never rerun the long VSEARCH search; resume only export/filter/genus processing.
7. Let the currently healthy PRJNA822685 taxonomy finish without interruption.
8. Complete taxonomy, genus tables and paired genus statistics for PRJNA666746 and PRJNA822685.
9. Then produce and preserve a two-cohort intermediate result: cross-cohort effects, candidate ranking, figures and preliminary report.
10. After that, process PRJNA813034 as an OTU-based sensitivity cohort using the existing strict processed paired FASTQs. Do not retry DADA2 and do not invent artificial quality scores.
11. Use a supported QIIME2/VSEARCH paired-read, dereplication and 99% OTU clustering workflow, SILVA 138.2 taxonomy, genus collapse and the existing exact Tumor/Matched_Normal pairing.
12. Transparently exclude low-depth samples and matched partners when required.
13. Then rerun three-cohort genus integration, clearly marking PRJNA813034 as OTU-based sensitivity validation.
14. Produce final candidate_ranking.tsv, consistent_genera.tsv, cohort_effects.tsv, figures and preliminary_report.md.
15. Run leave-one-cohort-out logistic regression only if technically valid and without leakage.

Scientific rules that must not be changed:
- cohort inclusion
- Tumor/Matched_Normal labels
- patient pairing
- primers
- fixed seed 20260901
- strict FR/RF filtering
- CLR approach
- paired Wilcoxon
- BH-FDR
- cross-cohort direction-consistency logic

Do not lower QC thresholds just to force success.
Do not invent metadata.
Do not claim causality.

Engineering rules:
- do not use sudo
- do not git push
- do not git reset --hard
- do not git clean
- do not reboot or shutdown
- do not delete raw FASTQ
- do not kill healthy unrelated processes
- check disk and memory before expensive work
- use at most 8 CPU threads
- avoid simultaneous memory-heavy QIIME jobs
- maintain combined/supervisor.log and combined/supervisor_status.tsv
- preserve resumability

Do not exit until:
A. the complete three-cohort analysis and final outputs are finished and validated; or
B. a genuinely unresolved scientific ambiguity requires human judgment.

For A, write combined/SUPERVISOR_COMPLETE.
For B, write combined/NEEDS_HUMAN_REVIEW.md after completing every independent stage that can still run.

Start by inspecting the current twotax process and PRJNA666746 taxonomy outputs. Do not interrupt healthy computation.