# OSCC 16S unattended overnight pipeline

Implement and launch a deterministic unattended preliminary analysis pipeline for this existing project.

Do not redesign the science. Inspect the existing metadata, manifests, scripts, README.md and docs/analysis_log.md first.

Hard rules:
- Do not change cohort inclusion, Tumor/Matched_Normal labels, or patient pairing.
- Do not delete raw FASTQ or overwrite original metadata.
- Do not use sudo, git push, git reset --hard, git clean, shutdown, reboot, pkill or killall.
- Never read PRJNA666746/logs/fastq_download.log because it is a multi-GB debug log.
- One cohort failing must not stop the others.
- All stages must be logged and resumable.
- The final pipeline must run independently of OpenCode after launch.

Environment:
- QIIME2 conda env: ~/miniconda3/envs/rachis-qiime2-2026.7
- cutadapt: ~/.local/bin/cutadapt
- FastQC, MultiQC, aria2c, pigz, Python3, tmux are available.
- Server: 8 vCPU, ~15 GiB RAM.
- Do not run multiple memory-heavy DADA2 jobs simultaneously.

Cohorts:

PRJNA666746:
- 50 paired patients, 100 samples.
- Download and MD5 already complete. logs/md5_status.txt = PASS.
- V3-V4 primers:
  Forward CCTACGGGNGGCWGCAG
  Reverse GACTACHVGGGTATCTAATCC
- Mixed orientation confirmed experimentally:
  FR ~38.8%, RF ~38.7%.
- Recognize both FR and RF, trim primers with cutadapt error rate 0.15 and no indels, discard untrimmed pairs, and normalize RF so final R1 is biological forward and R2 reverse.

PRJNA822685:
- 81 paired patients, 162 samples.
- fastq_urls.txt contains 324 FASTQ URLs.
- V3-V4 MiSeq paired-end.
- primers:
  Forward CCTACGGGNGGCWGCAG
  Reverse GGACTACNVGGGTWTCTAAT
- Automatically download, MD5 verify, inspect representative FASTQs for primer orientation, and normalize FR/RF if mixed.
- If <70% of sampled pairs have confidently recognized valid opposite primers, stop only this cohort and report failure.

PRJNA813034:
- Use ONLY the existing 40 samples in analysis_metadata.tsv: 20 strict Tumor/Matched_Normal pairs.
- fastq_urls.txt contains 80 FASTQ URLs.
- V3-V4 HiSeq 2500.
- primers:
  Forward CCTACGGGNBGCASCAG
  Reverse GACTACNVGGGTATCTAATCC
- Inspect/normalize primer orientation as above.
- After primer removal, trim 10 bases from the 3-prime end of both reads, matching the documented study workflow.
- If <70% valid primer orientation, stop only this cohort.

Preliminary fixed-depth strategy:
- fixed seed 20260901
- before expensive primer processing cap raw paired reads at 150000 pairs/sample
- preserve R1/R2 pairing
- after valid primer trimming/orientation normalization cap at 100000 valid pairs/sample
- record raw count, pretrim count, valid primer count and final retained count for every sample

Downloads:
- Start PRJNA822685 and PRJNA813034 downloads immediately.
- Use 2 aria2 concurrent downloads per cohort.
- Continue partial downloads.
- After download require expected FASTQ count, zero .aria2 files, and MD5 PASS.
- On MD5 failure redownload only failed files, maximum 3 attempts.

QIIME2:
- Process cohorts independently.
- Never merge raw reads or ASVs across studies.
- Import cleaned normalized paired FASTQs.
- Run q2-dada2 denoise-paired.
- Primers are already removed, trim-left=0.
- Do not invent arbitrary extra truncation.
- Produce table.qza, rep-seqs.qza, stats.qza and exported TSV/FASTA outputs.
- Record input, filtered, denoised, merged and non-chimeric reads.
- Samples with <2000 non-chimeric reads are LOW_DEPTH.
- If one member of a pair is excluded, remove the entire pair from downstream paired statistics and record the reason.

Taxonomy:
- Use a QIIME2-2026.7-compatible SILVA 138.2 taxonomy route.
- Prefer a current compatible SILVA classifier.
- If pretrained classifier acquisition is not reliable, use an appropriate supported QIIME2/RESCRIPt SILVA route.
- Never use an incompatible classifier.
- Remove chloroplast, mitochondria and inappropriate non-bacterial assignments.
- Collapse to genus level and export genus count and relative-abundance tables.
- Taxonomy failure must preserve successful DADA2 outputs.

Statistics:
- Use exact patient pairing from analysis_metadata.tsv.
- Genus-level paired exploratory analysis.
- CLR transform with one documented small pseudocount used consistently.
- For every genus report complete pairs, median Tumor abundance, median Matched_Normal abundance, median paired CLR Tumor-Normal effect, paired Wilcoxon p and BH-FDR q.
- Write paired_genus_stats.tsv for each cohort.

Cross-cohort:
If >=2 cohorts succeed, create in combined/:
- pipeline_status.tsv
- cohort_summary.tsv
- cohort_effects.tsv
- consistent_genera.tsv
- candidate_ranking.tsv
- preliminary_report.md
- effect_heatmap.png
- top_candidates.png
- pipeline_retention.png

Rank genera primarily by cross-cohort effect-direction consistency, statistical support and meaningful effect magnitude. Do not claim causality.

If all three genus tables succeed, additionally run exploratory leave-one-cohort-out regularized logistic regression:
- CLR genus features
- Tumor vs Matched_Normal target
- train on two cohorts, test on the third
- standardize on training data only
- no leakage
- report held-out ROC AUC and coefficient stability
- ML failure must not invalidate paired statistics.

Engineering:
- Create scripts under scripts/overnight/
- Create scripts/overnight_master.sh
- Use stage checkpoints and idempotent/resumable logic.
- Use a lock to prevent duplicate master instances.
- Write logs.
- Create combined/pipeline_status.tsv and combined/overnight_master.log.
- Fail closed on integrity problems.
- Do not wait for the overnight analysis to finish.

Execution priority:
1. inspect project quickly
2. verify PRJNA666746 MD5 PASS
3. immediately start downloads for the other two cohorts
4. immediately start PRJNA666746 preprocessing
5. build/test remaining scripts while those tasks run
6. launch deterministic master pipeline in tmux session named overnight
7. verify overnight exists and is running
8. write combined/automation_setup_summary.md
9. print AUTOMATION_LAUNCHED after successful launch
10. exit
