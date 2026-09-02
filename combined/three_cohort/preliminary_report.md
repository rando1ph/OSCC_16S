# OSCC 16S three-cohort preliminary report (genus level)

Generated: 2026-09-02 14:39:16
Cohorts integrated: PRJNA666746, PRJNA822685, PRJNA813034

## Cohorts and pipelines

- PRJNA666746: **DADA2 ASV cohort** (SILVA 138.2), 50 complete tumor/adjacent pairs.
- PRJNA822685: **DADA2 ASV cohort** (SILVA 138.2), 56 complete tumor/adjacent pairs after depth/pairing QC of 81.
- PRJNA813034: **OTU sensitivity cohort** (97% OTU workflow; Greengenes 13_8 taxonomy), 20 complete tumor/adjacent pairs.
- Two DADA2 (SILVA) cohorts + one Greengenes 13_8 OTU sensitivity cohort. PRJNA813034 is treated as sensitivity validation only, not as an equivalent DADA2 replication.

## Methods

- Each cohort was processed independently. Raw/ASV/OTU counts were **never pooled** across cohorts; cohort-specific effect estimates are retained throughout.
- Genus-level CLR transform with a single consistent pseudocount of 0.5.
- Per-genus paired Wilcoxon signed-rank test on paired Tumor-Matched_Normal CLR differences; BH-FDR across genera within each cohort.
- Genus harmonization uses exact genus-name matching after stripping database rank prefixes. SILVA 138.2 and Greengenes 13_8 largely agree at genus rank, but naming differences exist (e.g., SILVA lineage-group labels such as `Burkholderia-Caballeronia-Paraburkholderia`, `Christensenellaceae_R-7_group`, `[Eubacterium]_yurii_group` have no literal Greengenes counterpart; Greengenes has `Burkholderia`, `g__1-68`, candidate/environmental labels not present in SILVA). These naming differences are **not** interpreted as biological absence.
- Cross-cohort conclusions prefer genera detectable in >=2 cohorts (shared labels); unmatched/database-specific labels are catalogued separately in `harmonization_unmatched.tsv`.
- Ranking score weights (documented in `candidate_ranking.tsv`): direction consistency across cohorts, statistical support (q<0.1 and p<0.05), effect magnitude, prevalence/usable abundance, and robustness across the two DADA2 cohorts plus the Greengenes sensitivity cohort. This is candidate prioritisation, not causal proof.

## Cohort summary

| cohort      | analysis_route   | taxonomy        | design      | database   |   samples_processed |   complete_pairs |   excluded_samples |   low_depth_excluded |   n_genera |   n_genera_nominal_p0.05 |   n_genera_sig_q0.1 |   mean_raw_pairs |   mean_final_retained_pairs |   mean_bacterial_reads_postqc |   low_bacterial_reads_samples |
|:------------|:-----------------|:----------------|:------------|:-----------|--------------------:|-----------------:|-------------------:|---------------------:|-----------:|-------------------------:|--------------------:|-----------------:|----------------------------:|------------------------------:|------------------------------:|
| PRJNA666746 | DADA2_ASV        | SILVA_138.2     | primary     | SILVA      |                 100 |               50 |                  0 |                    0 |        480 |                      372 |                 389 |         974730   |                     99757.6 |                       37635.5 |                           nan |
| PRJNA822685 | DADA2_ASV        | SILVA_138.2     | primary     | SILVA      |                 162 |               56 |                 50 |                   27 |        239 |                      151 |                 168 |          41933.2 |                     33825.1 |                        7233.4 |                           nan |
| PRJNA813034 | OTU_VSEARCH_97   | Greengenes_13_8 | sensitivity | Greengenes |                  40 |               20 |                  0 |                    0 |        283 |                       60 |                   0 |         554618   |                     87029.6 |                       61781.9 |                             2 |

## Genus-label sharing / harmonization

- Shared genus labels detected in >=2 cohorts: 259 (see `cross_cohort_consistency.tsv`).
- Genus labels found in all three cohorts: 99 biological genera (plus the 'Unclassified' pseudo-label).
- Database-specific (unmatched) labels: 382 (see `harmonization_unmatched.tsv`).

## Consistent genera

- Genera with the same effect direction in >=2 cohorts and q<0.1 in >=1 cohort: 147 (see `consistent_genera.tsv`).

| genus                 |   n_cohorts |   n_DADA2_SILVA |   in_sensitivity_GG | direction   |   n_sig_q0.1 |       min_q |   mean_prevalence |   mean_abs_effect | cohort_effects                                               |
|:----------------------|------------:|----------------:|--------------------:|:------------|-------------:|------------:|------------------:|------------------:|:-------------------------------------------------------------|
| g__Campylobacter      |           3 |               2 |                   1 | up          |            2 | 0.0154164   |            0.8342 |         0.469237  | 666746(SILVA*+0.993);822685(SILVA*+0.342);813034(GG+0.0722)  |
| g__Peptostreptococcus |           3 |               2 |                   1 | up          |            2 | 0.0180647   |            0.7719 |         0.737597  | 666746(SILVA*+0.225);822685(SILVA*+0.815);813034(GG+1.17)    |
| g__Catonella          |           3 |               2 |                   1 | up          |            2 | 0.013465    |            0.6785 |         0.542621  | 666746(SILVA*+0.666);822685(SILVA*+0.58);813034(GG+0.382)    |
| g__Actinomyces        |           3 |               2 |                   1 | down        |            2 | 0.00412384  |            0.5789 |         0.623546  | 666746(SILVA*-1.04);822685(SILVA*-0.111);813034(GG-0.723)    |
| g__Megasphaera        |           3 |               2 |                   1 | down        |            2 | 0.0455351   |            0.5567 |         0.509823  | 666746(SILVA*-0.931);822685(SILVA*-0.088);813034(GG-0.51)    |
| g__Filifactor         |           3 |               2 |                   1 | up          |            2 | 0.000390308 |            0.4736 |         0.359407  | 666746(SILVA*+0.25);822685(SILVA*+0.197);813034(GG+0.631)    |
| g__Acidovorax         |           3 |               2 |                   1 | up          |            2 | 0.0338832   |            0.2058 |         0.60328   | 666746(SILVA*+0.142);822685(SILVA*+0.0836);813034(GG+1.58)   |
| g__Azospirillum       |           3 |               2 |                   1 | up          |            2 | 0.00157225  |            0.1774 |         0.191815  | 666746(SILVA*+0.135);822685(SILVA*+0.0746);813034(GG+0.366)  |
| g__Deinococcus        |           3 |               2 |                   1 | up          |            2 | 0.01315     |            0.1637 |         0.175007  | 666746(SILVA*+0.216);822685(SILVA*+0.0836);813034(GG+0.225)  |
| g__Brevibacterium     |           3 |               2 |                   1 | up          |            2 | 0.0179546   |            0.1454 |         0.168204  | 666746(SILVA*+0.139);822685(SILVA*+0.0836);813034(GG+0.282)  |
| g__Morganella         |           3 |               2 |                   1 | up          |            2 | 0.00156521  |            0.1348 |         0.178663  | 666746(SILVA*+0.139);822685(SILVA*+0.0836);813034(GG+0.313)  |
| g__Desulfovibrio      |           3 |               2 |                   1 | up          |            2 | 0.00148434  |            0.1208 |         0.0892981 | 666746(SILVA*+0.152);822685(SILVA*+0.0746);813034(GG+0.0412) |
| g__Bilophila          |           3 |               2 |                   1 | up          |            2 | 0.000496365 |            0.1154 |         0.118061  | 666746(SILVA*+0.139);822685(SILVA*+0.0836);813034(GG+0.131)  |
| g__Anaerovorax        |           3 |               2 |                   1 | up          |            2 | 0.000391488 |            0.1137 |         0.124766  | 666746(SILVA*+0.152);822685(SILVA*+0.0746);813034(GG+0.148)  |
| g__Anaerococcus       |           3 |               2 |                   1 | up          |            2 | 0.0172614   |            0.0944 |         0.186963  | 666746(SILVA*+0.135);822685(SILVA*+0.0746);813034(GG+0.351)  |
| g__Mycobacterium      |           3 |               2 |                   1 | up          |            2 | 0.0057271   |            0.0871 |         0.0962479 | 666746(SILVA*+0.13);822685(SILVA*+0.0836);813034(GG+0.0746)  |
| g__Finegoldia         |           3 |               2 |                   1 | up          |            2 | 0.0388238   |            0.0773 |         0.154691  | 666746(SILVA*+0.0911);822685(SILVA*+0.0912);813034(GG+0.282) |
| g__Dietzia            |           3 |               2 |                   1 | up          |            2 | 0.000566597 |            0.0691 |         0.156978  | 666746(SILVA*+0.171);822685(SILVA*+0.0746);813034(GG+0.225)  |
| g__Lysinibacillus     |           3 |               2 |                   1 | up          |            2 | 0.00541194  |            0.0671 |         0.213118  | 666746(SILVA*+0.135);822685(SILVA*+0.0836);813034(GG+0.421)  |
| g__Aerococcus         |           3 |               2 |                   1 | up          |            2 | 0.0230813   |            0.0571 |         0.165309  | 666746(SILVA*+0.13);822685(SILVA*+0.0836);813034(GG+0.282)   |

## Top 25 candidate genera (three-cohort ranking)

|   rank | tier                             | genus                 |   n_cohorts |   n_DADA2_SILVA |   in_sensitivity_GG | direction   |   n_sig_q0.1 |   n_nominal_p0.05 |       min_q |   mean_prevalence |   mean_abs_effect |   mean_median_effect |   score |
|-------:|:---------------------------------|:----------------------|------------:|----------------:|--------------------:|:------------|-------------:|------------------:|------------:|------------------:|------------------:|---------------------:|--------:|
|      1 | A_all3_consistent                | g__Neisseria          |           3 |               2 |                   1 | down        |            1 |                 2 | 0.0114911   |            0.8802 |          1.28739  |            -1.28739  | 5.99793 |
|      2 | A_all3_consistent                | g__Veillonella        |           3 |               2 |                   1 | down        |            1 |                 2 | 2.04004e-05 |            0.8374 |          1.17036  |            -1.17036  | 5.51309 |
|      3 | A_all3_consistent                | g__Peptostreptococcus |           3 |               2 |                   1 | up          |            2 |                 1 | 0.0180647   |            0.7719 |          0.737597 |             0.737597 | 5.14257 |
|      4 | A_all3_consistent                | g__Rothia             |           3 |               2 |                   1 | down        |            1 |                 2 | 0.00410632  |            0.7214 |          1.27412  |            -1.27412  | 5.12908 |
|      5 | A_all3_consistent                | g__Campylobacter      |           3 |               2 |                   1 | up          |            2 |                 2 | 0.0154164   |            0.8342 |          0.469237 |             0.469237 | 4.74906 |
|      6 | A_all3_consistent                | g__Parvimonas         |           3 |               2 |                   1 | up          |            1 |                 2 | 0.0388238   |            0.7419 |          0.937794 |             0.937794 | 4.56178 |
|      7 | A_all3_consistent                | g__Streptococcus      |           3 |               2 |                   1 | down        |            1 |                 1 | 0.000109322 |            0.978  |          0.874726 |            -0.874726 | 4.39262 |
|      8 | A_all3_consistent                | g__Catonella          |           3 |               2 |                   1 | up          |            2 |                 2 | 0.013465    |            0.6785 |          0.542621 |             0.542621 | 4.23263 |
|      9 | A_all3_consistent                | g__Actinomyces        |           3 |               2 |                   1 | down        |            2 |                 2 | 0.00412384  |            0.5789 |          0.623546 |            -0.623546 | 3.9416  |
|     10 | A_all3_consistent                | g__Megasphaera        |           3 |               2 |                   1 | down        |            2 |                 1 | 0.0455351   |            0.5567 |          0.509823 |            -0.509823 | 3.6131  |
|     11 | A_all3_consistent                | g__Granulicatella     |           3 |               2 |                   1 | down        |            1 |                 1 | 9.18705e-06 |            0.7974 |          0.672887 |            -0.672887 | 3.40413 |
|     12 | A_all3_consistent                | g__Fusobacterium      |           3 |               2 |                   1 | up          |            1 |                 1 | 0.0093789   |            0.8185 |          0.611101 |             0.611101 | 3.36576 |
|     13 | A_all3_consistent                | g__Filifactor         |           3 |               2 |                   1 | up          |            2 |                 2 | 0.000390308 |            0.4736 |          0.359407 |             0.359407 | 2.97643 |
|     14 | A_all3_consistent                | g__Treponema          |           3 |               2 |                   1 | up          |            1 |                 1 | 0.0308664   |            0.6512 |          0.52447  |             0.52447  | 2.71229 |
|     15 | A_all3_consistent                | g__Corynebacterium    |           3 |               2 |                   1 | down        |            1 |                 1 | 0.0230813   |            0.6065 |          0.478198 |            -0.478198 | 2.51355 |
|     16 | A_all3_consistent                | g__Acidovorax         |           3 |               2 |                   1 | up          |            2 |                 3 | 0.0338832   |            0.2058 |          0.60328  |             0.60328  | 2.50248 |
|     17 | A_all3_consistent                | g__Oribacterium       |           3 |               2 |                   1 | down        |            1 |                 1 | 0.0388238   |            0.6945 |          0.275595 |            -0.275595 | 2.46181 |
|     18 | A_all3_consistent                | g__Peptococcus        |           3 |               2 |                   1 | up          |            1 |                 1 | 0.0180647   |            0.5116 |          0.335763 |             0.335763 | 2.05932 |
|     19 | C_DADA2_consistent_GGnotdetected | g__Bergeyella         |           2 |               2 |                   0 | up          |            1 |                 1 | 0.0179546   |            0.5172 |          0.769624 |             0.769624 | 1.93607 |
|     20 | A_all3_consistent                | g__Azospirillum       |           3 |               2 |                   1 | up          |            2 |                 1 | 0.00157225  |            0.1774 |          0.191815 |             0.191815 | 1.58399 |
|     21 | A_all3_consistent                | g__Deinococcus        |           3 |               2 |                   1 | up          |            2 |                 2 | 0.01315     |            0.1637 |          0.175007 |             0.175007 | 1.51632 |
|     22 | C_DADA2_consistent_GGnotdetected | g__Lentimicrobium     |           2 |               2 |                   0 | up          |            2 |                 2 | 0.0388238   |            0.2948 |          0.199935 |             0.199935 | 1.50358 |
|     23 | A_all3_consistent                | g__Brevibacterium     |           3 |               2 |                   1 | up          |            2 |                 2 | 0.0179546   |            0.1454 |          0.168204 |             0.168204 | 1.44271 |
|     24 | A_all3_consistent                | g__Morganella         |           3 |               2 |                   1 | up          |            2 |                 2 | 0.00156521  |            0.1348 |          0.178663 |             0.178663 | 1.41381 |
|     25 | A_all3_consistent                | g__Gemella            |           3 |               2 |                   1 | up          |            0 |                 0 | 0.187167    |            0.9045 |          0.171493 |             0.171493 | 1.40723 |

Tiers: A = all 3 cohorts, consistent; B = all 3 cohorts, discordant; C = both DADA2 cohorts consistent, Greengenes label not detected; D = both DADA2 cohorts discordant, GG not detected; E = one SILVA + GG consistent; F = two cohorts consistent; Z = discordant/partial.

## Outputs

- `combined/three_cohort/pipeline_status.tsv` - pipeline state
- `combined/three_cohort/cohort_summary.tsv` - cohort summary
- `combined/three_cohort/cohort_effects.tsv` - long-format cohort effect estimates
- `combined/three_cohort/genus_harmonization.tsv` - union genus catalogue
- `combined/three_cohort/harmonization_unmatched.tsv` - database-specific labels
- `combined/three_cohort/cross_cohort_consistency.tsv` - shared-genus consistency
- `combined/three_cohort/consistent_genera.tsv` - consistent, q<0.1 in >=1
- `combined/three_cohort/candidate_ranking.tsv` - ranked candidates
- `combined/three_cohort/effect_heatmap.png`, `top_candidates.png`, `direction_dotplot.png`, `pipeline_retention.png`

## Limitations and caveats

- PRJNA813034 is an OTU sensitivity cohort using a different denoising/OTU approach (97% OTU vs DADA2 ASV) and a different taxonomy database (Greengenes 13_8 vs SILVA 138.2); effect signs there are corroborating, not equivalent.
- Taxonomy-database heterogeneity: a genus label absent from one database/cohort may reflect naming conventions, resolution, or primer/platform differences rather than biological absence.
- PRJNA813034 has 20 pairs only; **no genus reaches q<0.1** in that cohort (60 genera are nominal p<0.05). Its role is direction/magnitude corroboration.
- PRJNA666746 yields a very large number of q<0.1 genera (most very rare taxa); the CLR pseudocount can make sparse low-abundance taxa appear testable. Prevalence and magnitude are therefore reported for every candidate.
- In PRJNA813034 two tumor samples (SRR18236732, SRR18236733, T4_P04 and T4_P03) had essentially all reads unassignable at kingdom level by Greengenes (~123 and ~253 bacterial-mapped reads of ~35k); they remain in the 20 pairs as validated but add noise.
- Cross-study, cross-primer, cross-platform comparisons at genus level are exploratory. Associations are correlational; this analysis provides candidate prioritisation, not causal proof.
