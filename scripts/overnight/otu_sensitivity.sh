#!/usr/bin/env bash
# PRJNA813034 OTU-based sensitivity analysis.
#
# DADA2 cannot be used for this cohort: its discrete quality scores prevent
# error-model estimation. This route uses a supported QIIME2/VSEARCH
# paired-read merge, dereplication, chimera filtering and 97% de novo OTU
# clustering workflow, then collapses to genus and computes the same paired
# statistics as the DADA2 cohorts.
#
# Parameters are taken from the source study methodology:
#   Pratap Singh R, et al. Intratumoral Microbiota Changes with Tumor Stage
#   and Influences the Immune Signature of Oral Squamous Cell Carcinoma.
#   Microbiol Spectr 2023;11(4):e04596-22 (PRJNA813034).
# The study merged paired reads (Fastq-join, QIIME 1.9.1) after trimming
# 10 bases from the 3' ends of both reads, removed chimeras, clustered into
# OTUs at 97% sequence similarity with USEARCH, and removed OTUs with only
# one read. VSEARCH (merge-pairs) is the modern equivalent of the QIIME
# 1.9.1 Fastq-join/UCLUST stack used by the study.
#
# Taxonomy: the study used UCLUST against Greengenes v13_8. To keep
# cross-cohort integration at genus level comparable with the DADA2 cohorts,
# classification here uses consensus-VSEARCH against the same SILVA 138.2
# reference used for PRJNA666746 and PRJNA822685 (documented deviation).
#
# Resumable: stage markers per step; never reruns a completed step.
# usage: bash scripts/overnight/otu_sensitivity.sh
set -euo pipefail
source "$(dirname "$0")/common.sh"

COHORT="PRJNA813034"
BASE="$(cohort_dir "$COHORT")"
RES="$BASE/results"
LOGS="$BASE/logs"
LOCK="$BASE/.otu.lock"
MARKER="$BASE/.otu_done"
REF_READS="$ROOT/data/silva/silva-138.2-derep.qza"
REF_TAX="$ROOT/data/silva/silva-138.2-derep-tax.qza"

export PATH="$QIIME_ENV/bin:$PATH"

mkdir -p "$RES" "$LOGS" "$RES/export"

acquire_lock "$LOCK" || { echo "cannot lock otu_sensitivity for $COHORT"; exit 1; }

if stage_done "$MARKER"; then
    echo "otu sensitivity already complete for $COHORT"
    exit 0
fi

if [[ ! -f "$RES/demux.qza" ]]; then
    echo "$(NOW) [$COHORT] demux.qza missing (must import processed paired FASTQs first)" | tee -a "$LOGS/otu.log"
    exit 1
fi

# ---- 1. merge paired reads (VSEARCH merge-pairs, no error-model needed) ----
if [[ ! -s "$RES/joined.qza" ]]; then
    echo "$(NOW) [$COHORT] vsearch merge-pairs" | tee -a "$LOGS/otu.log"
nano ~/OSCC_16S/scripts/overnight/otu_sensitivity.sh    "$QIIME_ENV/bin/qiime" vsearch merge-pairs \
        --i-demultiplexed-seqs "$RES/demux.qza" \
        --p-threads 4 \
        --o-merged-sequences "$RES/joined.qza" \
        --o-unmerged-sequences "$RES/unjoined.qza" \
        >> "$LOGS/otu.log" 2>&1
fi

# ---- 2. dereplicate joined sequences ----
if [[ ! -s "$RES/derep-table.qza" ]]; then
    echo "$(NOW) [$COHORT] vsearch dereplicate-sequences" | tee -a "$LOGS/otu.log"
    "$QIIME_ENV/bin/qiime" vsearch dereplicate-sequences \
        --i-sequences "$RES/joined.qza" \
        --o-dereplicated-table "$RES/derep-table.qza" \
        --o-dereplicated-sequences "$RES/derep-seqs.qza" \
        >> "$LOGS/otu.log" 2>&1
fi

# ---- 3. de novo chimera filtering (source study removes chimeras) ----
if [[ ! -s "$RES/uchime-stats.qza" ]]; then
    echo "$(NOW) [$COHORT] vsearch uchime-denovo" | tee -a "$LOGS/otu.log"
    "$QIIME_ENV/bin/qiime" vsearch uchime-denovo \
        --i-table "$RES/derep-table.qza" \
        --i-sequences "$RES/derep-seqs.qza" \
        --p-method uchime \
        --o-chimeras "$RES/chimeras.qza" \
        --o-nonchimeras "$RES/nonchimeras.qza" \
        --o-stats "$RES/uchime-stats.qza" \
        >> "$LOGS/otu.log" 2>&1
fi

if [[ ! -s "$RES/nonchimera-table.qza" ]]; then
    echo "$(NOW) [$COHORT] keeping non-chimeric features in dereplicated table" | tee -a "$LOGS/otu.log"
    # Keep only features present in the non-chimeric sequence set. uchime-denovo
    # skips some features (e.g. singletons) that are absent from BOTH its
    # chimera and non-chimera outputs, so excluding by chimera-id alone would
    # leave table features with no corresponding sequences.
    "$QIIME_ENV/bin/qiime" feature-table filter-features \
        --i-table "$RES/derep-table.qza" \
        --m-metadata-file "$RES/nonchimeras.qza" \
        --o-filtered-table "$RES/nonchimera-table.qza" \
        >> "$LOGS/otu.log" 2>&1
fi

# ---- 4. de novo 97% OTU clustering (source study: USEARCH at 97%) ----
if [[ ! -s "$RES/table.qza" ]]; then
    echo "$(NOW) [$COHORT] vsearch cluster-features-de-novo 97%" | tee -a "$LOGS/otu.log"
    "$QIIME_ENV/bin/qiime" vsearch cluster-features-de-novo \
        --i-table "$RES/nonchimera-table.qza" \
        --i-sequences "$RES/nonchimeras.qza" \
        --p-perc-identity 0.97 \
        --p-strand both \
        --p-threads 4 \
        --o-clustered-table "$RES/clustered-table.qza" \
        --o-clustered-sequences "$RES/rep-seqs.qza" \
        >> "$LOGS/otu.log" 2>&1
fi

# ---- 5. remove singleton OTUs (source study: OTUs with one read removed) ----
if [[ ! -s "$RES/table.qza" ]]; then
    echo "$(NOW) [$COHORT] removing singleton OTUs (min total frequency 2)" | tee -a "$LOGS/otu.log"
    "$QIIME_ENV/bin/qiime" feature-table filter-features \
        --i-table "$RES/clustered-table.qza" \
        --p-min-frequency 2 \
        --o-filtered-table "$RES/table.qza" \
        >> "$LOGS/otu.log" 2>&1
fi


# synchronize representative sequences with filtered OTU table
if [[ ! -s "$RES/rep-seqs-filtered.qza" ]]; then
    echo "$(NOW) [$COHORT] filtering rep-seqs to retained OTUs" | tee -a "$LOGS/otu.log"
    "$QIIME_ENV/bin/qiime" feature-table filter-seqs \
        --i-data "$RES/rep-seqs.qza" \
        --i-table "$RES/table.qza" \
        --o-filtered-data "$RES/rep-seqs-filtered.qza" \
        >> "$LOGS/otu.log" 2>&1
fi

# ---- 6. taxonomy classification (consensus-vsearch, SILVA 138.2) ----
# ---- 6. taxonomy classification (Greengenes 13_8; study-compatible OTU taxonomy) ----
if [[ ! -s "$RES/taxonomy.qza" ]]; then
    echo "$(NOW) [$COHORT] classifying retained OTUs with Greengenes 13_8" | tee -a "$LOGS/otu.log"
    "$QIIME_ENV/bin/qiime" feature-classifier classify-consensus-vsearch \
        --i-query "$RES/rep-seqs-filtered.qza" \
        --i-reference-reads "$HOME/OSCC_16S/data/greengenes13_8/gg13_8_97_seqs.qza" \
        --i-reference-taxonomy "$HOME/OSCC_16S/data/greengenes13_8/gg13_8_97_tax.qza" \
        --p-perc-identity 0.90 \
        --p-strand both \
        --p-maxaccepts 3 \
        --p-min-consensus 0.51 \
        --p-threads 8 \
        --o-classification "$RES/taxonomy.qza" \
        --o-search-results "$RES/taxonomy-search-results.qza" \
        >> "$LOGS/otu.log" 2>&1
fi

# ---- 7. exports ----
echo "$(NOW) [$COHORT] exporting results" | tee -a "$LOGS/otu.log"
"$QIIME_ENV/bin/qiime" tools export \
    --input-path "$RES/table.qza" \
    --output-path "$RES/export" \
    >> "$LOGS/otu.log" 2>&1

"$QIIME_ENV/bin/qiime" tools export \
    --input-path "$RES/rep-seqs-filtered.qza" \
    --output-path "$RES/export" \
    >> "$LOGS/otu.log" 2>&1

if [[ -f "$RES/export/feature-table.biom" && ! -f "$RES/export/table.biom" ]]; then
    mv "$RES/export/feature-table.biom" "$RES/export/table.biom"
fi
if [[ -f "$RES/export/dna-sequences.fasta" && ! -f "$RES/export/rep-seqs.fasta" ]]; then
    mv "$RES/export/dna-sequences.fasta" "$RES/export/rep-seqs.fasta"
fi

"$QIIME_ENV/bin/python" "$(dirname "$0")/qiime2_helpers.py" \
    table-to-tsv "$RES/table.qza" "$RES/export/table.tsv" \
    >> "$LOGS/otu.log" 2>&1

"$QIIME_ENV/bin/python" "$(dirname "$0")/qiime2_helpers.py" \
    taxa-to-tsv "$RES/taxonomy.qza" "$RES/export/taxonomy.tsv" \
    >> "$LOGS/otu.log" 2>&1

# ---- 8. filter + genus collapse (same rules as DADA2 cohorts) ----
echo "$(NOW) [$COHORT] filtering and collapsing to genus" | tee -a "$LOGS/otu.log"
"$QIIME_ENV/bin/python" "$(dirname "$0")/taxonomy_processing.py" \
    "$BASE" \
    "$RES/export/table.tsv" \
    "$RES/export/taxonomy.tsv" \
    >> "$LOGS/otu.log" 2>&1

# ---- 9. per-sample depth from OTU counts (LOW_DEPTH <2000, same rule) ----
DEPTH="$RES/depth_status.tsv"
"$QIIME_ENV/bin/python" - "$RES/export/table.tsv" "$DEPTH" <<'PYEOF' >> "$LOGS/otu.log" 2>&1
import sys
import pandas as pd
t = pd.read_csv(sys.argv[1], sep="\t", index_col=0)  # features x samples
samples = t.index.tolist()
counts = t.sum(axis=1)
depth = pd.DataFrame({
    "sample-id": samples,
    "input_reads": counts.values,
    "filtered_reads": counts.values,
    "denoised_reads": counts.values,
    "merged_reads": counts.values,
    "non_chimeric_reads": counts.values,
})
depth["depth_status"] = depth["non_chimeric_reads"].apply(
    lambda x: "OK" if x >= 2000 else "LOW_DEPTH")
depth.to_csv(sys.argv[2], sep="\t", index=False)
print("depth_status written", depth.shape)
print(depth[["sample-id", "non_chimeric_reads", "depth_status"]].to_string(index=False))
PYEOF

# ---- 10. paired genus statistics (same script as DADA2 cohorts) ----
echo "$(NOW) [$COHORT] paired genus statistics" | tee -a "$LOGS/otu.log"
"$QIIME_ENV/bin/python" "$(dirname "$0")/genus_stats.py" "$BASE" >> "$LOGS/otu.log" 2>&1

mark_done "$MARKER"
echo "$(NOW) [$COHORT] otu sensitivity COMPLETE" | tee -a "$LOGS/otu.log"
