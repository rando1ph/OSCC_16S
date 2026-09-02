#!/usr/bin/env bash
# QIIME2 DADA2 pipeline for one cohort.
# Resumable: if all four DADA2 artifacts already exist and validate, the
# denoising step is NEVER rerun; only the (idempotent) exports are produced.
# Exports always use `qiime tools export` -- never `qiime feature-table export`.
# usage: qiime2_pipeline.sh <COHORT>
set -euo pipefail
source "$(dirname "$0")/common.sh"

COHORT="$1"
BASE="$(cohort_dir "$COHORT")"
OUT="$BASE/processed"
RES="$BASE/results"
LOGS="$BASE/logs"
LOCK="$BASE/.qiime2.lock"
MARKER="$BASE/.qiime2_done"

export PATH="$QIIME_ENV/bin:$PATH"

mkdir -p "$RES" "$LOGS" "$RES/export"

acquire_lock "$LOCK" || { echo "cannot lock qiime2 for $COHORT"; exit 1; }

if stage_done "$MARKER"; then
    echo "qiime2 already complete for $COHORT"
    exit 0
fi

# ---- build processed manifest from analysis metadata ----
PMAN="$OUT/manifest.tsv"
{
    printf 'sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n'
    sed 's/\r$//' "$BASE/metadata/analysis_metadata.tsv" \
        | tail -n +2 \
        | cut -f1 \
        | while read -r sample; do
            printf '%s\t%s\t%s\n' \
                "$sample" \
                "$OUT/${sample}_1.fastq.gz" \
                "$OUT/${sample}_2.fastq.gz"
        done
} > "$PMAN"

nfiles=$(tail -n +2 "$PMAN" | wc -l)
missing=0
while IFS=$'\t' read -r sample r1 r2; do
    [[ "$sample" == "sample-id" ]] && continue
    if [[ ! -s "$r1" || ! -s "$r2" ]]; then
        echo "MISSING $sample" | tee -a "$LOGS/qiime2.log"
        missing=$((missing + 1))
    fi
done < "$PMAN"

echo "$(NOW) [$COHORT] qiime2: manifest samples=$nfiles missing=$missing" \
    | tee -a "$LOGS/qiime2.log"
if [[ "$missing" -gt 0 ]]; then
    echo "$(NOW) [$COHORT] FAIL: $missing samples have no processed FASTQs" \
        | tee -a "$LOGS/qiime2.log"
    mark_fail "$BASE/.qiime2_failed"
    exit 1
fi

# ---- DADA2: skip entirely if all four artifacts already exist ----
NEED_DADA2=0
for a in table rep-seqs stats base-transition-stats; do
    if [[ ! -s "$RES/$a.qza" ]]; then
        NEED_DADA2=1
    fi
done

if [[ "$NEED_DADA2" -eq 1 ]]; then
    if [[ ! -s "$RES/demux.qza" ]]; then
        echo "$(NOW) [$COHORT] importing demux" | tee -a "$LOGS/qiime2.log"
        "$QIIME_ENV/bin/qiime" tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path "$PMAN" \
            --output-path "$RES/demux.qza" \
            --input-format PairedEndFastqManifestPhred33V2 \
            >> "$LOGS/qiime2.log" 2>&1
    fi

    echo "$(NOW) [$COHORT] running dada2 denoise-paired (trunc-len 0, n-threads 4)" \
        | tee -a "$LOGS/qiime2.log"
    "$QIIME_ENV/bin/qiime" dada2 denoise-paired \
        --i-demultiplexed-seqs "$RES/demux.qza" \
        --p-trunc-len-f 0 \
        --p-trunc-len-r 0 \
        --p-trim-left-f 0 \
        --p-trim-left-r 0 \
        --p-n-threads 4 \
        --p-chimera-method consensus \
        --o-table "$RES/table.qza" \
        --o-representative-sequences "$RES/rep-seqs.qza" \
        --o-denoising-stats "$RES/stats.qza" \
        --o-base-transition-stats "$RES/base-transition-stats.qza" \
        >> "$LOGS/qiime2.log" 2>&1
else
    echo "$(NOW) [$COHORT] all four DADA2 artifacts present; preserving and NOT rerunning DADA2" \
        | tee -a "$LOGS/qiime2.log"
fi

# ---- exports (always `qiime tools export`, idempotent) ----
echo "$(NOW) [$COHORT] exporting results" | tee -a "$LOGS/qiime2.log"
"$QIIME_ENV/bin/qiime" tools export \
    --input-path "$RES/table.qza" \
    --output-path "$RES/export" \
    >> "$LOGS/qiime2.log" 2>&1

"$QIIME_ENV/bin/qiime" tools export \
    --input-path "$RES/rep-seqs.qza" \
    --output-path "$RES/export" \
    >> "$LOGS/qiime2.log" 2>&1

"$QIIME_ENV/bin/qiime" tools export \
    --input-path "$RES/stats.qza" \
    --output-path "$RES/export" \
    >> "$LOGS/qiime2.log" 2>&1

if [[ -f "$RES/export/feature-table.biom" && ! -f "$RES/export/table.biom" ]]; then
    mv "$RES/export/feature-table.biom" "$RES/export/table.biom"
fi
if [[ -f "$RES/export/dna-sequences.fasta" && ! -f "$RES/export/rep-seqs.fasta" ]]; then
    mv "$RES/export/dna-sequences.fasta" "$RES/export/rep-seqs.fasta"
fi

"$QIIME_ENV/bin/python" "$(dirname "$0")/qiime2_helpers.py" \
    table-to-tsv "$RES/table.qza" "$RES/export/table.tsv" \
    >> "$LOGS/qiime2.log" 2>&1

# ---- parse denoising stats into retention file ----
DADA2_RET="$RES/dada2_retention.tsv"
"$QIIME_ENV/bin/python" - "$RES/export/stats.tsv" "$DADA2_RET" <<'PYEOF' >> "$LOGS/qiime2.log" 2>&1
import sys
import pandas as pd
stats = pd.read_csv(sys.argv[1], sep="\t", index_col=0, comment="#")
stats = stats.rename(columns={
    "input": "input_reads",
    "filtered": "filtered_reads",
    "percentage of input passed filter": "pct_passed_filter",
    "denoised": "denoised_reads",
    "merged": "merged_reads",
    "non-chimeric": "non_chimeric_reads",
    "percentage of input non-chimeric": "pct_non_chimeric",
})
keep = [c for c in ["input_reads","filtered_reads","denoised_reads",
                    "merged_reads","non_chimeric_reads"] if c in stats.columns]
stats[keep].to_csv(sys.argv[2], sep="\t")
print("wrote", sys.argv[2], stats.shape)
PYEOF

# ---- depth status ----
DEPTH="$RES/depth_status.tsv"
"$QIIME_ENV/bin/python" - "$DADA2_RET" "$DEPTH" <<'PYEOF' >> "$LOGS/qiime2.log" 2>&1
import sys
import pandas as pd
d = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
d["depth_status"] = d["non_chimeric_reads"].apply(
    lambda x: "OK" if x >= 2000 else "LOW_DEPTH")
d.to_csv(sys.argv[2], sep="\t")
print("depth_status written", d.shape)
PYEOF

mark_done "$MARKER"
echo "$(NOW) [$COHORT] qiime2 pipeline COMPLETE" | tee -a "$LOGS/qiime2.log"
