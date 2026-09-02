#!/usr/bin/env bash
# Preprocess ONE sample: deterministic cap, primer trimming, orientation
# normalization, optional 3' trim, final cap. Writes retention row.
# usage: preprocess_sample.sh <COHORT> <sample> <r1> <r2>
set -euo pipefail
source "$(dirname "$0")/common.sh"

COHORT="$1"
LINE="$2"
IFS=$'\t' read -r SAMPLE R1 R2 <<< "$LINE"

BASE="$(cohort_dir "$COHORT")"
RAW="$BASE/raw"
META="$BASE/metadata"
OUT="$BASE/processed"
LOGS="$BASE/logs"
WORK="$OUT/.work"
RETENTION="$OUT/retention.tsv"
PRETRIM_CAP=150000
FINAL_CAP=100000
TAIL_TRIM=0

case "$COHORT" in
    PRJNA666746)
        FWD='CCTACGGGNGGCWGCAG'
        REV='GACTACHVGGGTATCTAATCC'
        TAIL_TRIM=0
        ;;
    PRJNA822685)
        FWD='CCTACGGGNGGCWGCAG'
        REV='GGACTACNVGGGTWTCTAAT'
        TAIL_TRIM=0
        ;;
    PRJNA813034)
        FWD='CCTACGGGNBGCASCAG'
        REV='GACTACNVGGGTATCTAATCC'
        TAIL_TRIM=10
        ;;
    *)
        echo "unknown cohort $COHORT" >&2
        exit 1
        ;;
esac

mkdir -p "$OUT" "$WORK" "$LOGS"

if awk -F'\t' -v s="$SAMPLE" '$1==s{f=1} END{exit !f}' "$RETENTION" 2>/dev/null; then
    echo "[$COHORT] $SAMPLE already processed"
    exit 0
fi

if [[ ! -f "$R1" || ! -f "$R2" ]]; then
    echo "$(NOW) [$COHORT] ERROR missing raw files for $SAMPLE ($R1, $R2)" >> "$LOGS/preprocess.log"
    exit 1
fi

echo "$(NOW) [$COHORT] processing $SAMPLE" >> "$LOGS/preprocess.log"

# stage 1: cap raw pairs
"$QIIME_ENV/bin/python" "$(dirname "$0")/subsample_pairs.py" \
    --r1 "$R1" --r2 "$R2" \
    --out1 "$WORK/${SAMPLE}_pre1.fastq.gz" \
    --out2 "$WORK/${SAMPLE}_pre2.fastq.gz" \
    --n "$PRETRIM_CAP" --seed "$SEED" \
    > "$WORK/${SAMPLE}.subsample.log" 2>&1

raw_pairs=$(awk -F'\t' 'NR==2{print $1}' "$WORK/${SAMPLE}_pre1.fastq.gz.meta.tsv")
pretrim_pairs=$(awk -F'\t' 'NR==2{print $2}' "$WORK/${SAMPLE}_pre1.fastq.gz.meta.tsv")

# stage 2: FR and RF primer trimming
"$CUTADAPT" -j 3 -e 0.15 --no-indels \
    -g "^$FWD" -G "^$REV" --discard-untrimmed \
    -o "$WORK/${SAMPLE}_FR_R1.fastq.gz" -p "$WORK/${SAMPLE}_FR_R2.fastq.gz" \
    "$WORK/${SAMPLE}_pre1.fastq.gz" "$WORK/${SAMPLE}_pre2.fastq.gz" \
    > "$WORK/${SAMPLE}_FR.log" 2>&1

"$CUTADAPT" -j 3 -e 0.15 --no-indels \
    -g "^$REV" -G "^$FWD" --discard-untrimmed \
    -o "$WORK/${SAMPLE}_RF_R1.fastq.gz" -p "$WORK/${SAMPLE}_RF_R2.fastq.gz" \
    "$WORK/${SAMPLE}_pre1.fastq.gz" "$WORK/${SAMPLE}_pre2.fastq.gz" \
    > "$WORK/${SAMPLE}_RF.log" 2>&1

fr=$(grep 'Pairs written (passing filters)' "$WORK/${SAMPLE}_FR.log" | awk '{gsub(",","",$5); print $5}')
rf=$(grep 'Pairs written (passing filters)' "$WORK/${SAMPLE}_RF.log" | awk '{gsub(",","",$5); print $5}')
fr=${fr:-0}; rf=${rf:-0}
valid=$((fr + rf))

# stage 3: normalize orientation
{
    pigz -dc "$WORK/${SAMPLE}_FR_R1.fastq.gz"
    pigz -dc "$WORK/${SAMPLE}_RF_R2.fastq.gz"
} | pigz -p 2 > "$WORK/${SAMPLE}_n1.fastq.gz"

{
    pigz -dc "$WORK/${SAMPLE}_FR_R2.fastq.gz"
    pigz -dc "$WORK/${SAMPLE}_RF_R1.fastq.gz"
} | pigz -p 2 > "$WORK/${SAMPLE}_n2.fastq.gz"

# stage 4: optional 3' trim
if [[ "$TAIL_TRIM" -gt 0 ]]; then
    "$CUTADAPT" -j 3 -u -"$TAIL_TRIM" -U -"$TAIL_TRIM" \
        -o "$WORK/${SAMPLE}_t1.fastq.gz" -p "$WORK/${SAMPLE}_t2.fastq.gz" \
        "$WORK/${SAMPLE}_n1.fastq.gz" "$WORK/${SAMPLE}_n2.fastq.gz" \
        > "$WORK/${SAMPLE}_trim.log" 2>&1
    mv "$WORK/${SAMPLE}_t1.fastq.gz" "$WORK/${SAMPLE}_n1.fastq.gz"
    mv "$WORK/${SAMPLE}_t2.fastq.gz" "$WORK/${SAMPLE}_n2.fastq.gz"
fi

# stage 5: cap valid pairs
"$QIIME_ENV/bin/python" "$(dirname "$0")/subsample_pairs.py" \
    --r1 "$WORK/${SAMPLE}_n1.fastq.gz" --r2 "$WORK/${SAMPLE}_n2.fastq.gz" \
    --out1 "$OUT/${SAMPLE}_1.fastq.gz" \
    --out2 "$OUT/${SAMPLE}_2.fastq.gz" \
    --n "$FINAL_CAP" --seed "$SEED" \
    > "$WORK/${SAMPLE}.final.log" 2>&1

final=$(awk -F'\t' 'NR==2{print $2}' "$OUT/${SAMPLE}_1.fastq.gz.meta.tsv")

frac=$(awk -v v="$valid" -v p="$pretrim_pairs" 'BEGIN{ if(p>0) printf "%.4f", v/p; else print "0" }')

# single-line append is atomic enough across workers (O_APPEND, <PIPE_BUF)
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE" "$raw_pairs" "$pretrim_pairs" "$fr" "$rf" "$valid" "$frac" "$final" \
    >> "$RETENTION"

rm -f "$WORK/${SAMPLE}_pre1.fastq.gz" "$WORK/${SAMPLE}_pre2.fastq.gz" \
    "$WORK/${SAMPLE}_pre1.fastq.gz.meta.tsv" \
    "$WORK/${SAMPLE}_FR_R1.fastq.gz" "$WORK/${SAMPLE}_FR_R2.fastq.gz" \
    "$WORK/${SAMPLE}_RF_R1.fastq.gz" "$WORK/${SAMPLE}_RF_R2.fastq.gz" \
    "$WORK/${SAMPLE}_n1.fastq.gz" "$WORK/${SAMPLE}_n2.fastq.gz"

echo "$(NOW) [$COHORT] $SAMPLE done (raw=$raw_pairs valid=$valid final=$final)" \
    >> "$LOGS/preprocess.log"
