#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/OSCC_16S/PRJNA666746"
RAW="$BASE/raw"
META="$BASE/metadata"
OUT="$BASE/processed"
LOG="$BASE/logs/preprocess_PRJNA666746.log"

CUTADAPT="$HOME/.local/bin/cutadapt"

FWD='CCTACGGGNGGCWGCAG'
REV='GACTACHVGGGTATCTAATCC'

mkdir -p "$OUT"
: > "$LOG"

echo "=== PRJNA666746 preprocessing ===" | tee -a "$LOG"
echo "Started: $(date)" | tee -a "$LOG"

count=0

tail -n +2 "$META/manifest.tsv" | while IFS=$'\t' read -r sample r1 r2; do

    count=$((count+1))

    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "ERROR: missing FASTQ for $sample" | tee -a "$LOG"
        exit 1
    fi

    if [[ -f "${r1}.aria2" || -f "${r2}.aria2" ]]; then
        echo "ERROR: unfinished download for $sample" | tee -a "$LOG"
        exit 1
    fi

    echo "[$count/100] $sample" | tee -a "$LOG"

    tmp="$OUT/.tmp_${sample}"
    mkdir -p "$tmp"

    # Orientation 1:
    # raw R1 starts with forward primer
    # raw R2 starts with reverse primer
    "$CUTADAPT" \
        -j 4 \
        -e 0.15 \
        --no-indels \
        -g "^$FWD" \
        -G "^$REV" \
        --discard-untrimmed \
        -o "$tmp/FR_R1.fastq.gz" \
        -p "$tmp/FR_R2.fastq.gz" \
        "$r1" "$r2" \
        > "$tmp/FR.log"

    # Orientation 2:
    # raw R1 starts with reverse primer
    # raw R2 starts with forward primer
    "$CUTADAPT" \
        -j 4 \
        -e 0.15 \
        --no-indels \
        -g "^$REV" \
        -G "^$FWD" \
        --discard-untrimmed \
        -o "$tmp/RF_R1.fastq.gz" \
        -p "$tmp/RF_R2.fastq.gz" \
        "$r1" "$r2" \
        > "$tmp/RF.log"

    # Normalize orientation:
    # FR: R1 -> forward, R2 -> reverse
    # RF: R2 -> forward, R1 -> reverse

    {
        pigz -dc "$tmp/FR_R1.fastq.gz"
        pigz -dc "$tmp/RF_R2.fastq.gz"
    } | pigz -p 2 > "$OUT/${sample}_1.fastq.gz"

    {
        pigz -dc "$tmp/FR_R2.fastq.gz"
        pigz -dc "$tmp/RF_R1.fastq.gz"
    } | pigz -p 2 > "$OUT/${sample}_2.fastq.gz"

    fr=$(grep 'Pairs written (passing filters)' "$tmp/FR.log" | awk '{gsub(",","",$5); print $5}')
    rf=$(grep 'Pairs written (passing filters)' "$tmp/RF.log" | awk '{gsub(",","",$5); print $5}')

    echo -e "$sample\tFR=$fr\tRF=$rf" >> "$LOG"

    rm -rf "$tmp"

done

echo "Finished: $(date)" | tee -a "$LOG"
echo "=== COMPLETE ===" | tee -a "$LOG"
