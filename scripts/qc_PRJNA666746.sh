#!/usr/bin/env bash
set -u

BASE="$HOME/OSCC_16S/PRJNA666746"
RAW="$BASE/raw"
META="$BASE/metadata"
RESULTS="$BASE/results"
LOGS="$BASE/logs"

mkdir -p "$RESULTS/fastqc" "$RESULTS/multiqc" "$LOGS"

cd "$RAW"

fastq_count=$(find . -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)
partial_count=$(find . -maxdepth 1 -type f -name '*.aria2' | wc -l)

echo "FASTQ files: $fastq_count"
echo "Partial files: $partial_count"

if [ "$fastq_count" -ne 200 ]; then
    echo "ERROR: expected 200 FASTQ files."
    exit 1
fi

if [ "$partial_count" -ne 0 ]; then
    echo "ERROR: download is not finished."
    exit 1
fi

echo "=== MD5 CHECK ==="

if md5sum -c "$META/fastq_md5.txt" > "$LOGS/md5_check.log" 2>&1; then
    echo "MD5: PASS"
else
    echo "MD5: FAIL"
    tail -20 "$LOGS/md5_check.log"
    exit 1
fi

echo "=== FASTQC ==="

fastqc \
    -t 4 \
    -o "$RESULTS/fastqc" \
    ./*.fastq.gz \
    > "$LOGS/fastqc.log" 2>&1

echo "FastQC: DONE"

echo "=== MULTIQC ==="

"$HOME/.local/bin/multiqc" \
    "$RESULTS/fastqc" \
    -o "$RESULTS/multiqc" \
    -f \
    > "$LOGS/multiqc.log" 2>&1

echo "MultiQC: DONE"

echo "=== QC COMPLETE ==="
echo "FastQC results: $RESULTS/fastqc"
echo "MultiQC report: $RESULTS/multiqc/multiqc_report.html"
