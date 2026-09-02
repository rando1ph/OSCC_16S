#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/OSCC_16S/PRJNA666746"
RAW="$BASE/raw"
META="$BASE/metadata"
LOG="$BASE/logs/md5_check.log"
STATUS="$BASE/logs/md5_status.txt"

cd "$RAW"

echo "RUNNING" > "$STATUS"
echo "MD5 check started: $(date)" > "$LOG"

if md5sum -c "$META/fastq_md5.txt" >> "$LOG" 2>&1; then
    echo "PASS" > "$STATUS"
    echo "MD5 check PASSED: $(date)" >> "$LOG"
else
    echo "FAIL" > "$STATUS"
    echo "MD5 check FAILED: $(date)" >> "$LOG"
    exit 1
fi
