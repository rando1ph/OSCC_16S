#!/usr/bin/env bash
# Download + MD5 verify one cohort's raw FASTQs.
# usage: download.sh <COHORT>
set -euo pipefail
source "$(dirname "$0")/common.sh"

COHORT="$1"
BASE="$(cohort_dir "$COHORT")"
RAW="$BASE/raw"
META="$BASE/metadata"
LOGS="$BASE/logs"
URLS="$META/fastq_urls.txt"
MD5F="$META/fastq_md5.txt"
LOCK="$BASE/.download.lock"
STATUS="$LOGS/download_status.txt"
DONE_MARKER="$BASE/.download_done"

mkdir -p "$RAW" "$LOGS"

{
    acquire_lock "$LOCK" || { echo "cannot lock $LOCK"; exit 1; }

    if stage_done "$DONE_MARKER"; then
        echo "download already complete for $COHORT"
        exit 0
    fi

    expected=$(wc -l < "$URLS")
    echo "DOWNLOADING" > "$STATUS"

    attempt=1
    ok=0
    while [[ $attempt -le 3 && $ok -eq 0 ]]; do
        echo "$(NOW) [$COHORT] download attempt $attempt ($expected urls, 2 concurrent)" \
            | tee -a "$LOGS/download_run.log"
        input="$URLS"
        if [[ $attempt -gt 1 && -s "$LOGS/failed_urls.txt" ]]; then
            input="$LOGS/failed_urls.txt"
        fi
        aria2c \
            -i "$input" \
            -d "$RAW" \
            -j 2 \
            -x 4 -s 4 -k 1M \
            -c \
            --max-tries=5 \
            --retry-wait=5 \
            --timeout=60 \
            --console-log-level=warn \
            --summary-interval=0 \
            >> "$LOGS/download_run.log" 2>&1 || true

        present=$(find "$RAW" -maxdepth 1 -name '*.fastq.gz' | wc -l)
        partial=$(find "$RAW" -maxdepth 1 -name '*.aria2' | wc -l)
        echo "$(NOW) [$COHORT] present=$present partial=$partial expected=$expected" \
            | tee -a "$LOGS/download_run.log"

        if [[ "$present" -ne "$expected" || "$partial" -ne 0 ]]; then
            echo "$(NOW) [$COHORT] incomplete download (attempt $attempt)" \
                | tee -a "$LOGS/download_run.log"
            attempt=$((attempt + 1))
            continue
        fi

        # integrity check
        echo "VERIFYING" > "$STATUS"
        if ( cd "$RAW" && md5sum -c "$MD5F" ) > "$LOGS/md5_check.log" 2>&1; then
            echo "DONE" > "$STATUS"
            mark_done "$DONE_MARKER"
            rm -f "$LOGS/failed_urls.txt"
            echo "$(NOW) [$COHORT] MD5 PASS, download complete" \
                | tee -a "$LOGS/download_run.log"
            ok=1
        else
            grep -E 'FAILED|No such file' "$LOGS/md5_check.log" \
                | awk '{print $1}' | sed 's/:$//' | sort -u > "$LOGS/failed_files.txt"
            # build failed urls
            : > "$LOGS/failed_urls.txt"
            while IFS= read -r fname; do
                [ -z "$fname" ] && continue
                grep "/$fname$" "$URLS" >> "$LOGS/failed_urls.txt" || true
            done < "$LOGS/failed_files.txt"
            echo "$(NOW) [$COHORT] MD5 FAIL attempt $attempt (failed: $(wc -l < "$LOGS/failed_files.txt"))" \
                | tee -a "$LOGS/download_run.log"
            attempt=$((attempt + 1))
        fi
    done

    if [[ $ok -eq 0 ]]; then
        echo "FAILED" > "$STATUS"
        mark_fail "$BASE/.download_failed"
        echo "$(NOW) [$COHORT] DOWNLOAD FAILED after 3 attempts" \
            | tee -a "$LOGS/download_run.log"
        exit 1
    fi
}
