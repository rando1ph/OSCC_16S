#!/usr/bin/env bash

RAW="$HOME/OSCC_16S/PRJNA666746/raw"
LOG="$HOME/OSCC_16S/PRJNA666746/logs/download_watch.log"

while true; do
    total=$(find "$RAW" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)
    partial=$(find "$RAW" -maxdepth 1 -type f -name '*.aria2' | wc -l)
    size=$(du -sh "$RAW" | cut -f1)

    echo "$(date '+%Y-%m-%d %H:%M:%S') FASTQ=$total Partial=$partial Size=$size" >> "$LOG"

    if [ "$total" -eq 200 ] && [ "$partial" -eq 0 ]; then
        echo "$(date '+%Y-%m-%d %H:%M:%S') DOWNLOAD COMPLETE" >> "$LOG"
        break
    fi

    sleep 120
done
