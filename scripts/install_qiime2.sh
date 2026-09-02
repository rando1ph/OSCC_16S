#!/usr/bin/env bash
set -euo pipefail

source "$HOME/miniconda3/etc/profile.d/conda.sh"

LOG="$HOME/OSCC_16S/qiime2_install.log"

echo "===== QIIME2 INSTALL START $(date) =====" >> "$LOG"

conda env create \
  --name rachis-qiime2-2026.7 \
  --file https://raw.githubusercontent.com/qiime2/distributions/refs/heads/dev/2026.7/qiime2/released/rachis-qiime2-linux-64-conda.yml \
  >> "$LOG" 2>&1

conda activate rachis-qiime2-2026.7

echo "===== QIIME INFO =====" >> "$LOG"
qiime info >> "$LOG" 2>&1

echo "===== QIIME2 INSTALL COMPLETE $(date) =====" >> "$LOG"
