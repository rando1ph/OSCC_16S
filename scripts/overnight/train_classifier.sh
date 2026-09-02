#!/usr/bin/env bash
# Train a shared SILVA 138.2 V3-V4 Naive Bayes classifier with RESCRIPt.
# Runs once; result cached in data/classifiers/.
set -euo pipefail
source "$(dirname "$0")/common.sh"

DATA="$ROOT/data"
CLS_DIR="$DATA/classifiers"
SILVA="$DATA/silva"
CLASSIFIER="$CLS_DIR/silva-138.2-v3v4-nb-classifier.qza"
LOCK="$DATA/.train_classifier.lock"
MARKER="$CLS_DIR/.classifier_done"
LOG="$DATA/train_classifier.log"

export PATH="$QIIME_ENV/bin:$PATH"

mkdir -p "$CLS_DIR" "$SILVA"

acquire_lock "$LOCK" || { echo "cannot lock classifier training"; exit 1; }

if stage_done "$MARKER"; then
    echo "classifier already trained: $CLASSIFIER"
    exit 0
fi

trap 'mark_fail "$CLS_DIR/.classifier_failed"; echo "$(NOW) classifier training FAILED" >> "$LOG"' ERR

echo "$(NOW) classifier training: fetching SILVA 138.2 SSURef_NR99" | tee -a "$LOG"
"$QIIME_ENV/bin/qiime" rescript get-silva-data \
    --p-version '138.2' \
    --p-target 'SSURef_NR99' \
    --p-no-include-species-labels \
    --o-silva-sequences "$SILVA/silva-138.2-seqs.qza" \
    --o-silva-taxonomy "$SILVA/silva-138.2-tax.qza" \
    >> "$LOG" 2>&1

echo "$(NOW) culling sequences" | tee -a "$LOG"
"$QIIME_ENV/bin/qiime" rescript cull-seqs \
    --i-sequences "$SILVA/silva-138.2-seqs.qza" \
    --p-num-degenerates 5 \
    --p-homopolymer-length 8 \
    --o-clean-sequences "$SILVA/silva-138.2-culled.qza" \
    >> "$LOG" 2>&1

echo "$(NOW) filtering by length/taxon" | tee -a "$LOG"
"$QIIME_ENV/bin/qiime" rescript filter-seqs-length-by-taxon \
    --i-sequences "$SILVA/silva-138.2-culled.qza" \
    --i-taxonomy "$SILVA/silva-138.2-tax.qza" \
    --p-min-len 900 \
    --p-max-len 1900 \
    --o-filtered-sequences "$SILVA/silva-138.2-filtered.qza" \
    --o-filtered-taxonomy "$SILVA/silva-138.2-filtered-tax.qza" \
    >> "$LOG" 2>&1

echo "$(NOW) dereplicating" | tee -a "$LOG"
"$QIIME_ENV/bin/qiime" rescript dereplicate \
    --i-sequences "$SILVA/silva-138.2-filtered.qza" \
    --i-taxonomy "$SILVA/silva-138.2-filtered-tax.qza" \
    --p-mode uniq \
    --o-dereplicated-sequences "$SILVA/silva-138.2-derep.qza" \
    --o-dereplicated-taxonomy "$SILVA/silva-138.2-derep-tax.qza" \
    >> "$LOG" 2>&1

echo "$(NOW) extracting V3-V4 reads" | tee -a "$LOG"
"$QIIME_ENV/bin/qiime" feature-classifier extract-reads \
    --i-sequences "$SILVA/silva-138.2-derep.qza" \
    --p-f-primer 'CCTACGGGNGGCWGCAG' \
    --p-r-primer 'GACTACHVGGGTATCTAATCC' \
    --p-n-jobs 4 \
    --p-read-orientation both \
    --o-reads "$SILVA/silva-138.2-v3v4-reads.qza" \
    >> "$LOG" 2>&1

echo "$(NOW) fitting Naive Bayes classifier" | tee -a "$LOG"
"$QIIME_ENV/bin/qiime" feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$SILVA/silva-138.2-v3v4-reads.qza" \
    --i-reference-taxonomy "$SILVA/silva-138.2-derep-tax.qza" \
    --p-n-jobs 4 \
    --o-classifier "$CLASSIFIER" \
    >> "$LOG" 2>&1

mark_done "$MARKER"
echo "$(NOW) classifier training COMPLETE: $CLASSIFIER" | tee -a "$LOG"
