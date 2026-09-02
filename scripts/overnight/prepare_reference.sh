#!/usr/bin/env bash
# Prepare SILVA 138.2 reference reads + taxonomy for consensus-vsearch
# classification, using RESCRIPt on ALREADY-DOWNLOADED artifacts.
# Never downloads anything new; idempotent and resumable.
# usage: prepare_reference.sh
set -euo pipefail
source "$(dirname "$0")/common.sh"

SILVA="$ROOT/data/silva"
LOCK="$SILVA/.prepare_reference.lock"
MARKER="$SILVA/.reference_done"
LOG="$SILVA/prepare_reference.log"

export PATH="$QIIME_ENV/bin:$PATH"

mkdir -p "$SILVA"

acquire_lock "$LOCK" || { echo "cannot lock reference preparation"; exit 1; }

if stage_done "$MARKER"; then
    echo "reference already prepared: $SILVA/silva-138.2-derep.qza"
    exit 0
fi

trap 'mark_fail "$SILVA/.reference_failed"; echo "$(NOW) reference preparation FAILED" >> "$LOG"' ERR

SEQS="$SILVA/silva-138.2-seqs.qza"
TAX="$SILVA/silva-138.2-tax.qza"
CULLED="$SILVA/silva-138.2-culled.qza"

if [[ ! -f "$SEQS" || ! -f "$TAX" ]]; then
    echo "$(NOW) FATAL: downloaded SILVA artifacts missing ($SEQS, $TAX)" | tee -a "$LOG"
    mark_fail "$SILVA/.reference_failed"
    exit 1
fi

if [[ ! -f "$CULLED" ]]; then
    echo "$(NOW) culling sequences" | tee -a "$LOG"
    "$QIIME_ENV/bin/qiime" rescript cull-seqs \
        --i-sequences "$SEQS" \
        --p-num-degenerates 5 \
        --p-homopolymer-length 8 \
        --o-clean-sequences "$CULLED" \
        >> "$LOG" 2>&1
fi

# Dereplicate the culled reference against the downloaded taxonomy.
# The taxonomy artifact is a superset of the culled sequence IDs, which
# rescript dereplicate accepts (extra taxonomy entries are ignored).
echo "$(NOW) dereplicating reference" | tee -a "$LOG"
"$QIIME_ENV/bin/qiime" rescript dereplicate \
    --i-sequences "$CULLED" \
    --i-taxa "$TAX" \
    --p-mode uniq \
    --p-threads 4 \
    --o-dereplicated-sequences "$SILVA/silva-138.2-derep.qza" \
    --o-dereplicated-taxa "$SILVA/silva-138.2-derep-tax.qza" \
    >> "$LOG" 2>&1

mark_done "$MARKER"
echo "$(NOW) reference preparation COMPLETE: $SILVA/silva-138.2-derep.qza" | tee -a "$LOG"
