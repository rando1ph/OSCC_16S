#!/usr/bin/env bash
# Taxonomy classification + filtering + genus collapse for one cohort.
# Uses qiime feature-classifier classify-consensus-vsearch against a
# SILVA 138.2 reference prepared by prepare_reference.sh (no sklearn training).
# usage: taxonomy.sh <COHORT>
set -euo pipefail
source "$(dirname "$0")/common.sh"

COHORT="$1"
BASE="$(cohort_dir "$COHORT")"
RES="$BASE/results"
LOGS="$BASE/logs"
LOCK="$BASE/.taxonomy.lock"
MARKER="$BASE/.taxonomy_done"
REF_READS="$ROOT/data/silva/silva-138.2-derep.qza"
REF_TAX="$ROOT/data/silva/silva-138.2-derep-tax.qza"
REF_DONE="$ROOT/data/silva/.reference_done"

export PATH="$QIIME_ENV/bin:$PATH"

mkdir -p "$RES" "$LOGS"

acquire_lock "$LOCK" || { echo "cannot lock taxonomy for $COHORT"; exit 1; }

if stage_done "$MARKER"; then
    echo "taxonomy already complete for $COHORT"
    exit 0
fi

if ! stage_done "$REF_DONE"; then
    echo "$(NOW) [$COHORT] SILVA reference not prepared yet: $REF_DONE" | tee -a "$LOGS/taxonomy.log"
    mark_fail "$BASE/.taxonomy_failed"
    exit 1
fi

if [[ ! -f "$REF_READS" || ! -f "$REF_TAX" ]]; then
    echo "$(NOW) [$COHORT] reference artifacts missing: $REF_READS / $REF_TAX" | tee -a "$LOGS/taxonomy.log"
    mark_fail "$BASE/.taxonomy_failed"
    exit 1
fi

if [[ ! -f "$RES/rep-seqs.qza" ]]; then
    echo "$(NOW) [$COHORT] rep-seqs.qza missing (DADA2 incomplete)" | tee -a "$LOGS/taxonomy.log"
    mark_fail "$BASE/.taxonomy_failed"
    exit 1
fi

echo "$(NOW) [$COHORT] classifying rep-seqs with consensus-vsearch" | tee -a "$LOGS/taxonomy.log"
"$QIIME_ENV/bin/qiime" feature-classifier classify-consensus-vsearch \
    --i-query "$RES/rep-seqs.qza" \
    --i-reference-reads "$REF_READS" \
    --i-reference-taxonomy "$REF_TAX" \
    --p-perc-identity 0.8 \
    --p-strand both \
    --p-maxaccepts 3 \
    --p-min-consensus 0.51 \
    --p-threads 4 \
    --o-classification "$RES/taxonomy.qza" \
    --o-search-results "$RES/taxonomy-search-results.qza" \
    >> "$LOGS/taxonomy.log" 2>&1

echo "$(NOW) [$COHORT] exporting taxonomy" | tee -a "$LOGS/taxonomy.log"
"$QIIME_ENV/bin/python" "$(dirname "$0")/qiime2_helpers.py" \
    taxa-to-tsv "$RES/taxonomy.qza" "$RES/export/taxonomy.tsv" \
    >> "$LOGS/taxonomy.log" 2>&1

echo "$(NOW) [$COHORT] filtering and collapsing to genus" | tee -a "$LOGS/taxonomy.log"
"$QIIME_ENV/bin/python" "$(dirname "$0")/taxonomy_processing.py" \
    "$BASE" \
    "$RES/export/table.tsv" \
    "$RES/export/taxonomy.tsv" \
    >> "$LOGS/taxonomy.log" 2>&1

mark_done "$MARKER"
echo "$(NOW) [$COHORT] taxonomy COMPLETE" | tee -a "$LOGS/taxonomy.log"
