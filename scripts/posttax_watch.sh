#!/usr/bin/env bash
# Post-taxonomy watcher + downstream analysis for PRJNA666746 and PRJNA822685.
# Waits for combined/tax_666746.exit and combined/tax_822685.exit (written by
# the independent twotax tmux session), then runs the existing downstream
# entrypoints exactly as encoded in AUTOMATION_TASK.md and scripts/overnight/.
# Resumable: flock-protected, stage markers, idempotent downstream scripts.
# usage: bash scripts/posttax_watch.sh
set -u

ROOT="${POSTTAX_ROOT:-$HOME/OSCC_16S}"
QIIME_ENV="$HOME/miniconda3/envs/rachis-qiime2-2026.7"
PY="$QIIME_ENV/bin/python"
COHORTS="PRJNA666746 PRJNA822685"

COMBINED="$ROOT/combined"
POSTTAX_LOG="${POSTTAX_LOG:-$COMBINED/posttax.log}"
POSTTAX_STATUS="${POSTTAX_STATUS:-$COMBINED/posttax_status.tsv}"
WATCH_INTERVAL="${WATCH_INTERVAL:-60}"
LOCK="$COMBINED/.posttax.lock"

EXIT_666746="$COMBINED/tax_666746.exit"
EXIT_822685="$COMBINED/tax_822685.exit"

MARKER_STATS="$COMBINED/.posttax_stats_done"
MARKER_CROSS="$COMBINED/.posttax_cross_done"
MARKER_DONE="$COMBINED/.posttax_done"

mkdir -p "$COMBINED"

exec 9>"$LOCK"
if ! flock -n 9; then
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] [FATAL] another posttax watcher is running ($LOCK)" >> "$POSTTAX_LOG"
    exit 1
fi

NOW() { date '+%Y-%m-%d %H:%M:%S'; }

plog() {
    local level="$1"; shift
    echo "[$(NOW)] [$level] $*" >> "$POSTTAX_LOG"
}

pset() {
    local key="$1" value="$2"
    mkdir -p "$(dirname "$POSTTAX_STATUS")"
    if [[ ! -f "$POSTTAX_STATUS" ]]; then
        printf 'key\tvalue\tupdated\n' > "$POSTTAX_STATUS"
    fi
    local tmpf
    tmpf="$(mktemp "$POSTTAX_STATUS.tmp.XXXXXX")"
    grep -v "^$key\t" "$POSTTAX_STATUS" 2>/dev/null > "$tmpf" || true
    printf '%s\t%s\t%s\n' "$key" "$value" "$(NOW)" >> "$tmpf"
    mv "$tmpf" "$POSTTAX_STATUS"
}

fail() {
    plog "FAIL" "$*"
    pset "status" "FAILED"
    pset "failure" "$*"
    exit 1
}

run_stage() {
    local stage="$1"; shift
    plog "stage=$stage starting: $*"
    pset "$stage" "STARTED"
    "$@" >> "$POSTTAX_LOG" 2>&1
    local rc=$?
    if [[ $rc -eq 0 ]]; then
        plog "stage=$stage DONE exit=0"
        pset "$stage" "DONE"
        return 0
    fi
    plog "stage=$stage FAILED exit=$rc"
    pset "$stage" "FAILED exit=$rc"
    return $rc
}

echo "=== posttax watcher started: $(NOW) ===" >> "$POSTTAX_LOG"
pset "status" "WAITING"
pset "started" "$(NOW)"
pset "cohorts" "PRJNA666746,PRJNA822685"
pset "waiting_for" "tax_666746.exit,tax_822685.exit"

while true; do
    if [[ -f "$EXIT_666746" && -f "$EXIT_822685" ]]; then
        break
    fi
    sleep "$WATCH_INTERVAL"
done
pset "status" "EXIT_FILES_PRESENT"
plog "INFO" "both exit files present"

code_666746="$(tr -d '[:space:]' < "$EXIT_666746")"
code_822685="$(tr -d '[:space:]' < "$EXIT_822685")"
plog "INFO" "tax_666746.exit=$code_666746 tax_822685.exit=$code_822685"
pset "tax_666746_exit" "$code_666746"
pset "tax_822685_exit" "$code_822685"

if [[ "$code_666746" != "0" || "$code_822685" != "0" ]]; then
    fail "taxonomy failed (PRJNA666746 exit=$code_666746, PRJNA822685 exit=$code_822685); downstream analysis NOT run"
fi

MISSING=""
for c in $COHORTS; do
    [[ -f "$ROOT/$c/.taxonomy_done" ]] || MISSING="$MISSING $c/.taxonomy_done"
    [[ -f "$ROOT/$c/results/taxonomy.qza" ]] || MISSING="$MISSING $c/results/taxonomy.qza"
    [[ -s "$ROOT/$c/results/genus_count.tsv" ]] || MISSING="$MISSING $c/results/genus_count.tsv"
    [[ -s "$ROOT/$c/results/genus_relabund.tsv" ]] || MISSING="$MISSING $c/results/genus_relabund.tsv"
done
if [[ -n "$MISSING" ]]; then
    fail "taxonomy/genus outputs missing or empty:$MISSING"
fi
pset "taxonomy_verify" "DONE"
plog "INFO" "taxonomy/genus outputs verified for both cohorts"

if [[ ! -f "$MARKER_STATS" ]]; then
    for c in $COHORTS; do
        if ! run_stage "genus_stats:$c" "$PY" "$ROOT/scripts/overnight/genus_stats.py" "$ROOT/$c"; then
            fail "genus_stats failed for $c"
        fi
        [[ -s "$ROOT/$c/results/paired_genus_stats.tsv" ]] || fail "paired_genus_stats.tsv missing for $c"
    done
    : > "$MARKER_STATS"
    pset "genus_stats" "DONE"
else
    pset "genus_stats" "SKIPPED (marker present)"
    plog "INFO" "genus_stats skipped (marker present)"
fi

if [[ ! -f "$MARKER_CROSS" ]]; then
    if run_stage "cross_cohort" "$PY" "$ROOT/scripts/overnight/cross_cohort.py" PRJNA666746 PRJNA822685; then
        CROSS_FILES="candidate_ranking.tsv cohort_summary.tsv cohort_effects.tsv consistent_genera.tsv preliminary_report.md effect_heatmap.png top_candidates.png pipeline_retention.png"
        MISS=""
        for f in $CROSS_FILES; do
            [[ -s "$COMBINED/$f" ]] || MISS="$MISS $f"
        done
        if [[ -n "$MISS" ]]; then
            fail "cross-cohort outputs missing or empty:$MISS"
        fi
        : > "$MARKER_CROSS"
        pset "cross_cohort" "DONE"
    else
        fail "cross_cohort failed"
    fi
else
    pset "cross_cohort" "SKIPPED (marker present)"
    plog "INFO" "cross_cohort skipped (marker present)"
fi

: > "$MARKER_DONE"
pset "status" "COMPLETED"
pset "finished" "$(NOW)"
plog "posttax watcher COMPLETE"
exit 0
