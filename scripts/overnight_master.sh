#!/usr/bin/env bash
# OSCC_16S unattended overnight master pipeline.
# Deterministic, resumable, idempotent, cohort-isolated.
# usage: bash scripts/overnight_master.sh [--resume]
set -uo pipefail

ROOT="$HOME/OSCC_16S"
source "$ROOT/scripts/overnight/common.sh"

COMBINED="$ROOT/combined"
MASTER_LOG="$COMBINED/overnight_master.log"
MASTER_LOCK="$ROOT/.overnight_master.lock"
MASTER_MARKER="$ROOT/.overnight_master_started"

mkdir -p "$COMBINED"
export OVERNIGHT_MASTER_LOG="$MASTER_LOG"

# ---- prevent duplicate master instances ----
exec 8>"$MASTER_LOCK"
if ! flock -n 8; then
    echo "[$(NOW)] [FATAL] another master instance is running (lock: $MASTER_LOCK)" \
        >> "$MASTER_LOG"
    exit 1
fi

echo "=== overnight master started: $(NOW) ===" >> "$MASTER_LOG"
set_status "master_started" "$(NOW)"
set_status "status" "RUNNING"

# SILVA 138.2 reference (prepared once via RESCRIPt; no sklearn training) -----
REFERENCE_SH="$ROOT/scripts/overnight/prepare_reference.sh"
REFERENCE_DONE="$ROOT/data/silva/.reference_done"

run_cohort_stage() {
    # $1 stage name, $2 cohort, $3 script, remaining args
    local stage="$1" cohort="$2" script="$3"; shift 3
    log "stage=$stage cohort=$cohort starting"
    set_status "$stage:$cohort" "STARTED"
    local ok=1
    case "$script" in
        *.py) "$QIIME_ENV/bin/python" "$script" "$@" >> "$MASTER_LOG" 2>&1 || ok=0 ;;
        *)    bash "$script" "$@" >> "$MASTER_LOG" 2>&1 || ok=0 ;;
    esac
    if [[ $ok -eq 1 ]]; then
        set_status "$stage:$cohort" "DONE"
        log "stage=$stage cohort=$cohort DONE"
        return 0
    else
        set_status "$stage:$cohort" "FAILED"
        log "stage=$stage cohort=$cohort FAILED"
        return 1
    fi
}

ALL_COHORTS="PRJNA666746 PRJNA822685 PRJNA813034"

# ---------------------------------------------------------------------------
# 1) SILVA 138.2 reference preparation (needed before taxonomy)
# ---------------------------------------------------------------------------
log "preparing SILVA 138.2 reference"
if bash "$REFERENCE_SH" >> "$MASTER_LOG" 2>&1; then
    set_status "reference" "DONE"
    log "reference preparation DONE"
else
    set_status "reference" "FAILED"
    log "reference preparation FAILED"
fi

# ---------------------------------------------------------------------------
# 2) downloads
# ---------------------------------------------------------------------------
pids=()
for c in PRJNA822685 PRJNA813034; do
    ( run_cohort_stage "download" "$c" "$ROOT/scripts/overnight/download.sh" "$c" ) &
    pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || true; done
DL_OK=0
for c in PRJNA822685 PRJNA813034; do
    if stage_done "$ROOT/$c/.download_done"; then DL_OK=$((DL_OK+1)); fi
done
log "downloads: 2 expected, done=$DL_OK"
set_status "download:summary" "2 expected, done=$DL_OK"

# ---------------------------------------------------------------------------
# 3) preprocessing (parallel across cohorts, flock-protected)
#    Only starts a cohort once its download is complete.
# ---------------------------------------------------------------------------
ppids=()
for c in $ALL_COHORTS; do
    if [[ "$c" == "PRJNA822685" || "$c" == "PRJNA813034" ]]; then
        if ! stage_done "$ROOT/$c/.download_done"; then
            log "skipping preprocess for $c (download not complete)"
            set_status "preprocess:$c" "SKIPPED (download not complete)"
            continue
        fi
    fi
    ( run_cohort_stage "preprocess" "$c" "$ROOT/scripts/overnight/preprocess.sh" "$c" ) &
    ppids+=($!)
done
for p in "${ppids[@]}"; do wait "$p" || true; done

# ---------------------------------------------------------------------------
# 4) DADA2 (sequential, one memory-heavy job at a time)
# ---------------------------------------------------------------------------
for c in $ALL_COHORTS; do
    if ! stage_done "$ROOT/$c/.preprocess_done"; then
        log "skipping DADA2 for $c (preprocess not done)"
        set_status "dada2:$c" "SKIPPED"
        continue
    fi
    run_cohort_stage "dada2" "$c" "$ROOT/scripts/overnight/qiime2_pipeline.sh" "$c"
done

# ---------------------------------------------------------------------------
# 5) taxonomy (requires SILVA reference)
# ---------------------------------------------------------------------------
if stage_done "$REFERENCE_DONE"; then
    for c in $ALL_COHORTS; do
        if ! stage_done "$ROOT/$c/.qiime2_done"; then
            log "skipping taxonomy for $c (dada2 not done)"
            set_status "taxonomy:$c" "SKIPPED"
            continue
        fi
        run_cohort_stage "taxonomy" "$c" "$ROOT/scripts/overnight/taxonomy.sh" "$c"
    done
else
    log "no SILVA reference available; taxonomy skipped for all cohorts"
    for c in $ALL_COHORTS; do set_status "taxonomy:$c" "SKIPPED"; done
fi

# ---------------------------------------------------------------------------
# 7) paired genus statistics
# ---------------------------------------------------------------------------
GENUS_OK=0
for c in $ALL_COHORTS; do
    if ! stage_done "$ROOT/$c/.taxonomy_done"; then
        log "skipping stats for $c (taxonomy not done)"
        set_status "genus_stats:$c" "SKIPPED"
        continue
    fi
    run_cohort_stage "genus_stats" "$c" "$ROOT/scripts/overnight/genus_stats.py" "$c" || true
    if [[ -f "$ROOT/$c/results/paired_genus_stats.tsv" ]]; then
        GENUS_OK=$((GENUS_OK+1))
    fi
done
set_status "genus_stats:summary" "$GENUS_OK cohorts OK"

# ---------------------------------------------------------------------------
# 8) cross-cohort integration (>=2 cohorts)
# ---------------------------------------------------------------------------
GOOD_COHORTS=()
for c in $ALL_COHORTS; do
    [[ -f "$ROOT/$c/results/paired_genus_stats.tsv" ]] && GOOD_COHORTS+=("$c")
done
if [[ ${#GOOD_COHORTS[@]} -ge 2 ]]; then
    log "cross-cohort analysis with: ${GOOD_COHORTS[*]}"
    if "$QIIME_ENV/bin/python" "$ROOT/scripts/overnight/cross_cohort.py" \
            "${GOOD_COHORTS[@]}" >> "$MASTER_LOG" 2>&1; then
        set_status "cross_cohort" "DONE"
    else
        set_status "cross_cohort" "FAILED"
    fi
else
    set_status "cross_cohort" "SKIPPED (<2 cohorts with genus stats)"
    log "cross-cohort skipped: only ${#GOOD_COHORTS[@]} cohort(s) with genus stats"
fi

# ---------------------------------------------------------------------------
# 9) leave-one-cohort-out ML (only if all three genus tables succeed)
# ---------------------------------------------------------------------------
if [[ ${#GOOD_COHORTS[@]} -eq 3 ]]; then
    log "running leave-one-cohort-out ML"
    if "$QIIME_ENV/bin/python" "$ROOT/scripts/overnight/ml_analysis.py" \
            "${GOOD_COHORTS[@]}" >> "$MASTER_LOG" 2>&1; then
        set_status "ml_leave_one_cohort_out" "DONE"
    else
        set_status "ml_leave_one_cohort_out" "FAILED (paired stats unaffected)"
    fi
else
    set_status "ml_leave_one_cohort_out" "SKIPPED (need all 3 genus tables)"
fi

set_status "status" "COMPLETED"
set_status "master_finished" "$(NOW)"
echo "=== overnight master finished: $(NOW) ===" >> "$MASTER_LOG"
