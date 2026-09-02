#!/usr/bin/env bash
# Primer-aware preprocessing for one cohort (parallel, resumable):
#   - deterministic cap raw pairs at 150000/sample (seed fixed)
#   - FR/RF primer recognition + trimming (cutadapt e=0.15, no indels)
#   - orientation normalization (final R1 = biological forward)
#   - optional 3' 10bp trim (PRJNA813034)
#   - cap valid pairs at 100000/sample
#   - cohort-level validity gate (>=70% valid primer pairs) for gated cohorts
# usage: preprocess.sh <COHORT>
set -euo pipefail
source "$(dirname "$0")/common.sh"

COHORT="$1"
BASE="$(cohort_dir "$COHORT")"
META="$BASE/metadata"
OUT="$BASE/processed"
LOGS="$BASE/logs"
WORK="$OUT/.work"
LOCK="$BASE/.preprocess.lock"
MARKER="$BASE/.preprocess_done"
DIR="$(dirname "$0")"

GATE_VALID="no"
case "$COHORT" in
    PRJNA666746) GATE_VALID="no" ;;
    PRJNA822685) GATE_VALID="yes" ;;
    # 70% global primer-retention gate withdrawn for PRJNA813034:
    # only strict confidently-assigned FR/RF processed reads are used.
    PRJNA813034) GATE_VALID="no" ;;
    *) echo "unknown cohort $COHORT" >&2; exit 1 ;;
esac

PARALLELISM="${PREPROCESS_PARALLELISM:-3}"

mkdir -p "$OUT" "$WORK" "$LOGS"
RETENTION="$OUT/retention.tsv"
VALID_SUMMARY="$OUT/orientation_summary.tsv"

acquire_lock "$LOCK" || { echo "cannot lock preprocess for $COHORT"; exit 1; }

if stage_done "$MARKER"; then
    echo "preprocess already complete for $COHORT"
    exit 0
fi

if [[ ! -s "$RETENTION" ]]; then
    printf 'sample\traw_pairs\tpretrim_pairs\tFR_pairs\tRF_pairs\tvalid_primer_pairs\tvalid_fraction\tfinal_retained_pairs\n' > "$RETENTION"
fi

# fast resume: if every manifest sample already has a retention row AND its
# final processed FASTQs exist, mark the cohort complete without redoing work.
RESUME_OK=1
if [[ -s "$RETENTION" ]]; then
    n_manifest=$(sed 's/\r$//' "$META/manifest.tsv" | tail -n +2 | wc -l)
    n_ret=$(tail -n +2 "$RETENTION" | wc -l)
    if [[ "$n_manifest" -ne "$n_ret" ]]; then
        RESUME_OK=0
    else
        while IFS=$'\t' read -r sample raw pre fr rf valid frac final; do
            if [[ ! -s "$OUT/${sample}_1.fastq.gz" || ! -s "$OUT/${sample}_2.fastq.gz" ]]; then
                RESUME_OK=0
                break
            fi
        done < <(tail -n +2 "$RETENTION")
    fi
else
    RESUME_OK=0
fi
if [[ "$RESUME_OK" -eq 1 ]]; then
    echo "$(NOW) [$COHORT] retention complete for all $n_manifest samples; marking preprocess done" \
        | tee -a "$LOGS/preprocess.log"
    mark_done "$MARKER"
    echo "$(NOW) [$COHORT] preprocessing COMPLETE (fast resume)" \
        | tee -a "$LOGS/preprocess.log"
    exit 0
fi

sed 's/\r$//' "$META/manifest.tsv" | tail -n +2 \
    | xargs -d '\n' -P "$PARALLELISM" -I'{}' \
        bash "$DIR/preprocess_sample.sh" "$COHORT" '{}'

rc=$?
if [[ $rc -ne 0 ]]; then
    echo "$(NOW) [$COHORT] preprocess worker(s) failed" | tee -a "$LOGS/preprocess.log"
    mark_fail "$BASE/.preprocess_failed"
    exit 1
fi

# ---- cohort-level validity gate ----
echo "$(NOW) [$COHORT] orientation summary:" | tee -a "$LOGS/preprocess.log"
overall=$(awk -F'\t' 'NR>1{sv+=$6; sp+=$3} END{if(sp>0) printf "%.4f", sv/sp; else print "0"}' "$RETENTION")
printf 'sum_valid\t sum_pretrim\t overall_fraction\n' > "$VALID_SUMMARY"
awk -F'\t' -v o="$overall" 'NR>1{sv+=$6; sp+=$3} END{printf "%d\t%d\t%s\n", sv, sp, o}' "$RETENTION" >> "$VALID_SUMMARY"
tee -a "$LOGS/preprocess.log" < "$VALID_SUMMARY" > /dev/null

if [[ "$GATE_VALID" == "yes" ]]; then
    if awk -v o="$overall" 'BEGIN{ exit !(o < 0.70) }'; then
        echo "$(NOW) [$COHORT] FAIL: valid primer fraction $overall < 0.70. Stopping only this cohort." \
            | tee -a "$LOGS/preprocess.log"
        mark_fail "$BASE/.preprocess_failed"
        exit 1
    fi
fi

mark_done "$MARKER"
echo "$(NOW) [$COHORT] preprocessing COMPLETE (overall valid fraction $overall)" \
    | tee -a "$LOGS/preprocess.log"
