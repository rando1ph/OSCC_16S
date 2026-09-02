#!/usr/bin/env bash
# Shared helpers for the OSCC_16S overnight pipeline.
set -u

ROOT="$HOME/OSCC_16S"
QIIME_ENV="$HOME/miniconda3/envs/rachis-qiime2-2026.7"
export R_HOME="$QIIME_ENV/lib/R"
CUTADAPT="$HOME/.local/bin/cutadapt"
SEED="20260901"

NOW() { date '+%Y-%m-%d %H:%M:%S'; }

OVERNIGHT_MASTER_LOG="${OVERNIGHT_MASTER_LOG:-$ROOT/combined/overnight_master.log}"

log() {
    local level="$1"; shift
    local msg="$*"
    local stamp
    stamp="$(NOW)"
    echo "[$stamp] [$level] $msg" >> "$OVERNIGHT_MASTER_LOG"
    if [[ "${VERBOSE:-0}" == "1" ]]; then
        echo "[$stamp] [$level] $msg"
    fi
}

# Acquire an exclusive flock. Blocks until available.
# usage: acquire_lock <lockfile>
acquire_lock() {
    local lockfile="$1"
    exec 9>"$lockfile"
    flock -x 9 || return 1
    return 0
}

# Run a command while logging to a cohort log file. First arg is the logfile.
log_run() {
    local lf="$1"; shift
    mkdir -p "$(dirname "$lf")"
    {
        echo "=== command: $* ==="
        "$@"
        echo "=== exit=$? ==="
    } >> "$lf" 2>&1
    return $?
}

cohort_dir() {
    echo "$ROOT/$1"
}

# stage marker helpers -------------------------------------------------------
# Stages are idempotent: a marker file records completion.
stage_done() {  # <marker>
    [[ -f "$1" ]]
}

mark_done() {  # <marker>
    mkdir -p "$(dirname "$1")"
    : > "$1"
}

mark_fail() {  # <marker>
    mkdir -p "$(dirname "$1")"
    echo "FAILED $(NOW)" > "$1"
}

# Update combined/pipeline_status.tsv key=value line.
# usage: set_status <key> <value>
# NOTE: uses mktemp (not $$) so concurrent subshell calls cannot collide.
STATUS_FILE="$ROOT/combined/pipeline_status.tsv"
set_status() {
    local key="$1" value="$2"
    mkdir -p "$ROOT/combined"
    if [[ ! -f "$STATUS_FILE" ]]; then
        printf 'key\tvalue\tupdated\n' > "$STATUS_FILE"
    fi
    local updated
    updated="$(NOW)"
    local tmpf
    tmpf="$(mktemp "$STATUS_FILE.tmp.XXXXXX")"
    # remove existing key then append
    grep -v "^$key\t" "$STATUS_FILE" 2>/dev/null > "$tmpf" || true
    printf '%s\t%s\t%s\n' "$key" "$value" "$updated" >> "$tmpf"
    mv "$tmpf" "$STATUS_FILE"
}

# Count gzipped fastq pairs in a pair of files.
count_pairs() {
    local r1="$1" r2="$2"
    python3 - "$r1" "$r2" <<'PYEOF'
import gzip, sys
def count(path):
    n = 0
    with gzip.open(path, 'rt') as fh:
        while True:
            line = fh.readline()
            if not line:
                break
            n += 1
            for _ in range(3):
                if not fh.readline():
                    break
    return n // 4
print(count(sys.argv[1]))
PYEOF
}
