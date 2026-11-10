#!/usr/bin/env bash
# Experimental GAP decomposition. Build separately with `make neon`.
set -euo pipefail
export LC_ALL=C

BIN=./bminv
INPUT=
OUT=
DATASET=sorted
MAX_VECTORS=0
THREADS=10
BATCH=512
PROGRESS=100
REFERENCE=
DISCARD=0
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat <<'HELP'
Run from utils/inv_list after `make neon`:
  bash experiments/gap_xor/run_gov2_gap_xor.sh [options]
  --dataset sorted|unsorted    Gov2 input (default sorted)
  --input FILE                 override raw Gov2 input
  --out DIR                    fresh directory on USB_1TB USB
  --reference-dir DIR          verified control with auto.bv/auto.idx/report.txt
  --max-vectors N              first N vectors; 0 means full (default 0)
  --threads N                  worker threads (default 10)
  --batch-vectors N            pending vectors, 1..512 (default 512)
  --progress N                 vector progress interval (default 100)
  --bin FILE                   bminv binary (default ./bminv)
  --discard-streams            remove generated streams after verification/report
Without a reference, the harness writes, verifies, and decodes both control
and experimental streams. A reference saves space and time; its first N BLOB
offsets supply the control byte count. All streams remain unless discarded.
HELP
}
while (( $# )); do
    case "$1" in
        -h|--help) usage; exit 0;;
        --discard-streams) DISCARD=1; shift;;
        --dataset|--input|--out|--reference-dir|--max-vectors|--threads|--batch-vectors|--progress|--bin)
            (( $# >= 2 )) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --dataset) DATASET=$2;; --input) INPUT=$2;; --out) OUT=$2;;
                --reference-dir) REFERENCE=$2;; --max-vectors) MAX_VECTORS=$2;;
                --threads) THREADS=$2;; --batch-vectors) BATCH=$2;;
                --progress) PROGRESS=$2;; --bin) BIN=$2;;
            esac
            shift 2;;
        *) echo "Unknown option: $1" >&2; exit 2;;
    esac
done
[[ $DATASET == sorted || $DATASET == unsorted ]] || { echo 'Invalid dataset' >&2; exit 2; }
[[ -n $INPUT ]] || INPUT=/Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.$DATASET
for n in "$MAX_VECTORS" "$THREADS" "$BATCH" "$PROGRESS"; do
    [[ $n =~ ^(0|[1-9][0-9]*)$ && ${#n} -le 9 ]] || { echo "Invalid number: $n" >&2; exit 2; }
done
(( THREADS > 0 && THREADS <= 256 && BATCH > 0 && BATCH <= 512 )) || {
    echo 'Vector/thread/batch setting out of range' >&2; exit 2;
}
[[ -x $BIN && -s $INPUT ]] || { echo 'Missing binary or raw input' >&2; exit 2; }
[[ -n $OUT ]] || OUT="/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/gov2-$DATASET-gap-xor-$(date +%Y%m%d-%H%M%S)"
python3 - "$OUT" <<'PY'
import os, subprocess, sys
out = os.path.realpath(sys.argv[1])
if not out.startswith('/Volumes/USB_1TB/'):
    raise SystemExit('Output must be on USB_1TB')
if ' on /Volumes/USB_1TB (' not in subprocess.check_output(['mount'], text=True):
    raise SystemExit('USB_1TB is not mounted')
PY
[[ ! -e $OUT ]] || { echo "Output exists: $OUT" >&2; exit 2; }
if [[ -n $REFERENCE ]]; then
    [[ -s $REFERENCE/auto.bv && -s $REFERENCE/auto.idx && -s $REFERENCE/report.txt ]] || {
        echo 'Reference stream/index/report missing' >&2; exit 2;
    }
    rg -q '^Correctness: PASS' "$REFERENCE/report.txt" || {
        echo 'Reference is not verified' >&2; exit 2;
    }
    [[ $DATASET == sorted ]] || { echo 'This reference mode is for Gov2.sorted' >&2; exit 2; }
    rg -F -x -q "Input: $INPUT" "$REFERENCE/run.txt" || {
        echo 'Reference input differs' >&2; exit 2;
    }
    rg -F -q 'min_gain=16' "$REFERENCE/run.txt" || {
        echo 'Reference BWT threshold differs' >&2; exit 2;
    }
    rg -F -x -q 'mode=auto' "$REFERENCE/run.txt" || {
        echo 'Reference formula mode differs' >&2; exit 2;
    }
    rg -F -q -- '-block-select-step2' "$REFERENCE/run.txt" || {
        echo 'Reference selector differs' >&2; exit 2;
    }
    python3 "$SCRIPT_DIR/report_gap_xor.py" check-reference "$REFERENCE/auto.idx" "$MAX_VECTORS"
fi
mkdir -p -- "$(dirname -- "$OUT")"
mkdir -- "$OUT"
PHASE=initialization
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED phase=%s exit=%s\n" "$PHASE" "$rc" > "$OUT/FAILURE.txt"; printf "FAILED: %s; see %s\n" "$PHASE" "$OUT" >&2; fi; exit "$rc"' EXIT
printf 'Input: %s\nDataset: %s\nReference: %s\nVectors: %s\nThreads: %s\nBatch: %s\n' \
    "$INPUT" "$DATASET" "${REFERENCE:-fresh}" "$MAX_VECTORS" "$THREADS" "$BATCH" > "$OUT/run.txt"

check_space() {
    mount | rg -F -q ' on /Volumes/USB_1TB (' || { echo 'USB disconnected' >&2; return 1; }
    local free_bytes
    free_bytes=$(df -Pk /Volumes/USB_1TB | awk 'NR==2 {print $4 * 1024}')
    python3 - "$free_bytes" <<'PY'
import sys
if float(sys.argv[1]) < 256 * 1024**2:
    raise SystemExit('Less than 256 MiB free on USB')
PY
}
run_phase() {
    PHASE=$1; shift
    check_space
    printf '[%s] %s\n' "$(date +%H:%M:%S)" "$PHASE"
    { printf '%s:' "$PHASE"; printf ' %q' "$@"; printf '\n'; } >> "$OUT/run.txt"
    TIMEFORMAT='%R %U %S'
    if { time "$@" 2>&1 | python3 -u "$SCRIPT_DIR/filter_log.py" | tee "$OUT/$PHASE.log" | python3 "$SCRIPT_DIR/../current/progress_console.py"; } 2> "$OUT/$PHASE.time"; then
        :
    else
        rc=$?; tail -n 15 "$OUT/$PHASE.log" >&2; return "$rc"
    fi
}
run_one() {
    local mode=$1
    local stream="$OUT/$mode.bv" index="$OUT/$mode.idx"
    local args=(-t -s -level 6 -u32in "$INPUT" -bvout "$stream"
        -block-select-step2 -periodic-step1-interleaved 4 -bwt-min-gain 16
        -threads "$THREADS" -batch-vectors "$BATCH" -bv-index-out "$index")
    [[ $mode != candidate ]] || args+=(-gap-xor)
    local read_args=(-bvin "$stream" -bv-index-in "$index"
        -threads "$THREADS" -batch-vectors "$BATCH")
    local verify_args=(-t -s -verify -u32in "$INPUT" "${read_args[@]}")
    local decode_args=(-t -s -decode "${read_args[@]}")
    if (( MAX_VECTORS )); then
        args+=(-max-vectors "$MAX_VECTORS")
        verify_args+=(-max-vectors "$MAX_VECTORS")
    fi
    if (( PROGRESS )); then
        args+=(-progress-vectors "$PROGRESS")
        verify_args+=(-progress-vectors "$PROGRESS")
        decode_args+=(-progress-vectors "$PROGRESS")
    fi
    run_phase "$mode.compress" "$BIN" "${args[@]}"
    run_phase "$mode.verify" "$BIN" "${verify_args[@]}"
    run_phase "$mode.decode" "$BIN" "${decode_args[@]}"
}
if [[ -z $REFERENCE ]]; then run_one control; fi
run_one candidate
PHASE=report
python3 "$SCRIPT_DIR/report_gap_xor.py" report "$OUT" "$MAX_VECTORS" "$REFERENCE" | tee "$OUT/report.txt"
if (( DISCARD )); then
    rm -f -- "$OUT"/*.bv "$OUT"/*.idx
    printf 'Generated streams discarded after successful verify/decode/report.\n' | tee -a "$OUT/report.txt"
fi
printf '\nReports and logs: %s\n' "$OUT"
