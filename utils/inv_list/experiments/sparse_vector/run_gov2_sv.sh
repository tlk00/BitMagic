#!/usr/bin/env bash
# Complete Gov2 BV-versus-SV comparison. Build separately: make neon svexp-neon.
set -euo pipefail
export LC_ALL=C

INPUT=/Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.sorted
OUT=
BV_REFERENCE_DIR=
MAX_VECTORS=0
THREADS=10
BATCH=512
PROGRESS=1000
BIN=./svexp
BV_BIN=./bminv

usage() {
    cat <<'HELP'
Usage: bash experiments/sparse_vector/run_gov2_sv.sh [options]
Run from utils/inv_list after: make neon svexp-neon
  --input FILE             Gov2.sorted or Gov2.unsorted raw collection
  --out DIR                fresh directory on mounted USB_1TB USB
  --bv-reference-dir DIR   reuse a verified BV auto.bv/auto.idx (same input/count)
  --max-vectors N          0 means complete collection (default 0)
  --threads N              worker threads (default 10)
  --batch N                pending vectors, 1..512 (default 512)
  --progress N             progress interval, 0 disables (default 1000)
  --bin FILE               svexp binary (default ./svexp)
  --bv-bin FILE            bminv binary (default ./bminv)
The run writes six SV streams (raw IDs / reduced gaps, each plain / XOR /
BWT+chain selection), verifies and decodes each, and compares a BV reference.
HELP
}
while (( $# )); do
    case "$1" in
        -h|--help) usage; exit 0 ;;
        --input|--out|--bv-reference-dir|--max-vectors|--threads|--batch|--progress|--bin|--bv-bin)
            (( $# >= 2 )) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --input) INPUT=$2;; --out) OUT=$2;;
                --bv-reference-dir) BV_REFERENCE_DIR=$2;;
                --max-vectors) MAX_VECTORS=$2;; --threads) THREADS=$2;;
                --batch) BATCH=$2;; --progress) PROGRESS=$2;;
                --bin) BIN=$2;; --bv-bin) BV_BIN=$2;;
            esac
            shift 2 ;;
        *) echo "Unknown option: $1" >&2; exit 2 ;;
    esac
done
for n in "$MAX_VECTORS" "$THREADS" "$BATCH" "$PROGRESS"; do
    [[ $n =~ ^(0|[1-9][0-9]*)$ && ${#n} -le 9 ]] || {
        echo "Invalid number: $n" >&2; exit 2;
    }
done
(( THREADS >= 1 && THREADS <= 256 && BATCH >= 1 && BATCH <= 512 )) || {
    echo 'Thread or batch setting out of range' >&2; exit 2;
}
[[ -x "$BIN" && -x "$BV_BIN" && -s "$INPUT" ]] || {
    echo 'Missing binaries or Gov2 input; build with make neon svexp-neon' >&2; exit 2;
}
if [[ -z $OUT ]]; then
    case "${INPUT##*/}" in
        gov2.sorted) DATASET=gov2-sorted;;
        gov2.unsorted) DATASET=gov2-unsorted;;
        *) DATASET=collection;;
    esac
    OUT="/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/${DATASET}-svexp-$(date +%Y%m%d-%H%M%S)"
fi
python3 - "$OUT" <<'CHECK_MOUNT'
import os, subprocess, sys
path = os.path.realpath(sys.argv[1])
if not path.startswith('/Volumes/USB_1TB/'):
    raise SystemExit('Output must be under /Volumes/USB_1TB')
if ' on /Volumes/USB_1TB (' not in subprocess.check_output(['mount'], text=True):
    raise SystemExit('USB_1TB USB drive is not mounted')
CHECK_MOUNT
[[ ! -e $OUT ]] || { echo "Refusing existing output directory: $OUT" >&2; exit 2; }
if [[ -n $BV_REFERENCE_DIR ]]; then
    [[ -s $BV_REFERENCE_DIR/auto.bv && -s $BV_REFERENCE_DIR/auto.idx &&
       -s $BV_REFERENCE_DIR/report.txt && -s $BV_REFERENCE_DIR/run.txt ]] || {
        echo 'Missing verified BV reference artifacts' >&2; exit 2;
    }
    rg -q '^Correctness: PASS' "$BV_REFERENCE_DIR/report.txt" || {
        echo 'BV reference report did not pass' >&2; exit 2;
    }
    rg -F -x -q "Input: $INPUT" "$BV_REFERENCE_DIR/run.txt" || {
        echo 'BV reference input differs' >&2; exit 2;
    }
    rg -F -q "level=6 vectors=$MAX_VECTORS " "$BV_REFERENCE_DIR/run.txt" || {
        echo 'BV reference count or level differs' >&2; exit 2;
    }
    rg -F -q 'min_gain=16' "$BV_REFERENCE_DIR/run.txt" || {
        echo 'BV reference threshold differs' >&2; exit 2;
    }
    rg -F -x -q 'mode=auto' "$BV_REFERENCE_DIR/run.txt" || {
        echo 'BV reference must use adaptive gamma/fixed formula choice' >&2; exit 2;
    }
    rg -F -q -- '-block-select-step2' "$BV_REFERENCE_DIR/run.txt" || {
        echo 'BV reference must use the independent step-2 selector' >&2; exit 2;
    }
fi
mkdir -p -- "$(dirname -- "$OUT")"
mkdir -- "$OUT"
PHASE=initialization
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED phase=%s exit=%s\n" "$PHASE" "$rc" > "$OUT/FAILURE.txt"; printf "FAILED: %s; see %s\n" "$PHASE" "$OUT" >&2; fi; exit "$rc"' EXIT
printf 'Input: %s\nBV reference: %s\nvectors=%s threads=%s batch=%s progress=%s\n' \
    "$INPUT" "${BV_REFERENCE_DIR:-fresh}" "$MAX_VECTORS" "$THREADS" "$BATCH" "$PROGRESS" > "$OUT/run.txt"

check_space() {
    mount | rg -F -q ' on /Volumes/USB_1TB (' || {
        echo 'USB_1TB USB drive disconnected' >&2; return 1;
    }
    local free_bytes
    free_bytes=$(df -Pk /Volumes/USB_1TB | awk 'NR==2 {print $4 * 1024}')
    python3 - "$free_bytes" <<'CHECK_SPACE'
import sys
if float(sys.argv[1]) < 4 * 1024**3:
    raise SystemExit('Less than 4 GiB USB space remains; stop before next stream')
CHECK_SPACE
}

check_mount() {
    mount | rg -F -q ' on /Volumes/USB_1TB (' || {
        echo 'USB_1TB USB drive disconnected' >&2; return 1;
    }
}

run_phase() {
    PHASE=$1; shift
    case "$PHASE" in
        *.compress) check_space;;
        *) check_mount;;
    esac
    printf '[%s] %s\n' "$(date +%H:%M:%S)" "$PHASE"
    { printf '%s:' "$PHASE"; printf ' %q' "$@"; printf '\n'; } >> "$OUT/run.txt"
    TIMEFORMAT='%R %U %S'
    if { time "$@" 2>&1 | tee "$OUT/$PHASE.log" | \
         python3 experiments/current/progress_console.py; } 2> "$OUT/$PHASE.time"; then
        :
    else
        rc=$?; tail -n 12 "$OUT/$PHASE.log" >&2; return "$rc"
    fi
}

if [[ -z $BV_REFERENCE_DIR ]]; then
    BV_STREAM=$OUT/bv.auto.bv
    BV_INDEX=$OUT/bv.auto.idx
    bv_args=(-t -s -level 6 -u32in "$INPUT" -bvout "$BV_STREAM"
             -block-select-step2 -periodic-step1-interleaved 4 -bwt-min-gain 16
             -threads "$THREADS" -batch-vectors "$BATCH" -bv-index-out "$BV_INDEX")
    bv_verify=(-t -s -verify -u32in "$INPUT" -bvin "$BV_STREAM"
               -bv-index-in "$BV_INDEX" -threads "$THREADS" -batch-vectors "$BATCH")
    bv_decode=(-t -s -decode -bvin "$BV_STREAM" -bv-index-in "$BV_INDEX"
               -threads "$THREADS" -batch-vectors "$BATCH")
    if (( MAX_VECTORS )); then
        bv_args+=(-max-vectors "$MAX_VECTORS")
        bv_verify+=(-max-vectors "$MAX_VECTORS")
    fi
    if (( PROGRESS )); then
        bv_args+=(-progress-vectors "$PROGRESS")
        bv_verify+=(-progress-vectors "$PROGRESS")
        bv_decode+=(-progress-vectors "$PROGRESS")
    fi
    run_phase bv.compress "$BV_BIN" "${bv_args[@]}"
    run_phase bv.verify "$BV_BIN" "${bv_verify[@]}"
    run_phase bv.decode "$BV_BIN" "${bv_decode[@]}"
else
    BV_STREAM=$BV_REFERENCE_DIR/auto.bv
    BV_INDEX=$BV_REFERENCE_DIR/auto.idx
fi

for repr in raw gaps; do
    for codec in plain xor select; do
        name=$repr.$codec
        stream=$OUT/$name.sv
        index=$OUT/$name.idx
        common=(--stream "$stream" --index "$index" --representation "$repr"
                --codec "$codec" --threads "$THREADS" --batch "$BATCH"
                --progress "$PROGRESS")
        run_phase "$name.compress" "$BIN" encode --input "$INPUT" "${common[@]}" --max-vectors "$MAX_VECTORS"
        run_phase "$name.verify" "$BIN" verify --input "$INPUT" "${common[@]}" --max-vectors "$MAX_VECTORS"
        run_phase "$name.decode" "$BIN" decode "${common[@]}"
    done
done

PHASE=report
python3 experiments/sparse_vector/report_gov2_sv.py "$OUT" "$BV_STREAM" \
    "$INPUT" "$MAX_VECTORS" > "$OUT/report.txt"
cat "$OUT/report.txt"
printf '\nReports and logs: %s\n' "$OUT"
