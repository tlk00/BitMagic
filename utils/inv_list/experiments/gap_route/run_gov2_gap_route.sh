#!/usr/bin/env bash
# Build separately: cd utils/inv_list && make neon
set -euo pipefail
export LC_ALL=C

BIN=./bminv
INPUT=/Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.sorted
OUT=
THREADS=10
BATCH=512
PROGRESS=1000
MAX_VECTORS=0
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)

usage() {
    cat <<'HELP'
Run from utils/inv_list after make neon:
  bash experiments/gap_route/run_gov2_gap_route.sh [options]
  --out DIR          fresh directory on USB_1TB (required)
  --input FILE       raw Gov2.sorted (default /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.sorted)
  --bin FILE         bminv (default ./bminv)
  --threads N        compression and read workers (default 10)
  --batch-vectors N  pending vectors, 1..512 (default 512)
  --progress N       progress interval (default 1000)
  --max-vectors N    first N vectors; omit or use 0 for full collection
HELP
}
while (( $# )); do
    case "$1" in
        -h|--help) usage; exit 0;;
        --out|--input|--bin|--threads|--batch-vectors|--progress|--max-vectors)
            (( $# >= 2 )) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --out) OUT=$2;; --input) INPUT=$2;; --bin) BIN=$2;;
                --threads) THREADS=$2;; --batch-vectors) BATCH=$2;;
                --progress) PROGRESS=$2;; --max-vectors) MAX_VECTORS=$2;;
            esac
            shift 2;;
        *) echo "Unknown option: $1" >&2; exit 2;;
    esac
done
[[ -n $OUT && -x $BIN && -s $INPUT ]] || { echo 'Require --out, built bminv, and input' >&2; exit 2; }
for n in "$THREADS" "$BATCH" "$PROGRESS" "$MAX_VECTORS"; do
    [[ $n =~ ^(0|[1-9][0-9]*)$ && ${#n} -le 9 ]] || { echo "Invalid number: $n" >&2; exit 2; }
done
(( THREADS > 0 && THREADS <= 256 && BATCH > 0 && BATCH <= 512 )) || { echo 'Thread/batch setting out of range' >&2; exit 2; }
python3 - "$OUT" <<'PY'
import os, subprocess, sys
out = os.path.realpath(sys.argv[1])
if not out.startswith('/Volumes/USB_1TB/'):
    raise SystemExit('Output must be on USB_1TB')
if ' on /Volumes/USB_1TB (' not in subprocess.check_output(['mount'], text=True):
    raise SystemExit('USB_1TB is not mounted')
PY
[[ ! -e $OUT ]] || { echo "Output already exists: $OUT" >&2; exit 2; }
mkdir -p -- "$(dirname -- "$OUT")"
mkdir -- "$OUT"
PHASE=initialization
trap 'rc=$?; if (( rc != 0 )); then printf "FAILED phase=%s exit=%s\n" "$PHASE" "$rc" > "$OUT/FAILURE.txt"; printf "FAILED: %s; see %s\n" "$PHASE" "$OUT" >&2; fi; exit "$rc"' EXIT
printf 'Input: %s\nVectors: %s\nThreads: %s\n' "$INPUT" "$MAX_VECTORS" "$THREADS" > "$OUT/run.txt"

check_space() {
    mount | rg -F -q ' on /Volumes/USB_1TB (' || { echo 'USB disconnected' >&2; return 1; }
    local free_kib
    free_kib=$(df -Pk /Volumes/USB_1TB | awk 'NR==2 {print $4}')
    local needed_kib=262144
    [[ $PHASE != compress ]] || needed_kib=4194304
    (( free_kib >= needed_kib )) || { echo "Need at least $((needed_kib / 1024)) MiB free on USB" >&2; return 1; }
}
run_phase() {
    PHASE=$1; shift
    check_space
    printf '[%s] %s\n' "$(date +%H:%M:%S)" "$PHASE"
    { printf '%s:' "$PHASE"; printf ' %q' "$@"; printf '\n'; } >> "$OUT/run.txt"
    TIMEFORMAT='%R %U %S'
    if { time "$@" 2>&1 | tee "$OUT/$PHASE.log" | python3 "$SCRIPT_DIR/../current/progress_console.py"; } 2> "$OUT/$PHASE.time"; then
        :
    else
        rc=$?; tail -n 15 "$OUT/$PHASE.log" >&2; return "$rc"
    fi
}

stream="$OUT/survey.bv"
index="$OUT/survey.idx"
common=(-t -s)
compress=("${common[@]}" -level 6 -u32in "$INPUT" -bvout "$stream"
    -block-select-step2 -periodic-step1-interleaved 4 -bwt-min-gain 16
    -gap-route-stats -threads "$THREADS" -batch-vectors "$BATCH" -bv-index-out "$index")
read_args=(-bvin "$stream" -bv-index-in "$index" -threads "$THREADS" -batch-vectors "$BATCH")
verify=("${common[@]}" -verify -u32in "$INPUT" "${read_args[@]}")
decode=("${common[@]}" -decode "${read_args[@]}")
if (( MAX_VECTORS )); then
    compress+=(-max-vectors "$MAX_VECTORS")
    verify+=(-max-vectors "$MAX_VECTORS")
fi
if (( PROGRESS )); then
    compress+=(-progress-vectors "$PROGRESS")
    verify+=(-progress-vectors "$PROGRESS")
    decode+=(-progress-vectors "$PROGRESS")
fi
run_phase compress "$BIN" "${compress[@]}"
run_phase verify "$BIN" "${verify[@]}"
run_phase decode "$BIN" "${decode[@]}"
PHASE=report
python3 "$SCRIPT_DIR/report_gap_route.py" "$OUT" | tee "$OUT/report.txt"
printf '\nReports and logs: %s\n' "$OUT"
