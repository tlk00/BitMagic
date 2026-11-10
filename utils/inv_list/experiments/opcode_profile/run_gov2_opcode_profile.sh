#!/usr/bin/env bash
# Build separately with make neon in utils/inv_list.
set -euo pipefail
export LC_ALL=C
BIN=./bminv
INPUT_DIR=/Volumes/USB_1TB/int_set_bench/Gov2Flat
OUT=
DATASETS=sorted
THREADS=10
MAX_VECTORS=0
PROGRESS=1000
KEEP=0
LEVEL=
BWT_GAIN=
CORRECTION_SELECT=0
SBLOCK_LARGE=0
INVERSE_LIMIT=196607
DIAGNOSTICS=0
BLOCK_ENDPOINT_SELECT=0
BWT_SHARED=0
SBLOCK_WAVE_ENCODE=0
SBLOCK_WAVE_ADAPTIVE=0
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
usage() {
    cat <<'HELP'
Usage: bash experiments/opcode_profile/run_gov2_opcode_profile.sh --out DIR [options]
  --level N            use the library preset (0..10), without codec flags
  --bwt-min-gain N     override BWT threshold for preset 9/10
  --datasets both|sorted|unsorted (default sorted)
  --input-dir DIR       contains gov2.sorted and gov2.unsorted
  --bin FILE            default ./bminv; build separately
  --threads N          default 10; use 1 for controlled timing
  --max-vectors N      default 0 (full collection)
  --progress N         every N vectors or 30 seconds (default N=1000)
  --diagnostics         print worker timing diagnostics
  --keep-streams        retain verified .bv files (default: delete only this
                        run's .bv file after each mode report passes)
  --sblock-large         optional opcode 98 for positive/inverse groups
  --sblock-large-inverse-limit N  maximum zero-list count (default 196607)
  --sblock-correction-select compare corrected/uncorrected opcode-68 records
  --block-endpoint-select emit the smaller existing opcode-65 endpoint layout
  --bwt-shared            original-polarity shared BWT for 72/96 selection
  --sblock-wave-adaptive emit experimental adaptive-width opcode 97 using dry scores
  --sblock-wave         emit opcode 90 when it beats ordinary/adaptive superblock
Output must be a fresh directory on mounted /Volumes/USB_1TB.
One stream with current adaptive superblocks and an opcode-size profile.
Grouped flags and complete-record selection against ordinary BIC.
Without --level, uses the historical explicit configuration (16-byte BWT cutoff).
With --level, uses the library preset and optional BWT gain override.
Reports, indices and logs are always retained. Timings are end-to-end measured;
The stream is emitted, verified and independently decoded.
HELP
}
while (( $# )); do
    case "$1" in
        -h|--help) usage; exit 0;;
        --keep-streams) KEEP=1; shift;;
        --sblock-large) SBLOCK_LARGE=1; shift;;
        --sblock-correction-select) CORRECTION_SELECT=1; shift;;
        --diagnostics) DIAGNOSTICS=1; shift;;
        --block-endpoint-select) BLOCK_ENDPOINT_SELECT=1; shift;;
        --bwt-shared) BWT_SHARED=1; shift;;
        --sblock-wave-adaptive) SBLOCK_WAVE_ADAPTIVE=1; SBLOCK_WAVE_ENCODE=1; shift;;
        --sblock-wave) SBLOCK_WAVE_ENCODE=1; shift;;
        --out|--bin|--input-dir|--datasets|--threads|--max-vectors|--progress|--sblock-large-inverse-limit|--level|--bwt-min-gain)
            (( $# >= 2 )) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --level) LEVEL=$2;; --bwt-min-gain) BWT_GAIN=$2;;
                --out) OUT=$2;; --bin) BIN=$2;; --input-dir) INPUT_DIR=$2;;
                --datasets) DATASETS=$2;; --threads) THREADS=$2;;
                --max-vectors) MAX_VECTORS=$2;; --progress) PROGRESS=$2;; --sblock-large-inverse-limit) INVERSE_LIMIT=$2;;
            esac
            shift 2;;
        *) echo "Unknown argument: $1" >&2; exit 2;;
    esac
done
[[ -n $OUT && -x $BIN ]] || { echo 'Require --out and built binary' >&2; exit 2; }
case "$DATASETS" in both) sets=(sorted unsorted);; sorted|unsorted) sets=("$DATASETS");; *) echo 'Bad --datasets' >&2; exit 2;; esac
for n in "$THREADS" "$MAX_VECTORS" "$PROGRESS" "$INVERSE_LIMIT"; do
    [[ $n =~ ^(0|[1-9][0-9]*)$ && ${#n} -le 9 ]] || { echo "Bad number: $n" >&2; exit 2; }
done
(( THREADS > 0 && THREADS <= 256 && PROGRESS > 0 )) || { echo 'Bad thread/progress/coder setting' >&2; exit 2; }
if [[ -n $LEVEL ]]; then
    [[ $LEVEL =~ ^([0-9]|10)$ ]] || { echo 'Level must be 0..10' >&2; exit 2; }
    (( ! CORRECTION_SELECT && ! SBLOCK_LARGE && ! BLOCK_ENDPOINT_SELECT && ! BWT_SHARED && ! SBLOCK_WAVE_ENCODE )) || {
        echo 'Use --level without legacy codec-selection flags' >&2; exit 2;
    }
fi
if [[ -n $BWT_GAIN ]]; then
    [[ $LEVEL =~ ^(9|10)$ && $BWT_GAIN =~ ^[0-9]+$ && ${#BWT_GAIN} -le 5 ]] && (( BWT_GAIN <= 65535 )) || {
        echo '--bwt-min-gain requires level 9/10 and gain 0..65535' >&2; exit 2;
    }
fi
for dataset in "${sets[@]}"; do
    [[ -s $INPUT_DIR/gov2.$dataset ]] || { echo "Missing $INPUT_DIR/gov2.$dataset" >&2; exit 2; }
done
python3 - "$OUT" <<'PY'
import os, subprocess, sys
if not os.path.realpath(sys.argv[1]).startswith('/Volumes/USB_1TB/'):
    raise SystemExit('Output must be under /Volumes/USB_1TB')
if ' on /Volumes/USB_1TB (' not in subprocess.check_output(['mount'], text=True):
    raise SystemExit('USB_1TB is not mounted')
PY
[[ ! -e $OUT ]] || { echo "Output already exists: $OUT" >&2; exit 2; }
mkdir -p -- "$(dirname -- "$OUT")"
mkdir -- "$OUT"
PHASE=initialization
trap 'rc=$?; if (( rc )); then printf "FAILED phase=%s exit=%s\n" "$PHASE" "$rc" > "$OUT/FAILURE.txt"; printf "FAILED: %s; see %s\n" "$PHASE" "$OUT" >&2; fi; exit "$rc"' EXIT
printf 'Datasets: %s\nThreads: %s\nVectors: %s\nKeep streams: %s\nConfiguration: current best selection plus top-level opcode profile\n' "$DATASETS" "$THREADS" "$MAX_VECTORS" "$KEEP" > "$OUT/run.txt"
run_phase() {
    local dir=$1 phase=$2
    shift 2
    PHASE="$dataset.$mode.$phase"
    mount | rg -F -q ' on /Volumes/USB_1TB (' || { echo 'USB disconnected' >&2; return 1; }
    if [[ $phase == compress ]]; then
        local free_kib
        free_kib=$(df -Pk /Volumes/USB_1TB | awk 'NR==2 {print $4}')
        (( free_kib >= 4194304 )) || { echo 'Need at least 4 GiB free before compression' >&2; return 1; }
    fi
    printf '[%s] %s\n' "$(date +%H:%M:%S)" "$PHASE"
    { printf '%s:' "$PHASE"; printf ' %q' "$@"; printf '\n'; } >> "$OUT/run.txt"
    TIMEFORMAT='%R %U %S'
    if { time "$@" 2>&1 | tee "$dir/$phase.log" | python3 "$SCRIPT_DIR/../current/progress_console.py" --hide-opcode-profile; } 2> "$dir/$phase.time"; then
        :
    else
        rc=$?; tail -n 15 "$dir/$phase.log" >&2; return "$rc"
    fi
}
for dataset in "${sets[@]}"; do
    for mode in profile; do
        dir="$OUT/$dataset/$mode"
        mkdir -p -- "$dir"
        input="$INPUT_DIR/gov2.$dataset"
        compress=(-t -s -level 6 -u32in "$input" -bvout "$dir/candidate.bv"
            -block-select-step2 -periodic-step1-interleaved 4 -bwt-min-gain 16
            -gap-xor -sblock-gap-tolerance 100 -sblock-adaptive -profile-opcodes
            -threads "$THREADS" -batch-vectors 512 -bv-index-out "$dir/candidate.idx"
            -progress-vectors "$PROGRESS")
        if [[ -n $LEVEL ]]; then
            # Preset-only path: library API supplies every codec setting.
            compress=(-t -s -level "$LEVEL" -u32in "$input" -bvout "$dir/candidate.bv"
                -profile-opcodes -threads "$THREADS" -batch-vectors 512
                -bv-index-out "$dir/candidate.idx" -progress-vectors "$PROGRESS")
            if [[ -n $BWT_GAIN ]]; then compress+=(-bwt-min-gain "$BWT_GAIN"); fi
            if (( LEVEL >= 8 )); then compress+=(-sblock-large-inverse-limit "$INVERSE_LIMIT"); fi
        fi
        if (( SBLOCK_LARGE )); then compress+=(-sblock-large -sblock-large-inverse-limit "$INVERSE_LIMIT"); fi
        if (( CORRECTION_SELECT )); then compress+=(-sblock-correction-select); fi
        if (( BLOCK_ENDPOINT_SELECT )); then compress+=(-block-endpoint-select); fi
        if (( DIAGNOSTICS )); then compress+=(-compression-diagnostics); fi
        if (( BWT_SHARED )); then compress+=(-bwt-shared); fi
        if (( SBLOCK_WAVE_ENCODE )); then compress+=(-sblock-wave); fi
        if (( SBLOCK_WAVE_ADAPTIVE )); then compress+=(-sblock-wave-adaptive); fi
        read_args=(-t -s -bvin "$dir/candidate.bv" -bv-index-in "$dir/candidate.idx"
            -threads "$THREADS" -batch-vectors 512 -progress-vectors "$PROGRESS")
        verify=(-verify -u32in "$input" "${read_args[@]}")
        if (( MAX_VECTORS )); then
            compress+=(-max-vectors "$MAX_VECTORS")
            verify+=(-max-vectors "$MAX_VECTORS")
        fi
        run_phase "$dir" compress "$BIN" "${compress[@]}"
        run_phase "$dir" verify "$BIN" "${verify[@]}"
        run_phase "$dir" decode "$BIN" -decode "${read_args[@]}"
        PHASE="$dataset.$mode.report"
        python3 "$SCRIPT_DIR/report_opcode_profile.py" "$OUT/$dataset" "$MAX_VECTORS" "$mode" > "$dir/report.txt"
        shasum -a 256 "$dir/candidate.bv" > "$dir/stream.sha256"
        if (( ! KEEP )); then
            rm -- "$dir/candidate.bv"
        fi
    done
    cp "$OUT/$dataset/profile/report.txt" "$OUT/$dataset/report.txt"
done
PHASE=summary
: > "$OUT/report.txt"
for dataset in "${sets[@]}"; do
    printf '\n%s\n' "Gov2.$dataset" | tee -a "$OUT/report.txt"
    cat "$OUT/$dataset/report.txt" | tee -a "$OUT/report.txt"
done
printf '\nReports and logs: %s\n' "$OUT"
