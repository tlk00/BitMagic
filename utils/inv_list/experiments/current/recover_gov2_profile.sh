#!/usr/bin/env bash
# Finish verification/reporting for an existing profiled compression, without recompression.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BIN=./bminv
RUN=
DATASET=sorted
INPUT_DIR=/Volumes/USB_1TB/int_set_bench/Gov2Flat
THREADS=8
PROGRESS=1000
while (( $# )); do
    case "$1" in
        --run|--bin|--dataset|--input-dir|--threads|--progress)
            (( $# >= 2 )) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --run) RUN=$2;; --bin) BIN=$2;; --dataset) DATASET=$2;;
                --input-dir) INPUT_DIR=$2;; --threads) THREADS=$2;; --progress) PROGRESS=$2;;
            esac
            shift 2;;
        --help|-h)
            echo 'Usage: recover_gov2_profile.sh --run EXISTING_RUN --dataset sorted|unsorted [--bin FILE] [--input-dir DIR] [--threads 8] [--progress 1000]'
            exit 0;;
        *) echo "Unknown argument: $1" >&2; exit 2;;
    esac
done
[[ $DATASET == sorted || $DATASET == unsorted ]] || exit 2
[[ $THREADS =~ ^[1-9][0-9]*$ && $PROGRESS =~ ^[1-9][0-9]*$ ]] || exit 2
DIR="$RUN/$DATASET/profile"
[[ -n $RUN && -x $BIN && -s $RUN/run.txt && -s $DIR/compress.log &&
   -s $DIR/compress.time && -s $DIR/candidate.bv && -s $DIR/candidate.idx &&
   -s $INPUT_DIR/gov2.$DATASET ]] || { echo 'Missing binary, input, or saved profile artifacts' >&2; exit 2; }
# Confirm artifacts belong to the recorded completed compression before rerunning reads.
python3 - "$DIR" <<'CHECK'
from pathlib import Path
import re, sys
p = Path(sys.argv[1])
log = (p / 'compress.log').read_text(errors='replace')
size = re.search(r'^BV size = (\d+)$', log, re.M)
count = re.search(r'^Total vectors=(\d+)$', log, re.M)
if not size or not count or int(count[1]) == 0:
    raise SystemExit('Compression log has no completed stream summary')
if (p / 'candidate.bv').stat().st_size != int(size[1]):
    raise SystemExit('Saved stream size differs from compression log')
with (p / 'candidate.idx').open() as stream:
    header = stream.readline().split()
if header != ['BM_INV_BLOB_INDEX_V1', size[1], count[1]]:
    raise SystemExit('Saved index header differs from compression log')
CHECK
LIMIT=$(awk '/^Vectors: / {print $2; exit}' "$RUN/run.txt")
[[ $LIMIT =~ ^[0-9]+$ ]] || { echo 'Missing vector limit in run.txt' >&2; exit 2; }
STAMP=$(date +%Y%m%d-%H%M%S)
run_phase()
{
    local phase=$1
    shift
    # Preserve the original failure evidence.
    for extension in log time; do
        if [[ -e $DIR/$phase.$extension ]]; then
            cp -p -- "$DIR/$phase.$extension" "$DIR/$phase.$extension.before-recovery-$STAMP"
        fi
    done
    { printf 'recovery.%s.%s:' "$DATASET" "$phase"; printf ' %q' "$@"; printf '\n'; } >> "$RUN/run.txt"
    TIMEFORMAT='%R %U %S'
    { time "$@" 2>&1 | tee "$DIR/$phase.log" |
        python3 "$HERE/progress_console.py" --hide-opcode-profile; } 2> "$DIR/$phase.time"
}
read_args=(-t -s -bvin "$DIR/candidate.bv" -bv-index-in "$DIR/candidate.idx"
    -threads "$THREADS" -batch-vectors 512 -progress-vectors "$PROGRESS")
verify=(-verify -u32in "$INPUT_DIR/gov2.$DATASET" "${read_args[@]}")
if (( LIMIT )); then verify+=(-max-vectors "$LIMIT"); fi
run_phase verify "$BIN" "${verify[@]}"
run_phase decode "$BIN" -decode "${read_args[@]}"
python3 "$HERE/../opcode_profile/report_opcode_profile.py" "$RUN/$DATASET" "$LIMIT" profile > "$DIR/report.txt"
shasum -a 256 "$DIR/candidate.bv" > "$DIR/stream.sha256"
cp -- "$DIR/report.txt" "$RUN/$DATASET/report.txt"
python3 "$HERE/report_gov2_highest.py" "$RUN" | tee "$RUN/verification-report.txt"
printf 'Recovered %s: %s\nSaved stream retained. Original FAILURE.txt, if any, retained as history.\n' "$DATASET" "$RUN"
