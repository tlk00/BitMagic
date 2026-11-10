#!/usr/bin/env bash
# Compress, verify and decode one raw collection with explicit presets.
# Build separately. Input is never modified; output must be a fresh directory.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BIN="$SCRIPT_DIR/bminv"
INPUT=
OUT=
LEVELS=6
THREADS=
LIMIT=
PROGRESS=10000
usage() {
    cat <<'HELP'
Usage: run_collection.sh --input FILE --out NEW_DIR [options]
  --bin FILE       built bminv variant (default: bminv beside this script)
  --levels LIST    comma-separated unique presets 0..10 (default: 6)
  --threads N      use 1..256 workers and an index (default: serial)
  --max-vectors N  bound compression and raw verification to a positive prefix
  --progress N     progress interval in vectors (default: 10000)
Compresses, verifies against raw input, then independently decodes each stream.
Retains streams, indices, commands and logs under NEW_DIR/level-N.
Accepts Gov2, ClueWeb or any collection in the documented raw-list format.
HELP
}
while (( $# )); do
    case "$1" in
        --help|-h) usage; exit 0;;
        --input|--out|--bin|--levels|--threads|--max-vectors|--progress)
            (( $# >= 2 )) || { echo "Missing value for $1" >&2; exit 2; }
            [[ -n $2 ]] || { echo "Empty value for $1" >&2; exit 2; }
            case "$1" in
                --input) INPUT=$2;; --out) OUT=$2;; --bin) BIN=$2;;
                --levels) LEVELS=$2;; --threads) THREADS=$2;;
                --max-vectors) LIMIT=$2;; --progress) PROGRESS=$2;;
            esac
            shift 2;;
        *) echo "Unknown argument: $1" >&2; exit 2;;
    esac
done
[[ -f $INPUT && -s $INPUT && -n $OUT && -x $BIN && ! -d $BIN ]] || {
    echo 'Require nonempty --input, --out and executable --bin' >&2; exit 2;
}
[[ $LEVELS =~ ^(10|[0-9])(,(10|[0-9]))*$ ]] || {
    echo 'Levels must be comma-separated presets 0..10' >&2; exit 2;
}
IFS=, read -r -a levels <<< "$LEVELS"
seen=,
for level in "${levels[@]}"; do
    [[ $seen != *",$level,"* ]] || { echo "Duplicate level: $level" >&2; exit 2; }
    seen+="$level,"
done
for value in "$PROGRESS" ${LIMIT:+"$LIMIT"} ${THREADS:+"$THREADS"}; do
    [[ $value =~ ^[1-9][0-9]*$ && ${#value} -le 9 ]] || {
        echo 'Thread, progress and vector limits must be positive integers (at most 9 digits)' >&2; exit 2;
    }
done
[[ -z $THREADS ]] || (( THREADS <= 256 )) || { echo 'Threads must be 1..256' >&2; exit 2; }
# mkdir refuses existing directories, preventing accidental result replacement.
mkdir -- "$OUT"
trap 'rc=$?; if (( rc )); then printf "Collection run failed: exit=%s\n" "$rc" > "$OUT/FAILURE.txt"; fi' EXIT
prefix=()
if [[ -n $LIMIT ]]; then prefix=(-max-vectors "$LIMIT"); fi
for level in "${levels[@]}"; do
    destination="$OUT/level-$level"
    mkdir -- "$destination"
    stream="$destination/collection.bv"
    encode=(-s -t -level "$level" -u32in "$INPUT" -bvout "$stream" -progress-vectors "$PROGRESS")
    verify=(-s -t -u32in "$INPUT" -bvin "$stream" -verify -progress-vectors "$PROGRESS")
    decode=(-s -t -bvin "$stream" -decode -progress-vectors "$PROGRESS")
    if [[ -n $THREADS ]]; then
        encode+=(-threads "$THREADS" -bv-index-out "$stream.idx")
        verify+=(-threads "$THREADS" -bv-index-in "$stream.idx")
        decode+=(-threads "$THREADS" -bv-index-in "$stream.idx")
    fi
    encode+=(${prefix[@]+"${prefix[@]}"})
    verify+=(${prefix[@]+"${prefix[@]}"})
    for phase in compress verify decode; do
        case "$phase" in
            compress) arguments=("${encode[@]}");;
            verify) arguments=("${verify[@]}");;
            decode) arguments=("${decode[@]}");;
        esac
        printf '%q ' "$BIN" "${arguments[@]}" > "$destination/$phase.command"
        printf '\n' >> "$destination/$phase.command"
        "$BIN" "${arguments[@]}" 2>&1 | tee "$destination/$phase.log"
    done
done
printf 'Compression, raw verification and decoding completed: %s\n' "$OUT"
