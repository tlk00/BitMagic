#!/usr/bin/env bash
# Build separately with make neon in utils/inv_list.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
BIN=./bminv
OUT=
THREADS=10
INPUT_DIR=/Volumes/USB_1TB/int_set_bench/Gov2Flat
while (( $# )); do
    case "$1" in
        --bin|--out|--threads|--input-dir)
            (( $# >= 2 )) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --bin) BIN=$2;;
                --out) OUT=$2;;
                --threads) THREADS=$2;;
                --input-dir) INPUT_DIR=$2;;
            esac
            shift 2;;
        --help|-h)
            echo 'Usage: run_compare.sh --out FRESH_USB_DIR [--bin ./bminv] [--threads 10] [--input-dir DIR]'
            exit 0;;
        *) echo "Unknown argument: $1" >&2; exit 2;;
    esac
done
[[ -n $OUT && ! -e $OUT && -x $BIN ]] || {
    echo 'Require a fresh --out directory and executable --bin' >&2
    exit 2
}
[[ $THREADS =~ ^[1-9][0-9]*$ && ${#THREADS} -le 3 ]] && (( THREADS <= 256 )) || exit 2
[[ -d /Volumes/USB_1TB && $OUT == /Volumes/USB_1TB/* ]] || exit 2
for dataset in sorted unsorted; do
    [[ -s $INPUT_DIR/gov2.$dataset ]] || { echo "Missing gov2.$dataset" >&2; exit 2; }
done
mount | rg -F -q " on /Volumes/USB_1TB (" || { echo "USB_1TB is not mounted" >&2; exit 2; }
mkdir -p -- "$OUT"
for mode in 108 both 109; do
    bash "$HERE/../current/run_gov2_highest.sh" --bin "$BIN" \
        --out "$OUT/$mode" --datasets both --level 10 \
        --threads "$THREADS" --input-dir "$INPUT_DIR" --wave-mode "$mode"
done
python3 "$HERE/report_compare.py" "$OUT/108" "$OUT/both" "$OUT/109" > "$OUT/comparison.txt"
cat "$OUT/comparison.txt"
