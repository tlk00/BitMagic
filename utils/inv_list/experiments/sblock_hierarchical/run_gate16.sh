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
            (( $# >= 2 )) || exit 2
            case "$1" in
                --bin) BIN=$2;;
                --out) OUT=$2;;
                --threads) THREADS=$2;;
                --input-dir) INPUT_DIR=$2;;
            esac
            shift 2;;
        --help|-h)
            echo 'Usage: run_gate16.sh --out FRESH_USB_DIR [--bin ./bminv] [--threads 10] [--input-dir DIR]'
            exit 0;;
        *) echo "Unknown argument: $1" >&2; exit 2;;
    esac
done
[[ -x $BIN && -n $OUT && ! -e $OUT && $OUT == /Volumes/USB_1TB/* ]] || exit 2
[[ $THREADS =~ ^[1-9][0-9]*$ && ${#THREADS} -le 3 ]] && (( THREADS <= 256 )) || exit 2
mount | rg -F -q ' on /Volumes/USB_1TB (' || { echo 'USB_1TB is not mounted' >&2; exit 2; }
mkdir -p -- "$OUT"
for mode in exact gate16; do
    selector=both
    if [[ $mode == gate16 ]]; then selector=gate16; fi
    bash "$HERE/../current/run_gov2_highest.sh" --bin "$BIN" --datasets both \
        --level 8 --threads "$THREADS" --wave-mode "$selector" \
        --input-dir "$INPUT_DIR" --keep-streams --out "$OUT/$mode"
done
python3 "$HERE/report_gate16.py" "$OUT" > "$OUT/comparison.txt"
cat "$OUT/comparison.txt"
