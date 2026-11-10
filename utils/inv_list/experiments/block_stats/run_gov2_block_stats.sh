#!/usr/bin/env bash
# Report optimized bvector block composition of the raw Gov2 collection.
set -euo pipefail
export LC_ALL=C

BIN=./gov2_block_stats
INPUT=/Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.sorted
OUT=/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/gov2-sorted-block-stats-$(date +%Y%m%d-%H%M%S)
MAX_VECTORS=0
PROGRESS=1000
while (( $# )); do
    case "$1" in
        --bin|--input|--out|--max-vectors|--progress)
            (( $# >= 2 )) || { echo "Missing value for $1" >&2; exit 2; }
            case "$1" in
                --bin) BIN=$2;; --input) INPUT=$2;; --out) OUT=$2;;
                --max-vectors) MAX_VECTORS=$2;; --progress) PROGRESS=$2;;
            esac
            shift 2;;
        *) echo "Unknown option: $1" >&2; exit 2;;
    esac
done
for n in "$MAX_VECTORS" "$PROGRESS"; do
    [[ $n =~ ^(0|[1-9][0-9]*)$ && ${#n} -le 12 ]] || {
        echo "Invalid count: $n" >&2; exit 2;
    }
done
[[ -x $BIN && -s $INPUT ]] || { echo 'Missing binary or input; run make gov2-block-stats-neon' >&2; exit 2; }
[[ $OUT == /Volumes/USB_1TB/* && ! -e $OUT ]] || {
    echo 'Output must be a new directory under /Volumes/USB_1TB' >&2; exit 2;
}
mount | rg -F -q ' on /Volumes/USB_1TB (' || {
    echo 'USB_1TB is not mounted' >&2; exit 2;
}
mkdir -p -- "$OUT"
printf 'Input: %s\nMax vectors: %s\nProgress: %s\n' "$INPUT" "$MAX_VECTORS" "$PROGRESS" > "$OUT/run.txt"
{ printf 'Command:'; printf ' %q' "$BIN" --input "$INPUT" --out "$OUT" --max-vectors "$MAX_VECTORS" --progress "$PROGRESS"; printf '\n'; } >> "$OUT/run.txt"
if "$BIN" --input "$INPUT" --out "$OUT" --max-vectors "$MAX_VECTORS" \
        --progress "$PROGRESS" > "$OUT/scan.log" 2>&1; then
    cat "$OUT/report.txt"
    printf 'Reports and logs: %s\n' "$OUT"
else
    tail -n 20 "$OUT/scan.log" >&2
    echo "FAILED; see $OUT" >&2
    exit 1
fi
