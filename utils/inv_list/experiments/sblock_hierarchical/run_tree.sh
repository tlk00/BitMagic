#!/usr/bin/env bash
# Full L8 experimental tree versus exact 108/109 selection.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUT=${1:?Usage: run_tree.sh FRESH_USB_DIR [BINARY] [THREADS]}
BIN=${2:-./bminv}
THREADS=${3:-10}
[[ $OUT == /Volumes/USB_1TB/* && ! -e $OUT && -x $BIN ]] || exit 2
mkdir -p "$OUT"
for dataset in sorted unsorted; do
    for mode in exact tree; do
        selector=$mode
        if [[ $mode == exact ]]; then selector=both; fi
        args=(--out "$OUT/$dataset-$mode" --bin "$BIN" --level 8
              --threads "$THREADS" --progress 1000 --datasets "$dataset"
              --keep-streams --wave-mode "$selector")
        if [[ $mode == exact ]]; then
            args+=(--wave-leaf-out "$OUT/$dataset-leaves.tsv")
        fi
        bash "$HERE/../opcode_profile/run_gov2_opcode_profile.sh" "${args[@]}"
    done
done
python3 "$HERE/report_tree.py" "$OUT" > "$OUT/comparison.txt"
cat "$OUT/comparison.txt"
