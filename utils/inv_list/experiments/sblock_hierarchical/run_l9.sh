#!/usr/bin/env bash
# Full level-9 wave policy comparison; build bminv separately.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUT=${1:?Usage: run_l9.sh FRESH_USB_DIR [BINARY] [THREADS]}
BIN=${2:-./bminv}
THREADS=${3:-10}
[[ $OUT == /Volumes/USB_1TB/* && ! -e $OUT && -x $BIN ]] || exit 2
mkdir -p "$OUT"
for dataset in sorted unsorted; do
    for mode in preset exact tree; do
        args=(--out "$OUT/$dataset-$mode" --bin "$BIN" --level 9
              --threads "$THREADS" --progress 1000 --datasets "$dataset"
              --keep-streams)
        if [[ $mode == exact ]]; then
            args+=(--wave-mode both --wave-leaf-out "$OUT/$dataset-leaves.tsv")
        elif [[ $mode == tree ]]; then
            args+=(--wave-mode tree)
        fi
        bash "$HERE/../opcode_profile/run_gov2_opcode_profile.sh" "${args[@]}"
    done
done
python3 "$HERE/report_l9.py" "$OUT" > "$OUT/comparison.txt"
cat "$OUT/comparison.txt"
