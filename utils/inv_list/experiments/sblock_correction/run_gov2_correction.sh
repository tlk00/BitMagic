#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUT=
LIMIT=10000
opts=()
while (( $# )); do
    case "$1" in
        --out|--max-vectors)
            (( $# >= 2 )) || exit 2
            case "$1" in --out) OUT=$2;; --max-vectors) LIMIT=$2;; esac
            shift 2;;
        --bin|--input-dir|--threads|--progress)
            (( $# >= 2 )) || exit 2; opts+=("$1" "$2"); shift 2;;
        --help|-h) echo 'Usage: run_gov2_correction.sh --out FRESH_USB_DIR [--bin FILE] [--input-dir DIR] [--threads N] [--progress N] [--max-vectors N (default 10000)]'; exit 0;;
        *) echo "Unknown argument: $1" >&2; exit 2;;
    esac
done
[[ -n $OUT && ! -e $OUT ]] || { echo 'Require fresh --out directory' >&2; exit 2; }
trap 'rc=$?; if (( rc )) && [[ -d $OUT ]]; then printf "Correction comparison failed: exit=%s\n" "$rc" > "$OUT/FAILURE.txt"; fi; exit "$rc"' EXIT
# Child harness validates USB mount/output and builds no binaries.
for mode in baseline selected; do
    extra=(); if [[ $mode == selected ]]; then extra+=(--sblock-correction-select); fi
    bash "$SCRIPT_DIR/../opcode_profile/run_gov2_opcode_profile.sh" --datasets both \
        --threads 10 --progress 1000 --max-vectors "$LIMIT" \
        --out "$OUT/$mode" ${opts[@]+"${opts[@]}"} ${extra[@]+"${extra[@]}"}
    python3 "$SCRIPT_DIR/../current/report_gov2_highest.py" "$OUT/$mode" \
        | tee "$OUT/$mode/verification-report.txt"
done
python3 "$SCRIPT_DIR/report_correction.py" "$OUT" | tee "$OUT/report.txt"
printf 'Reports and logs: %s\n' "$OUT"
