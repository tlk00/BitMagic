#!/usr/bin/env bash
# Current highest-compression configuration; build separately with make neon.
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUT=
REFERENCE=
LEVEL=10
opts=()
while (( $# )); do
    case "$1" in
        --level) (( $# >= 2 )) || exit 2; LEVEL=$2; shift 2;;
        --out|--reference)
            (( $# >= 2 )) || exit 2
            case "$1" in --out) OUT=$2;; --reference) REFERENCE=$2;; esac
            shift 2;;
        --bin|--input-dir|--datasets|--threads|--max-vectors|--progress|--sblock-large-inverse-limit|--bwt-min-gain)
            (( $# >= 2 )) || exit 2;opts+=("$1" "$2");shift 2;;
        --diagnostics|--keep-streams) opts+=("$1");shift;;
        --help|-h) echo 'Usage: run_gov2_highest.sh --out FRESH_USB_DIR [--level N (default 10)] [--bwt-min-gain N] [--reference VERIFIED_RUN_ROOT] [--datasets both|sorted|unsorted] [--threads N] [--max-vectors N] [--progress N] [--bin FILE] [--input-dir DIR] [--diagnostics] [--keep-streams]';exit 0;;
        *) echo "Unknown argument: $1" >&2;exit 2;;
    esac
done
[[ -n $OUT ]] || { echo 'Require --out' >&2;exit 2; }
bash "$SCRIPT_DIR/../opcode_profile/run_gov2_opcode_profile.sh" \
    --datasets both --threads 10 --progress 1000 --level "$LEVEL" --out "$OUT" ${opts[@]+"${opts[@]}"}
extra=();if [[ -n $REFERENCE ]]; then extra+=(--reference "$REFERENCE"); fi
if ! python3 "$SCRIPT_DIR/report_gov2_highest.py" "$OUT" ${extra[@]+"${extra[@]}"} | tee "$OUT/verification-report.txt"; then
    echo 'FAILED: highest-compression verification report' > "$OUT/FAILURE.txt"
    exit 1
fi
printf 'Reports and logs: %s\n' "$OUT"
