#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUT=
INPUT=/Volumes/USB_1TB/int_set_bench/Gov2Flat
COUNT=10000
opts=()
while (( $# )); do
    case "$1" in
        --out|--input-dir|--sample-vectors)
            (( $# >= 2 )) || exit 2
            case "$1" in --out) OUT=$2;; --input-dir) INPUT=$2;; --sample-vectors) COUNT=$2;; esac
            shift 2;;
        --bin|--threads|--progress)
            (( $# >= 2 )) || exit 2; opts+=("$1" "$2"); shift 2;;
        --help|-h) echo 'Usage: run_gov2_correction_sample.sh --out FRESH_USB_DIR [--input-dir DIR] [--sample-vectors 10000] [--bin FILE] [--threads 10] [--progress 1000]'; exit 0;;
        *) echo "Unknown argument: $1" >&2; exit 2;;
    esac
done
[[ -n $OUT && ! -e $OUT ]] || { echo 'Require fresh --out directory' >&2; exit 2; }
[[ $COUNT =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid sample count' >&2; exit 2; }
python3 - "$OUT" <<'PY'
import os,subprocess,sys
if not os.path.realpath(sys.argv[1]).startswith('/Volumes/USB_1TB/'):
    raise SystemExit('Output must be on USB_1TB')
if ' on /Volumes/USB_1TB (' not in subprocess.check_output(['mount'],text=True):
    raise SystemExit('USB_1TB is not mounted')
PY
mkdir -p -- "$OUT"
trap 'rc=$?; if (( rc )); then printf "Sampled correction experiment failed: exit=%s\n" "$rc" > "$OUT/FAILURE.txt"; fi; exit "$rc"' EXIT
python3 "$SCRIPT_DIR/sample_gov2.py" --input-dir "$INPUT" --out "$OUT/sample" --count "$COUNT" | tee "$OUT/sampling.log"
bash "$SCRIPT_DIR/run_gov2_correction.sh" --input-dir "$OUT/sample" \
    --max-vectors "$COUNT" --out "$OUT/comparison" ${opts[@]+"${opts[@]}"}
cp "$OUT/comparison/report.txt" "$OUT/report.txt"
printf 'Reports and logs: %s\n' "$OUT"
