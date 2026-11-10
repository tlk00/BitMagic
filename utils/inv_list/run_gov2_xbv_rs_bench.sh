#!/bin/sh
# Common driver: sh run_gov2_xbv_rs_bench.sh shard|stratified
set -eu
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$SCRIPT_DIR"
BENCH_SET=${1:?Specify shard or stratified}
case "$BENCH_SET" in shard|stratified) ;; *) echo "Unknown benchmark set: $BENCH_SET" >&2; exit 2;; esac
BMINV=${BMINV:-./bminv}
SOURCE=${SOURCE:-gov2.unsorted.bv.cm6}
PAIR_INDEX=${PAIR_INDEX:-gov2.unsorted.pair-shards.csv}
RS_PREFIX=${RS_PREFIX:-gov2.rs}
RS_SIZE=${RS_SIZE:-7000}
RS_THREADS=${RS_THREADS:-8}
RS_SEED=${RS_SEED:-42}
RS_THRESHOLDS=${RS_THRESHOLDS:-0.1,0.25,0.5,1}
RS_REPEATS=${RS_REPEATS:-2}
RS_MIN_BYTES=${RS_MIN_BYTES:-4096}
RS_BYTE_COVERAGE=${RS_BYTE_COVERAGE:-90}
PYTHON=${PYTHON:-python3}
command -v "$BMINV" >/dev/null 2>&1 || { echo "Missing bminv: $BMINV" >&2; exit 127; }
command -v "$PYTHON" >/dev/null 2>&1 || { echo "Missing Python: $PYTHON" >&2; exit 127; }
[ -r "$SOURCE" ] && [ -r "$PAIR_INDEX" ] || {
  echo "Source or shard index missing. Create an index with -pair-index first." >&2; exit 1;
}
printf '\n========================================================================\n'
printf 'RS experiment: %s set; %s vectors; %s threads\n' "$BENCH_SET" "$RS_SIZE" "$RS_THREADS"
printf '========================================================================\n'
set --
if [ "$BENCH_SET" = shard ] && [ -n "${RS_SHARD:-}" ]; then
  set -- -rs-shard "$RS_SHARD"
fi
if "$BMINV" -rs-bench -bvin "$SOURCE" -pair-index-in "$PAIR_INDEX" \
  -rs-out "$RS_PREFIX.$BENCH_SET" -rs-set "$BENCH_SET" -rs-size "$RS_SIZE" \
  -rs-seed "$RS_SEED" -rs-thresholds "$RS_THRESHOLDS" -rs-repeats "$RS_REPEATS" \
  -rs-min-bytes "$RS_MIN_BYTES" -rs-byte-coverage "$RS_BYTE_COVERAGE" \
  -pair-threads "$RS_THREADS" -pair-min-pct 5 -pair-min-gain 256 -level 6 "$@"
then
  printf '\nRS benchmark complete. Analyzing detailed pair outcomes...\n'
else
  bench_rc=$?
  echo "ERROR: RS $BENCH_SET benchmark failed (exit $bench_rc)." >&2
  exit "$bench_rc"
fi
if "$PYTHON" "$SCRIPT_DIR/analyze_xbv_rs_bench.py" "$RS_PREFIX.$BENCH_SET"; then
  echo "RS $BENCH_SET analysis complete."
else
  analysis_rc=$?
  echo "ERROR: RS $BENCH_SET analysis failed (exit $analysis_rc)." >&2
  exit "$analysis_rc"
fi
