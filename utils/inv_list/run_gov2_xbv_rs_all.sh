#!/bin/sh
set -eu
SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
for bench_set in shard stratified; do
  printf '\n========================================================================\n'
  printf 'Benchmark set: %s\n' "$bench_set"
  printf '========================================================================\n'
  if sh "$SCRIPT_DIR/run_gov2_xbv_rs_bench.sh" "$bench_set"; then
    echo "Completed $bench_set."
  else
    experiment_rc=$?
    echo "ERROR: stopping RS experiments after $bench_set failure (exit $experiment_rc)." >&2
    exit "$experiment_rc"
  fi
done
echo "Both RS experiments and analyses completed."
