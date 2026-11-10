#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$SCRIPT_DIR"

BMINV=${BMINV:-./bminv}
SOURCE=${SOURCE:-gov2.unsorted.bv.cm6}
PAIR_INDEX=${PAIR_INDEX:-gov2.unsorted.pair-shards.csv}
RS_CLUSTER_FEATURE=${RS_CLUSTER_FEATURE:-regions3}
RS_CLUSTER_VECTORS=${RS_CLUSTER_VECTORS:-7000}
RS_CLUSTER_MEMORY_BYTES=${RS_CLUSTER_MEMORY_BYTES:-8589934592}
RS_MIN_BYTES=${RS_MIN_BYTES:-4096}
RS_BYTE_COVERAGE=${RS_BYTE_COVERAGE:-90}
RS_SEED=${RS_SEED:-42}
PAIR_THREADS=${PAIR_THREADS:-8}
RS_CLUSTER_PREFIX=${RS_CLUSTER_PREFIX:-gov2.rs.cluster.$RS_CLUSTER_FEATURE}
MANIFEST=$RS_CLUSTER_PREFIX.manifest.csv
PAIRS=$RS_CLUSTER_PREFIX.pairs.csv

stage()
{
  printf '\n========================================================================\n'
  printf 'Stage %s/2: %s\n' "$1" "$2"
  printf '========================================================================\n\n'
}

command -v "$BMINV" >/dev/null 2>&1 || { echo "ERROR: bminv executable missing: $BMINV" >&2; exit 127; }
[ -r "$SOURCE" ] && [ -r "$PAIR_INDEX" ] || {
  echo "ERROR: collection or pair index missing: $SOURCE / $PAIR_INDEX" >&2
  exit 1
}

if [ "${RS_CLEAN:-0}" = 1 ]; then
  rm -f "$MANIFEST" "$MANIFEST.report.txt" "$PAIRS" "$PAIRS.best.csv" "$PAIRS.report.txt"
fi

stage 1 "Build compact RS profiles and plan bounded clusters ($RS_CLUSTER_FEATURE)"
if "$BMINV" -rs-cluster-plan -bvin "$SOURCE" -pair-index-in "$PAIR_INDEX" \
  -rs-cluster-out "$MANIFEST" -rs-cluster-feature "$RS_CLUSTER_FEATURE" \
  -rs-cluster-seed "$RS_SEED" -rs-cluster-vectors "$RS_CLUSTER_VECTORS" \
  -rs-cluster-memory "$RS_CLUSTER_MEMORY_BYTES" -rs-min-bytes "$RS_MIN_BYTES" \
  -rs-byte-coverage "$RS_BYTE_COVERAGE" -level 6
then
  echo "Cluster manifest complete: $MANIFEST"
else
  stage_rc=$?
  echo "ERROR: RS cluster planning failed (exit $stage_rc)." >&2
  exit "$stage_rc"
fi

stage 2 "Search XOR pairs within RS clusters using $PAIR_THREADS threads"
if "$BMINV" -pair-search -bvin "$SOURCE" -pair-index-in "$PAIR_INDEX" \
  -pair-cluster-in "$MANIFEST" -csv "$PAIRS" -pair-max-ratio 0 \
  -pair-min-pct 5 -pair-min-gain 256 -pair-threads "$PAIR_THREADS" -level 6
then
  echo "RS clustered pair search complete: $PAIRS"
else
  stage_rc=$?
  echo "ERROR: RS clustered pair search failed (exit $stage_rc)." >&2
  exit "$stage_rc"
fi
