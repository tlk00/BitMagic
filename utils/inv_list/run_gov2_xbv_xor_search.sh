#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$SCRIPT_DIR"

BMINV=${BMINV:-./bminv}
SOURCE=${SOURCE:-gov2.unsorted.bv.cm6}
PAIRS=${PAIRS:-gov2.unsorted.pairs.5pct.2x.csv}
PAIR_INDEX=${PAIR_INDEX:-gov2.unsorted.pair-shards.csv}
PAIR_SHARD_VECTORS=${PAIR_SHARD_VECTORS:-7000}
PAIR_SHARD_BYTES=${PAIR_SHARD_BYTES:-3221225472}
PAIR_SHARD_TAIL_PCT=${PAIR_SHARD_TAIL_PCT:-5}

print_stage()
{
  printf '\n========================================================================\n'
  printf 'Stage %s/4: %s\n' "$1" "$2"
  printf '========================================================================\n\n'
}

if ! command -v "$BMINV" >/dev/null 2>&1; then
  echo "ERROR: bminv executable not found: $BMINV" >&2
  exit 127
fi

print_stage 1 "Index collection by serialized size and plan greedy shards"

if "$BMINV" \
  -pair-index \
  -bvin "$SOURCE" \
  -pair-index-out "$PAIR_INDEX" \
  -pair-max-ratio 2 \
  -pair-shard-vectors "$PAIR_SHARD_VECTORS" \
  -pair-shard-bytes "$PAIR_SHARD_BYTES" \
  -pair-shard-tail-pct "$PAIR_SHARD_TAIL_PCT" \
  -level 6
then
  echo "Pair indexing and shard planning completed successfully."
else
  status=$?
  echo "ERROR: pair indexing and shard planning failed with exit code $status." >&2
  exit "$status"
fi

print_stage 2 "Search XOR candidates within each size shard"

if "$BMINV" \
  -pair-search \
  -bvin "$SOURCE" \
  -pair-index-in "$PAIR_INDEX" \
  -csv "$PAIRS" \
  -pair-min-pct 5 \
  -pair-max-ratio 2 \
  -pair-threads 8 \
  -level 6
then
  echo "XOR search phase completed successfully."
else
  status=$?
  echo "ERROR: XOR search phase failed with exit code $status." >&2
  exit "$status"
fi
