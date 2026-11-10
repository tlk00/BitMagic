#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$SCRIPT_DIR"

BMINV=${BMINV:-./bminv}
SOURCE=${SOURCE:-gov2.unsorted.bv.cm6}
PAIRS=${PAIRS:-gov2.unsorted.pairs.5pct.2x.csv}
PAIR_INDEX=${PAIR_INDEX:-gov2.unsorted.pair-shards.csv}
XOR_MEMORY_BYTES=${XOR_MEMORY_BYTES:-8589934592}
XOR_PAYLOAD=${XOR_PAYLOAD:-gov2.unsorted.xor.bv}
XOR_REFS=${XOR_REFS:-gov2.unsorted.xor.refs}
XOR_INDEX=${XOR_INDEX:-gov2.unsorted.xor.offsets}

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

print_stage 3 "Build XOR payload, reference metadata, and offset index"

if "$BMINV" \
  -xor-build \
  -bvin "$SOURCE" \
  -pair-index-in "$PAIR_INDEX" \
  -pairs-in "$PAIRS" \
  -xorout "$XOR_PAYLOAD" \
  -xor-refs "$XOR_REFS" \
  -xor-index "$XOR_INDEX" \
  -ref-top 4 \
  -pair-threads 8 \
  -xor-memory "$XOR_MEMORY_BYTES" \
  -level 6
then
  echo "XOR build phase completed successfully."
else
  status=$?
  echo "ERROR: XOR build phase failed with exit code $status." >&2
  exit "$status"
fi
