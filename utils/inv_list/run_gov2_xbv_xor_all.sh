#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$SCRIPT_DIR"

BMINV=${BMINV:-./bminv}
SOURCE=${SOURCE:-gov2.unsorted.bv.cm6}
PAIRS=${PAIRS:-gov2.unsorted.pairs.5pct.2x.csv}
PAIR_INDEX=${PAIR_INDEX:-gov2.unsorted.pair-shards.csv}
XOR_PAYLOAD=${XOR_PAYLOAD:-gov2.unsorted.xor.bv}
XOR_REFS=${XOR_REFS:-gov2.unsorted.xor.refs}
XOR_INDEX=${XOR_INDEX:-gov2.unsorted.xor.offsets}

export BMINV SOURCE PAIRS PAIR_INDEX XOR_PAYLOAD XOR_REFS XOR_INDEX

if ! BMINV_COMMAND=$(command -v "$BMINV" 2>/dev/null); then
  echo "ERROR: bminv executable not found: $BMINV" >&2
  exit 127
fi
if [ ! -r "$SOURCE" ]; then
  echo "ERROR: source collection is not readable: $SOURCE" >&2
  exit 1
fi

canonical_path()
(
  canonical_input=$1
  canonical_dir=$(dirname -- "$canonical_input")
  canonical_base=$(basename -- "$canonical_input")
  cd "$canonical_dir" || exit 1
  printf '%s/%s\n' "$(pwd -P)" "$canonical_base"
)

SOURCE_PROTECTED=$(canonical_path "$SOURCE") || {
  echo "ERROR: cannot resolve source collection path: $SOURCE" >&2
  exit 1
}
BMINV_PROTECTED=$(canonical_path "$BMINV_COMMAND") || {
  echo "ERROR: cannot resolve bminv path: $BMINV_COMMAND" >&2
  exit 1
}

cleanup_file()
{
  cleanup_path=$1
  cleanup_resolved=$(canonical_path "$cleanup_path") || {
    echo "ERROR: cannot resolve output path: $cleanup_path" >&2
    exit 1
  }
  if [ "$cleanup_resolved" = "$SOURCE_PROTECTED" ] ||
     [ "$cleanup_resolved" = "$BMINV_PROTECTED" ]; then
    echo "ERROR: refusing to remove protected input: $cleanup_path" >&2
    exit 2
  fi
  if [ -e "$cleanup_path" ]; then
    echo "Removing previous output: $cleanup_path"
    if ! rm -f -- "$cleanup_path"; then
      echo "ERROR: cannot remove previous output: $cleanup_path" >&2
      exit 1
    fi
  fi
}

printf '\n========================================================================\n'
printf 'Preparation: Remove outputs from previous Gov2 XOR runs\n'
printf '========================================================================\n\n'

cleanup_file "$PAIR_INDEX"
cleanup_file "$PAIR_INDEX.report.txt"
cleanup_file "$PAIRS"
cleanup_file "$PAIRS.best.csv"
cleanup_file "$PAIRS.measured.csv"
cleanup_file "$PAIRS.report.txt"
cleanup_file "$XOR_PAYLOAD"
cleanup_file "$XOR_PAYLOAD.report.txt"
cleanup_file "$XOR_REFS"
cleanup_file "$XOR_INDEX"

run_phase()
{
  phase_name=$1
  phase_script=$2

  if [ ! -x "$phase_script" ]; then
    echo "ERROR: $phase_name script is not executable: $phase_script" >&2
    exit 126
  fi
  if "$phase_script"; then
    return 0
  else
    status=$?
    echo "ERROR: stopping after $phase_name failure (exit code $status)." >&2
    exit "$status"
  fi
}

run_phase "XOR search" "$SCRIPT_DIR/run_gov2_xbv_xor_search.sh"
run_phase "XOR build" "$SCRIPT_DIR/run_gov2_xbv_xor_build.sh"
run_phase "XOR verification" "$SCRIPT_DIR/run_gov2_xbv_xor_verify.sh"

printf '\n========================================================================\n'
echo "All Gov2 XOR phases completed successfully."
printf '========================================================================\n\n'
