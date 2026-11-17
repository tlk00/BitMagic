#!/bin/sh
# Benchmark four serialized NED dictionary variants.
set -eu

cd "$(dirname "$0")"
case "$(uname -s):$(uname -m)" in
    Darwin:arm64|Darwin:aarch64) default_prog=./xsample05_neon ;;
    Linux:x86_64|Linux:amd64) default_prog=./xsample05_release ;;
    *) default_prog=./xsample05_release ;;
esac
PROG=${PROG:-$default_prog}

if [ ! -x "$PROG" ]; then
    echo "Missing executable: $PROG (run ./build_all.sh first)" >&2
    exit 1
fi
for input in ned.sv ned_remap.sv ned_xor.sv ned_xor_remap.sv
do
    if [ ! -r "$input" ]; then
        echo "Missing serialized input: $input (run ./run_all_compr.sh first)" >&2
        exit 1
    fi
done

echo Running $PROG to benchmark compressed vector
echo ---------------------------------------------

echo Test default serialization
time "$PROG" -svin ned.sv -t -bench

echo ---------------------------------------------
echo Test char-set remapping
time "$PROG" -svin ned_remap.sv -t -bench

echo ---------------------------------------------
echo Test XOR compression
time "$PROG" -svin ned_xor.sv -t -bench

echo ---------------------------------------------
echo Test with XOR compression and char-set remapping
time "$PROG" -svin ned_xor_remap.sv -t -bench

echo
echo DONE

