#!/bin/sh
# Build four serialization variants from the NED catalog report.
set -eu

cd "$(dirname "$0")"
case "$(uname -s):$(uname -m)" in
    Darwin:arm64|Darwin:aarch64) default_prog=./xsample05_neon ;;
    Linux:x86_64|Linux:amd64) default_prog=./xsample05_release ;;
    *) default_prog=./xsample05_release ;;
esac
PROG=${PROG:-$default_prog}
DICT=${DICT:-ned.txt}

if [ ! -x "$PROG" ]; then
    echo "Missing executable: $PROG (run ./build_all.sh first)" >&2
    exit 1
fi
if [ ! -r "$DICT" ]; then
    echo "Missing NED report: $DICT (see readme)" >&2
    exit 1
fi

echo Running $PROG to construct sorted compressed vector
echo

echo Default serialization
time "$PROG" -idict "$DICT" -svout ned.sv -t

echo
echo Serialization with char-set remapping
time "$PROG" -remap -idict "$DICT" -svout ned_remap.sv -t

echo Serialization with XOR compression
time "$PROG" -xor -idict "$DICT" -svout ned_xor.sv -t

echo
echo Serialization with XOR compression and char-set remapping
time "$PROG" -xor -remap -idict "$DICT" -svout ned_xor_remap.sv -t

echo
echo DONE

