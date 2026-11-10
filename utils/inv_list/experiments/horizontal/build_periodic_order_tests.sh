#!/usr/bin/env bash
# Build only; benchmark execution and bminv builds are separate.
set -euo pipefail
[[ $# == 2 ]] || { echo "Usage: bash $0 BITMAGIC_ROOT BUILD_DIR" >&2; exit 2; }
REPO=$(cd "$1" && pwd)
mkdir -p "$2"
BUILD=$(cd "$2" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
CXX=${CXX:-c++}
for bits in 32 64; do
    DEFINES=()
    if [[ $bits == 64 ]]; then DEFINES=(-DTEST_BM64); fi
    "$CXX" -std=c++17 -O1 -g -fsanitize=address -fno-omit-frame-pointer \
        ${DEFINES[@]+"${DEFINES[@]}"} -I "$REPO/src" \
        "$HERE/periodic_order_test.cpp" -o "$BUILD/order-test$bits"
done
printf 'Built %s/order-test32 and order-test64; run them separately.\n' "$BUILD"
