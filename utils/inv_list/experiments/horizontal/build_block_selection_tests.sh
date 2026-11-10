#!/usr/bin/env bash
# Build only; execute the resulting tests separately.
set -euo pipefail
[[ $# == 2 ]] || { echo "Usage: bash $0 BITMAGIC_ROOT BUILD_DIR" >&2; exit 2; }
REPO=$(cd "$1" && pwd)
mkdir -p "$2"
BUILD=$(cd "$2" && pwd)
HERE=$(cd "$(dirname "$0")" && pwd)
for bits in 32 64; do
    DEFINES=()
    if [[ $bits == 64 ]]; then DEFINES=(-DTEST_BM64); fi
    "${CXX:-c++}" -std=c++17 -O1 -g -fsanitize=address -fno-omit-frame-pointer \
        ${DEFINES[@]+"${DEFINES[@]}"} -I "$REPO/src" \
        "$HERE/block_selection_test.cpp" -o "$BUILD/block-selection-test$bits"
done
printf 'Built %s/block-selection-test32 and block-selection-test64.\n' "$BUILD"
