#!/bin/sh
#
# Run all benchmarking builds
#

echo "Regular build"
./strsvsample09_release

echo "SSE4.2 build"
./strsvsample09_sse42

echo "AVX2 build"
./strsvsample09_avx2

