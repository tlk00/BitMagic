#!/bin/sh
#
# Run all benchmarking builds
#

echo "Regular build"
./strsvsample07_release -nodiag

echo "SSE4.2 build"
./strsvsample07_sse42 -nodiag

echo "AVX2 build"
./strsvsample07_avx2 -nodiag

