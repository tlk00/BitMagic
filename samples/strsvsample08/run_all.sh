#!/bin/sh
#
# Run all benchmarking builds
#

echo "Regular build"
./strsvsample08_release -nodiag -bits 20

echo "SSE4.2 build"
./strsvsample08_sse42 -nodiag -bits 20

echo "AVX2 build"
./strsvsample08_avx2 -nodiag -bits 20

