#!/bin/sh
#


echo Regular build:
echo
./xsample10_release

echo
echo SSE4.2 build:
echo
./xsample10_sse42

echo
echo AVX2 build:
echo
./xsample10_avx2

