#!/bin/sh
#


echo Regular build:
echo
./xsample09_release

echo
echo SSE4.2 build:
echo
./xsample09_sse42

echo
echo AVX2 build:
echo
./xsample09_avx2

