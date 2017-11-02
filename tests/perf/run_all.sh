#!/bin/sh

#  build_all.sh
#  
#
#  Created by Anatoliy Kuznetsov on 9/2/17.
#

echo "DEBUG"
#./perf_debug || exit 1

echo
echo
echo RELEASE

./perf_release || exit 1

echo
echo
echo SSE2

./perf_release_sse2 || exit 1

echo
echo
echo SSE42

./perf_release_sse42 || exit 1

echo
echo
echo AVX2

./perf_release_avx2 || exit 1


echo
echo
echo 64-bit


./perf_release_64 || exit 1
echo
echo "BitMagic Perf Test DONE"
