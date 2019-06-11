#!/bin/sh

#  build_all.sh
#  
#
#  Created by Anatoliy Kuznetsov on 9/2/17.
#

echo "DEBUG"
#./stress_debug || exit 1

echo
echo
echo RELEASE

./stress64_release || exit 1

echo
echo
echo SSE2

./stress64_release_sse2 || exit 1

echo
echo
echo SSE42

./stress64_release_sse42 || exit 1

echo
echo
echo AVX2

./stress64_release_avx2 || exit 1


echo
echo
echo 64-bit


./stress64_release_64 || exit 1

echo
echo "BitMagic Stress Test DONE"
