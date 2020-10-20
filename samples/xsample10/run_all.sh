#!/bin/sh
#


echo "-----------------------------------------------"
echo Regular build:
echo
echo "NO bookmarks"
./xsample10_release
echo "Bookmarks:"
./xsample10_release -bookm 128
echo
./xsample10_release -bookm 64
echo
./xsample10_release -bookm 32
echo
./xsample10_release -bookm 16

echo "-----------------------------------------------"
echo                                          SSE4.2
echo
echo "NO bookmarks"
./xsample10_sse42
echo "Bookmarks:"
./xsample10_sse42 -bookm 128
echo
./xsample10_sse42 -bookm 64
echo
./xsample10_sse42 -bookm 32
echo
./xsample10_sse42 -bookm 16

echo "-----------------------------------------------"
echo                                             AVX2
echo
echo "NO bookmarks"
./xsample10_avx2
echo "Bookmarks:"
./xsample10_avx2 -bookm 128
echo
./xsample10_avx2 -bookm 64
echo
./xsample10_avx2 -bookm 32
echo
./xsample10_avx2 -bookm 16


