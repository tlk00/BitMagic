#!/bin/sh
#

#CHECK=-check

echo "-----------------------------------------------"
echo Regular build:
echo
echo "NO bookmarks "
./xsample10_release -bookm 0 $CHECK
echo "Bookmarks:"
./xsample10_release -bookm 128 $CHECK
echo
./xsample10_release -bookm 64 $CHECK
echo
./xsample10_release -bookm 32 $CHECK
echo
./xsample10_release -bookm 16 $CHECK

echo "-----------------------------------------------"
echo                                          SSE4.2
echo
echo "NO bookmarks"
./xsample10_sse42 -bookm 0 $CHECK
echo "Bookmarks:"
./xsample10_sse42 -bookm 128 $CHECK
echo
./xsample10_sse42 -bookm 64 $CHECK
echo
./xsample10_sse42 -bookm 32 $CHECK
echo
./xsample10_sse42 -bookm 16 $CHECK

echo "-----------------------------------------------"
echo                                             AVX2
echo
echo "NO bookmarks"
./xsample10_avx2 -bookm 0 $CHECK
echo "Bookmarks:"
./xsample10_avx2 -bookm 128 $CHECK
echo
./xsample10_avx2 -bookm 64 $CHECK
echo
./xsample10_avx2 -bookm 32 $CHECK
echo
./xsample10_avx2 -bookm 16 $CHECK


