#!/bin/sh
#
# Script for in-place build for regular realease builds, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./strsvsample07
make clean

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./strsvsample07 ./strsvsample07_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./strsvsample07 ./strsvsample07_avx2

make rebuild
mv ./strsvsample07 ./strsvsample07_release

