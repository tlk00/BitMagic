#!/bin/sh
#
# Script for in-place build for regular realease builds, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./strsvsample09
make clean

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./strsvsample09 ./strsvsample09_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./strsvsample09 ./strsvsample09_avx2

make rebuild
mv ./strsvsample09 ./strsvsample09_release

