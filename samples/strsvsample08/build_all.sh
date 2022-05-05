#!/bin/sh
#
# Script for in-place build for regular realease builds, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./strsvsample08
make clean

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./strsvsample08 ./strsvsample08_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./strsvsample08 ./strsvsample08_avx2

make rebuild
mv ./strsvsample08 ./strsvsample08_release

