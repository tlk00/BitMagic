#!/bin/sh
#
# Script for in-place build for regular, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./xsample05_release ./xsample05 ./xsample05_sse42 ./xsample05_avx2
make clean

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./xsample05 ./xsample05_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./xsample05 ./xsample05_avx2

make rebuild
cp ./xsample05 ./xsample05_release

