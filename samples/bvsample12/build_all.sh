#!/bin/sh
#
# Script for in-place build for regular, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./sample12_* ./sample12
make clean

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./sample12 ./sample12_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./sample12 ./sample12_avx2

make rebuild
cp ./xsample12 ./sample12_release

