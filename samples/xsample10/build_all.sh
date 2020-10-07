#!/bin/sh
#
# Script for in-place build for regular, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./xsample10
make clean

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./xsample10 ./xsample10_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./xsample10 ./xsample10_avx2

make rebuild
cp ./xsample10 ./xsample10_release

