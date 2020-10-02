#!/bin/sh
#
# Script for in-place build for regular, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./xsample09
make clean

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./xsample09 ./xsample09_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./xsample09 ./xsample09_avx2

make rebuild
cp ./xsample09 ./xsample09_release

