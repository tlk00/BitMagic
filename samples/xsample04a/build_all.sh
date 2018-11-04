#!/bin/sh
#
# Script for in-place build for regular, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./xsample04a_* ./xsample04
make clean

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./xsample04a ./xsample04a_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./xsample04a ./xsample04a_avx2

make rebuild
cp ./xsample04a ./xsample04a_release

