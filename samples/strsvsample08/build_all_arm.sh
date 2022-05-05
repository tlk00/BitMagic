#!/bin/sh
#
# Script for in-place build for regular realease builds, SSE4.2 and AVX2
#

cd ../../
. ./bmenv.sh
cd -

rm  ./strsvsample07
make clean

make BMOPTFLAGS=-DBMNEONOPT 
mv ./strsvsample07 ./strsvsample07_neon

make clean
make 
mv ./strsvsample07 ./strsvsample07_release

