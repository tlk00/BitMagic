#!/bin/sh
#
#  build_all.sh
#  

make clean

rm -rf ./ptest_*
make DEBUG=YES rebuild || exit 1
mv ./ptest ./ptest_debug

make rebuild || exit 1
mv ./ptest ./ptest_release

make BMOPTFLAGS=-DBMSSE42OPT rebuild || exit 1
mv ./ptest ./ptest_release_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild || exit 1
mv ./ptest ./ptest_release_avx2
