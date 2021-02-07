#!/bin/sh
#
#  build_all.sh
#  

make clean

rm -rf ./ptest_*
make DEBUG=YES rebuild
mv ./ptest ./ptest_debug

make rebuild
mv ./ptest ./ptest_release

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./ptest ./ptest_release_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./ptest ./ptest_release_avx2
