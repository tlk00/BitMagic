#!/bin/sh

#  build_all.sh
#  
#
#  Created by Anatoliy Kuznetsov on 9/2/17.
#

make clean

rm -rf ./perf_*
make DEBUG=YES rebuild
mv ./perf ./perf_debug

make rebuild
mv ./perf ./perf_release

make BMOPTFLAGS=-DBMSSE2OPT rebuild
mv ./perf ./perf_release_sse2

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./perf ./perf_release_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./perf ./perf_release_avx2

make BMOPTFLAGS=-DBMAVX512OPT rebuild
mv ./perf ./perf_release_avx512

make BMOPTFLAGS=-DBM64OPT rebuild
mv ./perf ./perf_release_64
