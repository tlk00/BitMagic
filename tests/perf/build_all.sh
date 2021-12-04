#!/bin/sh

#  build_all.sh
#  
#
#  Created by Anatoliy Kuznetsov on 9/2/17.
#

make clean

rm -rf ./perf_*
make DEBUG=YES rebuild || exit 1
mv ./bmperf ./perf_debug

make rebuild
mv ./bmperf ./perf_release || exit 1

make BMOPTFLAGS=-DBMSSE2OPT rebuild || exit 1
mv ./bmperf ./perf_release_sse2

make BMOPTFLAGS=-DBMSSE42OPT rebuild || exit 1
mv ./bmperf ./perf_release_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild || exit 1
mv ./bmperf ./perf_release_avx2

make BMOPTFLAGS=-DBMAVX512OPT rebuild || exit 1
mv ./bmperf ./perf_release_avx512

make BMOPTFLAGS=-DBM64OPT rebuild || exit 1
mv ./bmperf ./perf_release_64
