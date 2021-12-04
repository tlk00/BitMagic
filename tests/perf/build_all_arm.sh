#!/bin/sh

#  build_all.sh
#  
#
#  Created by Anatoliy Kuznetsov on 11/20/21.
#

make clean

rm -rf ./perf_*
rm -rf ./bmperf*

make DEBUG=YES || exit 1
mv ./bmperf ./perf_debug
make clean
make 
mv ./bmperf ./perf_release
make clean
make BMOPTFLAGS=-DBMNEONOPT || exit 1
mv ./bmperf ./perf_release_neon
make clean
