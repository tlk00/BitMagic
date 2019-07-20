#!/bin/sh

#  build_all.sh
#  
#
#  Created by Anatoliy Kuznetsov on 9/2/17.
#

make clean

rm -rf ./bminv_* ./test
make rebuild
mv ./bminv ./bminv_reg

make BMOPTFLAGS=-DBMSSE42OPT rebuild
mv ./bminv ./bminv_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild
mv ./bminv ./bminv_avx2

