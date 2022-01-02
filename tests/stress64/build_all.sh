#!/bin/sh

#  build_all.sh
#  
#
#  Created by Anatoliy Kuznetsov on 9/2/17.
#

make clean

rm -rf ./stress_* ./test
make DEBUG=YES rebuild || exit 1
mv ./bmtest64 ./stress64_debug

make rebuild || exit 1
mv ./bmtest64 ./stress64_release

make BMOPTFLAGS=-DBMSSE2OPT DEBUG=YES rebuild || exit 1
mv ./bmtest64 ./stress64_debug_sse2

make BMOPTFLAGS=-DBMSSE2OPT rebuild || exit 1
mv ./bmtest64 ./stress64_release_sse2

make BMOPTFLAGS=-DBMSSE42OPT rebuild || exit 1
mv ./bmtest64 ./stress64_release_sse42

make BMOPTFLAGS=-DBMSSE42OPT DEBUG=YES rebuild || exit 1
mv ./bmtest64 ./stress64_debug_sse42

make BMOPTFLAGS=-DBMAVX2OPT rebuild || exit 1
mv ./bmtest64 ./stress64_release_avx2

make BMOPTFLAGS=-DBMAVX2OPT DEBUG=YES rebuild || exit 1
mv ./bmtest64 ./stress64_debug_avx2

#make BMOPTFLAGS=-DBMAVX512OPT rebuild
#mv ./test64 ./stress64_release_avx512

#make BMOPTFLAGS=-DBMAVX512OPT DEBUG=YES rebuild
#mv ./test64 ./stress64_debug_avx512


make BMOPTFLAGS=-DBM64OPT rebuild || exit 1
mv ./bmtest64 ./stress64_release_64
