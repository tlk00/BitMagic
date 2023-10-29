#!/bin/sh

#  build_all.sh
#  
#
#  Created by Anatoliy Kuznetsov on Oct-2023.
#

make clean

rm -rf ./stress_*
rm -rf ./bmtest*

make DEBUG=YES || exit 1
mv ./bmtest ./stress_debug
make clean
make 
mv ./bmtest ./stress_release
make clean
make BMOPTFLAGS=-DBMNEONOPT || exit 1
mv ./bmtest ./stress_release_neon
make clean
