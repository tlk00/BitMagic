#!/bin/sh

#  build_all_arm.sh
#  

make clean

rm -rf ./bminv_* ./test
make rebuild
mv ./bminv ./bminv_reg


