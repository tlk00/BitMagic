#!/bin/bash

#  run_all_tests.sh
#  
# run all tests in parallel
#

./run_all_parallel.sh sse42 || exit 1
./run_all_parallel.sh || exit 1
./run_all_parallel.sh avx2 || exit 1
./run_all_parallel.sh gcc || exit 1
./run_all_parallel.sh sse2 || exit 1
