#!/bin/sh
#
# Script for sample benchmarking
#

THREADS=4
K_MIN=7
K_MAX=8
SIMD=_sse42

for ((K = $K_MIN ; K <= K_MAX ; K++));
do
  echo "----------------------------------------------------------------- $K"
  ./xsample07${SIMD} -kd test_${K}.kd -kdc test_${K}.kdc -fa NC_000001.11.fa -k $K  -t -j $THREADS | tee -a run_${K}.log
done

