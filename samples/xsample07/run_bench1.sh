#!/bin/bash
#
# Script for sample benchmarking
#

THREADS=3
K_MIN=7
K_MAX=18

#SIMD=
#SIMD=_sse42
SIMD=_avx2
for ((K = $K_MIN ; K <= K_MAX ; K++));
do
  echo "----------------------------------------------------------------- $K"
  ./xsample07${SIMD} -kd test_${K}.kd -kdc test_${K}.kdc -fa NC_000001.11.fa -k $K  -t -j $THREADS | tee -a run_${K}${SIMD}.log
done

