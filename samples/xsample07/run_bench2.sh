#!/bin/bash
#
# Script for sample benchmarking
#


J_MIN=1
J_MAX=8

#SIMD=
#SIMD=_sse42
SIMD=_avx2

K=16

for ((THREADS = $J_MIN ; THREADS <= J_MAX ; THREADS+=1));
do
  echo "----------------------------------------------------------- ${THREADS}"
  ./xsample07${SIMD} -kd test_${K}.kd -kdc test_${K}.kdc -fa NC_000001.11.fa -k $K  -t -j $THREADS | tee -a run_${K}${SIMD}_p${THREADS}.log
done

