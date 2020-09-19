#!/bin/sh
#
# Script to run benchmarks on vector compression
#

PROG=./xsample05_release
#PROG=./xsample05_sse42
#PROG=./xsample05_avx2

echo Running $PROG to benchmark compressed vector
echo ---------------------------------------------

echo Test default serialization
time $PROG -svin ned.sv -t -bench

echo ---------------------------------------------
echo Test char-set remapping
time $PROG -svin ned_remap.sv -t -bench

echo ---------------------------------------------
echo Test XOR compression
time $PROG -svin ned_xor.sv -t -bench

echo ---------------------------------------------
echo Test with XOR compression and char-set remapping
time $PROG -svin ned_xor_remap.sv -t -bench

echo
echo DONE

