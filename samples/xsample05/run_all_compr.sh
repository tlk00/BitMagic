cp ./run#!/bin/sh
#
# Script to run variants of vector compressions
#

#PROG=./xsample05_release
#PROG=./xsample05_sse42
PROG=./xsample05_avx2

echo Running $PROG to construct sorted compressed vector
echo

echo Default serialization
time $PROG -idict ned.txt -svout ned.sv -t

echo
echo Serialization with char-set remapping
time $PROG -remap -idict ned.txt -svout ned_remap.sv -t

echo Serialization with XOR compression
$PROG -xor -idict ned.txt -svout ned_xor.sv -t

echo
echo Serialization with XOR compression and char-set remapping
time $PROG -xor -remap -idict ned.txt -svout ned_xor_remap.sv -t

echo
echo DONE

