#!/bin/sh

echo
echo "BitMagic inverted collections experiments"
echo

OUT_PATH="./"
#IN_PATH="./"
IN_PATH=/Volumes/WD-MacOS/int_set_bench/Gov2Flat/

echo "Compressing UnSorted Gov2 set (regular build). BIC-CM (6)"
./bminv_reg -t -s -level 6 -u32in ${IN_PATH}gov2.unsorted -bvout ${OUT_PATH}gov2.unsorted.bv.cm6

echo
echo "Verify UnSorted Gov2 set (regular build). BIC-CM (6)"
./bminv_reg -t -s -verify -u32in ${IN_PATH}gov2.unsorted -bvin ${OUT_PATH}gov2.unsorted.bv.cm6

echo
echo "Compressing UnSorted Gov2 set (SSE42). BIC-CM (6)"
./bminv_sse42 -t -s -level 6 -u32in ${IN_PATH}gov2.unsorted -bvout ${OUT_PATH}gov2.unsorted.bv.cm6

echo
echo "Verify UnSorted Gov2 set (SSE42). BIC-CM (6)"
./bminv_sse42 -t -s -verify -u32in ${IN_PATH}gov2.unsorted -bvin ${OUT_PATH}gov2.unsorted.bv.cm6

echo
echo "Compressing UnSorted Gov2 set (AVX2). BIC-CM (6)"
./bminv_avx2 -t -s -level 6 -u32in ${IN_PATH}gov2.unsorted -bvout ${OUT_PATH}gov2.unsorted.bv.cm6

echo
echo "Verify UnSorted Gov2 set (AVX2). BIC-CM (6)"
./bminv_avx2 -t -s -verify -u32in ${IN_PATH}gov2.unsorted -bvin ${OUT_PATH}gov2.unsorted.bv.cm6




echo "Compressing UnSorted Gov2 set (regular build). BIC-CM (5)"
./bminv_reg -t -s -level 5 -u32in ${IN_PATH}gov2.unsorted -bvout ${OUT_PATH}gov2.unsorted.bv.cm5

echo
echo "Verify UnSorted Gov2 set (regular build). BIC-CM (5)"
./bminv_reg -t -s -verify -u32in ${IN_PATH}gov2.unsorted -bvin ${OUT_PATH}gov2.unsorted.bv.cm5

echo
echo "Compressing UnSorted Gov2 set (SSE42). BIC-CM (5)"
./bminv_sse42 -t -s -level 5 -u32in ${IN_PATH}gov2.unsorted -bvout ${OUT_PATH}gov2.unsorted.bv.cm5

echo
echo "Verify UnSorted Gov2 set (SSE42). BIC-CM (5)"
./bminv_sse42 -t -s -verify -u32in ${IN_PATH}gov2.unsorted -bvin ${OUT_PATH}gov2.unsorted.bv.cm5

echo
echo "Compressing UnSorted Gov2 set (AVX2). BIC-CM (5)"
./bminv_avx2 -t -s -level 5 -u32in ${IN_PATH}gov2.unsorted -bvout ${OUT_PATH}gov2.unsorted.bv.cm5

echo
echo "Verify UnSorted Gov2 set (AVX2). BIC-CM (5)"
./bminv_avx2 -t -s -verify -u32in ${IN_PATH}gov2.unsorted -bvin ${OUT_PATH}gov2.unsorted.bv.cm5



echo "Compressing UnSorted Gov2 set (regular build). Elias Gamma."
./bminv_reg -t -s -level 4 -u32in ${IN_PATH}gov2.unsorted -bvout ${OUT_PATH}gov2.unsorted.bv.eg

echo
echo "Verify UnSorted Gov2 set (regular build). Elias Gamma."
./bminv_reg -t -s -verify -u32in ${IN_PATH}gov2.unsorted -bvin ${OUT_PATH}gov2.unsorted.bv.eg

echo
echo "Compressing UnSorted Gov2 set (SSE42). Elias Gamma."
./bminv_sse42 -t -s -level 4 -u32in ${IN_PATH}gov2.unsorted -bvout ${OUT_PATH}gov2.unsorted.bv.eg

echo
echo "Verify UnSorted Gov2 set (SSE42). Elias Gamma."
./bminv_sse42 -t -s -verify -u32in ${IN_PATH}gov2.unsorted -bvin ${OUT_PATH}gov2.unsorted.bv.eg

echo
echo "Compressing UnSorted Gov2 set (AVX2). Elias Gamma."
./bminv_avx2 -t -s -level 4 -u32in ${IN_PATH}gov2.unsorted -bvout ${OUT_PATH}gov2.unsorted.bv.eg

echo
echo "Verify UnSorted Gov2 set (AVX2). Elias Gamma."
./bminv_avx2 -t -s -verify -u32in ${IN_PATH}gov2.unsorted -bvin ${OUT_PATH}gov2.unsorted.bv.eg

echo
ls -lap gov2.*
echo

