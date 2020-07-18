#!/bin/sh

echo
echo "BitMagic inverted collections experiments"
echo

echo "Compressing Sorted Gov2 set (regular build). BIC-CM (6)"
./bminv_reg -t -s -level 6 -u32in gov2.sorted -bvout gov2.sorted.bv.cm6

echo
echo "Verify Sorted Gov2 set (regular build). BIC-CM (6)"
./bminv_reg -t -s -verify -u32in gov2.sorted -bvin gov2.sorted.bv.cm6

echo
echo "Compressing Sorted Gov2 set (SSE42). BIC-CM (6)"
./bminv_sse42 -t -s -level 6 -u32in gov2.sorted -bvout gov2.sorted.bv.cm6

echo
echo "Verify Sorted Gov2 set (SSE42). BIC-CM (6)"
./bminv_sse42 -t -s -verify -u32in gov2.sorted -bvin gov2.sorted.bv.cm6

echo
echo "Compressing Sorted Gov2 set (AVX2). BIC-CM (6)"
./bminv_avx2 -t -s -level 6 -u32in gov2.sorted -bvout gov2.sorted.bv.cm6

echo
echo "Verify Sorted Gov2 set (AVX2). BIC-CM (6)"
./bminv_avx2 -t -s -verify -u32in gov2.sorted -bvin gov2.sorted.bv.cm6



echo "Compressing Sorted Gov2 set (regular build). BIC-CM (5)"
./bminv_reg -t -s -level 5 -u32in gov2.sorted -bvout gov2.sorted.bv.cm5

echo
echo "Verify Sorted Gov2 set (regular build). BIC-CM (5)"
./bminv_reg -t -s -verify -u32in gov2.sorted -bvin gov2.sorted.bv.cm5

echo
echo "Compressing Sorted Gov2 set (SSE42). BIC-CM (5)"
./bminv_sse42 -t -s -level 5 -u32in gov2.sorted -bvout gov2.sorted.bv.cm5

echo
echo "Verify Sorted Gov2 set (SSE42). BIC-CM (5)"
./bminv_sse42 -t -s -verify -u32in gov2.sorted -bvin gov2.sorted.bv.cm5

echo
echo "Compressing Sorted Gov2 set (AVX2). BIC-CM (5)"
./bminv_avx2 -t -s -level 5 -u32in gov2.sorted -bvout gov2.sorted.bv.cm5

echo
echo "Verify Sorted Gov2 set (AVX2). BIC-CM (5)"
./bminv_avx2 -t -s -verify -u32in gov2.sorted -bvin gov2.sorted.bv.cm5



echo "Compressing Sorted Gov2 set (regular build). Elias Gamma."
./bminv_reg -t -s -level 4 -u32in gov2.sorted -bvout gov2.sorted.bv.eg

echo
echo "Verify Sorted Gov2 set (regular build). Elias Gamma."
./bminv_reg -t -s -verify -u32in gov2.sorted -bvin gov2.sorted.bv.eg

echo
echo "Compressing Sorted Gov2 set (SSE42). Elias Gamma."
./bminv_sse42 -t -s -level 4 -u32in gov2.sorted -bvout gov2.sorted.bv.eg

echo
echo "Verify Sorted Gov2 set (SSE42). Elias Gamma."
./bminv_sse42 -t -s -verify -u32in gov2.sorted -bvin gov2.sorted.bv.eg

echo
echo "Compressing Sorted Gov2 set (AVX2). Elias Gamma."
./bminv_avx2 -t -s -level 4 -u32in gov2.sorted -bvout gov2.sorted.bv.eg

echo
echo "Verify Sorted Gov2 set (AVX2). Elias Gamma."
./bminv_avx2 -t -s -verify -u32in gov2.sorted -bvin gov2.sorted.bv.cm

echo
ls -lap gov2.*
echo

