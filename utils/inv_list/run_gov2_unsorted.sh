#!/bin/sh

echo
echo "BitMagic inverted collections experiments"
echo

echo "Compressing Sorted Gov2 set (regular build). Binary Interpolated (CM)"
./bminv_reg -t -s -u32in gov2.unsorted -bvout gov2.unsorted.bv.cm

echo
echo "Verify Sorted Gov2 set (regular build). Binary Interpolated (CM)"
./bminv_reg -t -s -verify -u32in gov2.unsorted -bvin gov2.unsorted.bv.cm

echo
echo "Compressing Sorted Gov2 set (SSE42). Binary Interpolated (CM)"
./bminv_sse42 -t -s -u32in gov2.unsorted -bvout gov2.unsorted.bv.cm

echo
echo "Verify Sorted Gov2 set (SSE42). Binary Interpolated (CM)"
./bminv_sse42 -t -s -verify -u32in gov2.unsorted -bvin gov2.unsorted.bv.cm

echo
echo "Compressing Sorted Gov2 set (AVX2). Binary Interpolated (CM)"
./bminv_avx2 -t -s -u32in gov2.unsorted -bvout gov2.unsorted.bv.cm

echo
echo "Verify Sorted Gov2 set (AVX2). Binary Interpolated (CM)"
./bminv_avx2 -t -s -verify -u32in gov2.unsorted -bvin gov2.unsorted.bv.cm



echo "Compressing Sorted Gov2 set (regular build). Elias Gamma."
./bminv_reg -t -s -level 4 -u32in gov2.unsorted -bvout gov2.unsorted.bv.eg

echo
echo "Verify Sorted Gov2 set (regular build). Elias Gamma."
./bminv_reg -t -s -verify -u32in gov2.unsorted -bvin gov2.unsorted.bv.eg

echo
echo "Compressing Sorted Gov2 set (SSE42). Elias Gamma."
./bminv_sse42 -t -s -level 4 -u32in gov2.unsorted -bvout gov2.unsorted.bv.eg

echo
echo "Verify Sorted Gov2 set (SSE42). Elias Gamma."
./bminv_sse42 -t -s -verify -u32in gov2.unsorted -bvin gov2.unsorted.bv.eg

echo
echo "Compressing Sorted Gov2 set (AVX2). Elias Gamma."
./bminv_avx2 -t -s -level 4 -u32in gov2.unsorted -bvout gov2.unsorted.bv.eg

echo
echo "Verify Sorted Gov2 set (AVX2). Elias Gamma."
./bminv_avx2 -t -s -verify -u32in gov2.unsorted -bvin gov2.unsorted.bv.cm

echo
ls -lap gov2.*
echo

