#!/bin/sh
emcc -std=c++17 -s ALLOW_MEMORY_GROWTH=1 -O2 -msse4.2 -msimd128 -D BMWASMSIMDOPT -s WASM=1 -s DISABLE_EXCEPTION_CATCHING=0  -fno-rtti -I./../../src perf.cpp -o perf.html

#emcc -std=c++11 -g3 --closure 0 --profiling --emit-symbol-map -s ALLOW_MEMORY_GROWTH=1 -O2 -msse4.2 -msimd128 -D BMWASMSIMDOPT -s WASM=1 -s DISABLE_EXCEPTION_CATCHING=0  -fno-rtti -I./../../src perf.cpp -o perf.html
#emcc -std=c++11 -s ALLOW_MEMORY_GROWTH=1 -O2 -s WASM=1 -s -fno-rtti -I./../../src perf.cpp -o perf.html
