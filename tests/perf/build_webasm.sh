#!/bin/sh
emcc -std=c++17 -s ALLOW_MEMORY_GROWTH=1 -O2 -s WASM=1 -s DISABLE_EXCEPTION_CATCHING=0 -fno-rtti -I./../../src perf.cpp -o perf.html
#emcc -std=c++11 -s ALLOW_MEMORY_GROWTH=1 -O2 -s WASM=1 -s -fno-rtti -I./../../src perf.cpp -o perf.html
