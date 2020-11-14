#!/bin/sh
emcc -std=c++11 -pthread -O2 -s USE_PTHREADS=1 -s TOTAL_MEMORY=268435456 -s PTHREAD_POOL_SIZE=4 -s WASM=1 -s DISABLE_EXCEPTION_CATCHING=0 -fno-rtti -I./../../src wasmtest.cpp -o wasmtest.html
#emcc -std=c++11 -s ALLOW_MEMORY_GROWTH=1 -s PTHREAD_POOL_SIZE=2  -O2 -s WASM=1 -s -fno-rtti -I./../../src perf.cpp -o perf.html
