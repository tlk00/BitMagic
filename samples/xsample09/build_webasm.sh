#!/bin/sh
emcc -std=c++11 -s ALLOW_MEMORY_GROWTH=1 -O2 -s WASM=1 -s DISABLE_EXCEPTION_CATCHING=0 -fno-rtti -I./../../src xsample09.cpp -o xsample09.html

