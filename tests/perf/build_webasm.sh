#!/bin/sh
emcc -std=c++11 -s ALLOW_MEMORY_GROWTH=1 -O2  --emrun -I./../../src perf.cpp -o perf.html
