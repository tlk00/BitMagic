#!/bin/sh
#
# Run all benchmarking builds
#

echo "Regular build"
./strsvsample07_release -nodiag

echo "NEON build"
./strsvsample07_neon -nodiag

