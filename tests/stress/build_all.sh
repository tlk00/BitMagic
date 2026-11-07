#!/bin/sh
# Build platform-specific test binaries. Optional: --avx512 on x86 Linux.
set -eu

case "${1:-}" in
    '') avx512=0 ;;
    --avx512) avx512=1 ;;
    -h|--help) echo "Usage: $0 [--avx512]"; exit 0 ;;
    *) echo "Usage: $0 [--avx512]" >&2; exit 2 ;;
esac
if [ "$#" -gt 1 ]; then
    echo "Usage: $0 [--avx512]" >&2
    exit 2
fi
cd "$(dirname "$0")"
os=$(uname -s)
arch=$(uname -m)
case "$os:$arch" in
    Darwin:arm64|Darwin:aarch64) platform=arm ;;
    Linux:x86_64|Linux:amd64) platform=x86 ;;
    *) echo "Unsupported platform: $os $arch" >&2; exit 1 ;;
esac
if [ "$platform" = arm ] && [ "$avx512" -eq 1 ]; then
    echo "--avx512 requires x86 Linux" >&2
    exit 2
fi

build()
{
    output=$1
    shift
    echo "Building $output"
    make "$@" rebuild
    mv -f ./bmtest "./$output"
}

build stress_debug DEBUG=YES
build stress_release
if [ "$platform" = arm ]; then
    build stress_debug_neon DEBUG=YES BMOPTFLAGS=-DBMNEONOPT
    build stress_release_neon BMOPTFLAGS=-DBMNEONOPT
else
    build stress_release_gcc BMOPTFLAGS=-DBM_USE_GCC_BUILD
    build stress_debug_sse2 DEBUG=YES BMOPTFLAGS=-DBMSSE2OPT
    build stress_release_sse2 BMOPTFLAGS=-DBMSSE2OPT
    build stress_debug_sse42 DEBUG=YES BMOPTFLAGS=-DBMSSE42OPT
    build stress_release_sse42 BMOPTFLAGS=-DBMSSE42OPT
    build stress_debug_avx2 DEBUG=YES BMOPTFLAGS=-DBMAVX2OPT
    build stress_release_avx2 BMOPTFLAGS=-DBMAVX2OPT
    if [ "$avx512" -eq 1 ]; then
        build stress_debug_avx512 DEBUG=YES BMOPTFLAGS=-DBMAVX512OPT
        build stress_release_avx512 BMOPTFLAGS=-DBMAVX512OPT
    fi
fi
