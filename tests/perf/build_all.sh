#!/bin/sh
# Build platform-specific test binaries. Optional: --avx512 on x86 Linux.
set -eu

avx512=0
pedantic=0
for arg in "$@"; do
    case "$arg" in
        --avx512) avx512=1 ;;
        --pedantic) pedantic=1 ;;
        -h|--help) echo "Usage: $0 [--avx512] [--pedantic]"; exit 0 ;;
        *) echo "Unknown option: $arg" >&2; exit 2 ;;
    esac
done
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
    mv -f ./bmperf "./$output"
}

build perf_debug DEBUG=YES
build perf_release
if [ "$platform" = arm ]; then
    build perf_debug_neon DEBUG=YES BMOPTFLAGS=-DBMNEONOPT
    build perf_release_neon BMOPTFLAGS=-DBMNEONOPT
else
    build perf_release_sse2 BMOPTFLAGS=-DBMSSE2OPT
    build perf_release_sse42 BMOPTFLAGS=-DBMSSE42OPT
    build perf_release_avx2 BMOPTFLAGS=-DBMAVX2OPT
    build perf_release_64 BMOPTFLAGS=-DBM64OPT
    if [ "$avx512" -eq 1 ]; then
        build perf_release_avx512 BMOPTFLAGS=-DBMAVX512OPT
    fi
fi

# Pedantic TAgent coverage includes debug versions of every x86 kernel.
if [ "$platform" = x86 ] && [ "$pedantic" -eq 1 ]; then
    build perf_debug_sse2 DEBUG=YES BMOPTFLAGS=-DBMSSE2OPT
    build perf_debug_sse42 DEBUG=YES BMOPTFLAGS=-DBMSSE42OPT
    build perf_debug_avx2 DEBUG=YES BMOPTFLAGS=-DBMAVX2OPT
    build perf_debug_64 DEBUG=YES BMOPTFLAGS=-DBM64OPT
    if [ "$avx512" -eq 1 ]; then
        build perf_debug_avx512 DEBUG=YES BMOPTFLAGS=-DBMAVX512OPT
    fi
fi
