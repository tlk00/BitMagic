#!/bin/sh
# Build native and architecture-specific variants. AVX-512 is opt-in.
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
    mv -f ./bminv "./$output"
}

# Preserve historical names used by existing scripts while exposing the same
# debug/release naming convention as the perf and stress build families.
build bminv_debug DEBUG=YES
build bminv_release
cp -f ./bminv_release ./bminv_reg
if [ "$platform" = arm ]; then
    build bminv_debug_neon DEBUG=YES BMOPTFLAGS=-DBMNEONOPT
    build bminv_release_neon BMOPTFLAGS=-DBMNEONOPT
    cp -f ./bminv_release_neon ./bminv_neon
else
    for kernel in sse2 sse42 avx2; do
        case "$kernel" in
            sse2) flags=-DBMSSE2OPT ;;
            sse42) flags=-DBMSSE42OPT ;;
            avx2) flags=-DBMAVX2OPT ;;
        esac
        build "bminv_debug_$kernel" DEBUG=YES "BMOPTFLAGS=$flags"
        build "bminv_release_$kernel" "BMOPTFLAGS=$flags"
        cp -f "./bminv_release_$kernel" "./bminv_$kernel"
    done
    if [ "$avx512" -eq 1 ]; then
        build bminv_debug_avx512 DEBUG=YES BMOPTFLAGS=-DBMAVX512OPT
        build bminv_release_avx512 BMOPTFLAGS=-DBMAVX512OPT
        cp -f ./bminv_release_avx512 ./bminv_avx512
    fi
fi
