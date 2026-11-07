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
    mv -f ./strsvsample08 "./$output"
}


if [ "$platform" = arm ]; then
    build strsvsample08_neon BMOPTFLAGS=-DBMNEONOPT
else
    build strsvsample08_sse2 BMOPTFLAGS=-DBMSSE2OPT
    build strsvsample08_sse42 BMOPTFLAGS=-DBMSSE42OPT
    build strsvsample08_avx2 BMOPTFLAGS=-DBMAVX2OPT
    if [ "$avx512" -eq 1 ]; then
        build strsvsample08_avx512 BMOPTFLAGS=-DBMAVX512OPT
    fi
fi
# Build the scalar executable last for callers that use its ordinary name.
echo "Building strsvsample08_release"
make rebuild
mv -f ./strsvsample08 ./strsvsample08_release
