#!/bin/sh
# Compile independent stress variants with isolated Make outputs and logs.
set -eu
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
exec "${PYTHON:-python3}" "$HERE/build_stress.py" "$@"
