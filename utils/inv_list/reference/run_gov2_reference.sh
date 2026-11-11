#!/usr/bin/env bash
# Build separately; this entry point only runs and validates the reference matrix.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$HERE/gov2_reference.py" "$@"
