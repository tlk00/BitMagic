#!/usr/bin/env bash
# Test a chosen platform binary against immutable Gov2 reference collections.
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$HERE/reference/inv_reference_test.py" "$@"
