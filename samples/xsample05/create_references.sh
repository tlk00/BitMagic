#!/usr/bin/env bash
# Create or resume the xsample05 compression reference matrix.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$here/reference_matrix.py" create "$@"
