#!/usr/bin/env bash
# Validate every selected xsample05 reference using the supplied reader.
set -euo pipefail
here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
exec python3 "$here/reference_matrix.py" validate "$@"
