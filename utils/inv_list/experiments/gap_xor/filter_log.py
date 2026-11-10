#!/usr/bin/env python3
"""Keep collection statistics while omitting per-BIT-block trial lines."""
import sys

for line in sys.stdin:
    if not line.startswith(('Selection block ', 'BWT block ', 'Periodic block ')):
        sys.stdout.write(line)
