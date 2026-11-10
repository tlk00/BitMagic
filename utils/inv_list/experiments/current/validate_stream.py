#!/usr/bin/env python3
"""Validate collection scope, raw verification, independent decode and index framing."""
from pathlib import Path
import re

def one(text, pattern):
    m = re.search(pattern, text)
    if not m:
        raise SystemExit(f'Missing result: {pattern}')
    return int(m.group(1))

def read(mode, root, limit):
    d = root / mode
    logs = {p: (d/f'{p}.log').read_text(errors='replace')
            for p in ('compress', 'verify', 'decode')}
    vectors = one(logs['compress'], r'Total vectors=(\d+)')
    ids = one(logs['compress'], r'Total ints=(\d+)')
    size = one(logs['compress'], r'BV size = (\d+)')
    if vectors == 0 or ids == 0 or (limit and vectors != limit):
        raise SystemExit(f'{mode}: missing vectors/IDs or vector-limit mismatch')
    for phase, marker in [('verify', 'Verification complete.'), ('decode', 'Decode complete.')]:
        if marker not in logs[phase] or one(logs[phase], r'Total vectors=(\d+)') != vectors:
            raise SystemExit(f'{mode}.{phase}: incomplete')
    if one(logs['verify'], r'Total ints=(\d+)') != ids:
        raise SystemExit('Verification ID count differs')
    if (d/'candidate.bv').stat().st_size != size:
        raise SystemExit('Stream size differs from log')
    lines = (d/'candidate.idx').read_text().splitlines()
    if lines[0].split() != ['BM_INV_BLOB_INDEX_V1', str(size), str(vectors)]:
        raise SystemExit('Index header mismatch')
    offsets = list(map(int, lines[1:]))
    if len(offsets) != vectors + 1 or offsets[0] != 0 or offsets[-1] != size or any(
            b <= a for a, b in zip(offsets, offsets[1:])):
        raise SystemExit('Invalid index offsets')
    timings = [float((d/f'{p}.time').read_text().split()[0])
               for p in ('compress', 'verify', 'decode')]
    return dict(vectors=vectors, ids=ids, size=size, timings=timings, log=logs['compress'])
