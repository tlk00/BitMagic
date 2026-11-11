#!/usr/bin/env python3
"""Independent canonical-block oracle, cache equivalence, and induced failures."""
import argparse
from array import array
import hashlib
from pathlib import Path
import struct
import sys
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bin', required=True, type=Path)
parser.add_argument('--inv-bin', required=True, type=Path)
a = parser.parse_args()


def invoke(*flags, fail=None):
    """Run the experiment, require success or a contextual failure."""
    p = subprocess.run([str(a.bin.resolve()), *map(str, flags)], capture_output=True, text=True, timeout=90)
    assert (p.returncode == 0) == (fail is None), (p.stdout, p.stderr)
    if not fail:
        fields = dict(field.split('=', 1) for field in p.stdout.split() if '=' in field)
        assert int(fields['peak_staged_bytes']) <= int(fields['budget_bytes']), p.stdout
    if fail:
        assert fail in p.stderr, p.stderr
        assert 'PASS' not in p.stdout
    return p.stdout


def raw_file(path, rows):
    """Write the explicit LE collection used by the independent oracle."""
    with path.open('wb') as out:
        for row in rows:
            out.write(struct.pack('<I', len(row)))
            ids = array('I', row)
            if sys.byteorder != 'little':
                ids.byteswap()
            out.write(ids.tobytes())


def encode(path, level, raw):
    """Generate real indexed BLOBs through inv_list's ordinary entry point."""
    subprocess.run([str(a.inv_bin.resolve()), '-s', '-threads', '2', '-level', str(level),
                    '-u32in', str(raw), '-bvout', str(path), '-bv-index-out', str(path.with_suffix('.idx'))],
                   capture_output=True, check=True)


def oracle(row):
    """Hash block bytes independently using Python integer bit placement."""
    blocks = {}
    for value in row:
        number, bit = divmod(value, 65536)
        block = blocks.setdefault(number, bytearray(8192))
        block[bit // 8] |= 1 << (bit % 8)
    h = hashlib.sha256(b'BMBLOCKSHA256V1')
    for number, block in sorted(blocks.items()):
        h.update(block)
        h.update(struct.pack('<Q', number))
    return h.digest() + struct.pack('<Q', len(row))


with tempfile.TemporaryDirectory(prefix='block-sha-fixture-') as directory:
    root = Path(directory)
    rows = [[1, 3, 65539], list(range(65536, 131072)), [2**31 + 1, 2**31 + 3, 2**32 - 2],
            list(range(0, 65536, 2)), list(range(160000)), [1, 65537, 196609]] * 7
    raw = root / 'raw'
    raw_file(raw, rows)
    paths = [root / 'l6.bv', root / 'l10.bv']
    for path, level in zip(paths, [6, 10]):
        encode(path, level, raw)
    candidate = ['--candidate', paths[1], '--candidate-index', paths[1].with_suffix('.idx')]
    golden = ['--golden', paths[0], '--golden-index', paths[0].with_suffix('.idx')]
    common = ['--threads', 4, '--memory-mib', 1, '--batch-vectors', 3]
    for mode in ['exact', 'pair-sha']:
        for backend in [[], ['--mmap']] if sys.platform != 'win32' else [[]]:
            invoke('--mode', mode, *candidate, *golden, *common, *backend)
    caches = []
    for mode in ['cache-raw', 'cache-bv']:
        cache = root / mode
        source = ['--raw', raw] if mode == 'cache-raw' else ['--candidate', paths[0], '--candidate-index', paths[0].with_suffix('.idx')]
        invoke('--mode', mode, *source, '--cache', cache, *common)
        content = cache.read_bytes()
        assert content[:8] == b'BMSHA001'
        assert content[-32:] == hashlib.sha256(content[:-32]).digest()
        assert content[64:-32] == b''.join(map(oracle, rows))
        caches.append(cache)
        invoke('--mode', 'cached', *candidate, '--cache', cache, *common)
        if sys.platform != 'win32':
            invoke('--mode', 'cached', *candidate, '--cache', cache, *common, '--mmap')
        invoke('--mode', mode, *source, '--cache', cache, *common, fail='Refusing to overwrite')
    assert caches[0].read_bytes()[64:-32] == caches[1].read_bytes()[64:-32]
    # Independent raw oracle also covers empty vectors, which inv_list rejects on compression.
    empty_raw = root / 'empty-raw'
    raw_file(empty_raw, [[], [1, 3]])
    empty_cache = root / 'empty-cache'
    invoke('--mode', 'cache-raw', '--raw', empty_raw, '--cache', empty_cache, *common)
    assert empty_cache.read_bytes()[64:-32] == oracle([]) + oracle([1, 3])
    # No golden reads are possible after moving its BLOB and index away.
    paths[0].rename(root / 'hidden.bv')
    paths[0].with_suffix('.idx').rename(root / 'hidden.idx')
    invoke('--mode', 'cached', *candidate, '--cache', caches[0], *common)
    # Same cardinality, valid encoding, later vector substitution.
    changed = [row[:] for row in rows]
    changed[-1][0] = 2
    raw_file(raw, changed)
    encode(paths[1], 10, raw)
    invoke('--mode', 'cached', *candidate, '--cache', caches[0], *common, fail='vector=' + str(len(rows)-1))
    for mode in ['exact', 'pair-sha']:
        invoke('--mode', mode, *candidate, '--golden', root / 'hidden.bv',
               '--golden-index', root / 'hidden.idx', *common, fail='vector=' + str(len(rows)-1))
    # Exact consumption catches a valid BLOB followed by extra indexed bytes.
    saved_blob = paths[1].read_bytes()
    index_path = paths[1].with_suffix('.idx')
    saved_index = index_path.read_text()
    lines = saved_index.splitlines()
    header = lines[0].split()
    header[1] = str(int(header[1]) + 1)
    lines[0] = ' '.join(header)
    lines[-1] = str(int(lines[-1]) + 1)
    paths[1].write_bytes(saved_blob + b'\0')
    index_path.write_text('\n'.join(lines) + '\n')
    invoke('--mode', 'cached', *candidate, '--cache', caches[0], *common, fail='BLOB decode/length mismatch')
    paths[1].write_bytes(saved_blob)
    index_path.write_text(saved_index)
    # A one-record cache cannot validate a full candidate collection.
    prefix_cache = root / 'prefix'
    invoke('--mode', 'cache-bv', *candidate, '--cache', prefix_cache, *common, '--max-vectors', 1)
    invoke('--mode', 'cached', *candidate, '--cache', prefix_cache, *common, fail='Cache format/count/length mismatch')
    invoke('--mode', 'cached', *candidate, '--cache', prefix_cache, *common, '--max-vectors', 1)
    # Damaged cache, wrong count, and partial publication on raw read failure.
    damaged = root / 'damaged'
    data = bytearray(caches[0].read_bytes())
    data[64] ^= 1
    damaged.write_bytes(data)
    invoke('--mode', 'cached', *candidate, '--cache', damaged, *common, fail='Cache integrity mismatch')
    damaged.write_bytes(data[:-1])
    invoke('--mode', 'cached', *candidate, '--cache', damaged, *common, fail='Cache format/count/length mismatch')
    raw.write_bytes(raw.read_bytes()[:-1])
    failed_cache = root / 'failed-cache'
    invoke('--mode', 'cache-raw', '--raw', raw, '--cache', failed_cache, *common, fail='Truncated input')
    assert not failed_cache.exists()
    # Reject oversized vector lengths before payload allocation.
    raw.write_bytes(struct.pack('<I', 0xffffffff))
    invoke('--mode', 'cache-raw', '--raw', raw, '--cache', root / 'oversized', *common, fail='Vector exceeds staging budget')
print('PASS independent canonical oracle, exact/SHA comparisons, raw/BV caches, mmap, cache-only reads and induced failures')
