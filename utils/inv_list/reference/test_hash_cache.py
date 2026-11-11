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
    flags = list(map(str, flags))
    mode_at = flags.index('--mode')
    mode = flags[mode_at + 1]
    del flags[mode_at:mode_at + 2]
    replacements = {'--candidate': '-bvin', '--candidate-index': '-bv-index-in',
                    '--raw': '-u32in', '--threads': '-threads', '--memory-mib': '-read-memory-mib',
                    '--batch-vectors': '-batch-vectors', '--max-vectors': '-max-vectors', '--mmap': '-read-mmap'}
    cache_at = flags.index('--cache')
    flags[cache_at] = '-hash-cache-in' if mode == 'cached' else '-hash-cache-out'
    flags = [replacements.get(flag, flag) for flag in flags]
    p = subprocess.run([str(a.bin.resolve()), *(['-verify'] if mode == 'cached' else []), *flags],
                       capture_output=True, text=True, timeout=90)
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
    common = ['--threads', 4, '--memory-mib', 1, '--batch-vectors', 3]
    caches = []
    for mode in ['cache-raw', 'cache-bv']:
        cache = root / mode
        source = ['--raw', raw] if mode == 'cache-raw' else ['--candidate', paths[0], '--candidate-index', paths[0].with_suffix('.idx')]
        invoke('--mode', mode, *source, '--cache', cache, *common)
        assert not Path(str(cache) + '.lock').exists()
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
    # Normal in-memory round trips exercise all levels with exact restored equality.
    roundtrip = [str(a.bin.resolve()), '-verify', '-roundtrip', '-threads', '4',
                 '-read-memory-mib', '4', '-batch-vectors', '3', '-bvin', str(paths[1]),
                 '-bv-index-in', str(paths[1].with_suffix('.idx'))]
    for flags in [[], ['-hash-cache-in', str(caches[0])],
                  *([['-read-mmap', '-hash-cache-in', str(caches[0])]] if sys.platform != 'win32' else [])]:
        result = subprocess.run([*roundtrip, *flags], capture_output=True, text=True, timeout=90)
        assert result.returncode == 0, (result.stdout, result.stderr)
        assert 'Roundtrip complete.' in result.stdout
        summary = next(line for line in result.stdout.splitlines() if line.startswith('PASS '))
        fields = dict(field.split('=', 1) for field in summary.split()[1:])
        assert int(fields['peak_staged_bytes']) <= 2 * 1024 * 1024
        assert 'Roundtrip buffer budget per worker=524288' in result.stdout
        assert f'Total round trips={len(rows) * 10}' in result.stdout
        levels = [line for line in result.stdout.splitlines() if line.startswith('Roundtrip level=')]
        assert len(levels) == 10
        for number, line in enumerate(levels, 1):
            fields = dict(field.split('=', 1) for field in line.split()[1:])
            assert int(fields['level']) == number and int(fields['vectors']) == len(rows)
    # Conservative output bounds are checked before allocating a worker BLOB.
    wide_raw, wide_blob = root / 'wide.raw', root / 'wide.bv'
    raw_file(wide_raw, [list(range(0, 16 * 65536, 2))])
    encode(wide_blob, 6, wide_raw)
    result = subprocess.run([str(a.bin.resolve()), '-verify', '-roundtrip', '-threads', '4',
                             '-read-memory-mib', '1', '-bvin', str(wide_blob),
                             '-bv-index-in', str(wide_blob.with_suffix('.idx'))],
                            capture_output=True, text=True, timeout=90)
    assert result.returncode == 1 and 'serialization bound exceeds worker budget' in result.stderr, result.stderr
    assert 'Roundtrip complete.' not in result.stdout
    # A correct SHA digest with a wrong cached POPCOUNT must fail the cheap check.
    wrong_count = root / 'wrong-count'
    content = bytearray(caches[0].read_bytes())
    struct.pack_into('<Q', content, 64 + 32, len(rows[0]) + 1)
    content[-32:] = hashlib.sha256(content[:-32]).digest()
    wrong_count.write_bytes(content)
    invoke('--mode', 'cached', *candidate, '--cache', wrong_count, *common,
           fail='Cached population count mismatch')
    # A concurrent writer lock must fail without touching the protected output.
    locked = root / 'locked-cache'
    lock = Path(str(locked) + '.lock')
    lock.mkdir()
    invoke('--mode', 'cache-raw', '--raw', raw, '--cache', locked, *common,
           fail='Hash cache writer lock exists')
    assert not locked.exists()
    lock.rmdir()
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
    assert not Path(str(failed_cache) + '.lock').exists()
    # Reject oversized vector lengths before payload allocation.
    raw.write_bytes(struct.pack('<I', 0xffffffff))
    invoke('--mode', 'cache-raw', '--raw', raw, '--cache', root / 'oversized', *common, fail='Vector exceeds staging budget')
print('PASS canonical oracle, raw/BV caches, all-level RAM round trips, buffer bounds, mmap and induced failures')
