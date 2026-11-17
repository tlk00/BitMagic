#!/usr/bin/env python3
"""Focused compressed sidecar tests, including 64-bit values in a 32-bit BM build."""
from array import array
from pathlib import Path
import subprocess
import sys
import tempfile

from index_io import index_tool, read_offsets, write_offsets


def main():
    """Exercise levels 0–10, malformed inputs, indexed reads and threaded writes."""
    binary = Path(sys.argv[1]).resolve()
    tool = index_tool()

    def run(*args, good=True):
        """Require the requested exit outcome for one bounded fixture command."""
        result = subprocess.run(list(map(str, args)), capture_output=True)
        assert (result.returncode == 0) == good, (args, result.stdout, result.stderr)
        return result

    with tempfile.TemporaryDirectory(prefix='bm-index-test-') as temporary:
        root = Path(temporary)
        index = root / 'data.idx'
        # Both cumulative offsets and individual lengths exceed the 32-bit range.
        offsets = [0, 2**32 + 9, 2**63 + 17, 2**64 - 1]
        for level in range(0, 11):
            write_offsets(index, offsets, level)
            assert read_offsets(index, offsets[-1]) == offsets
        payload = index.read_bytes()
        header, body = payload.split(b'\n', 1)
        for damaged in (payload[:-1], payload + b'x', header + b'\n',
                        header.replace(b' 3 ', b' 4 ') + b'\n' + body,
                        header.rsplit(b' ', 1)[0] + b' 11\n' + body):
            index.write_bytes(damaged)
            run(tool, 'dump', index, offsets[-1], good=False)
        index.write_bytes(payload)
        run(tool, 'dump', index, offsets[-1] - 1, good=False)
        # Smaller claimed size reaches checked cumulative-length arithmetic.
        fields = header.split()
        fields[1] = b'8'
        index.write_bytes(b' '.join(fields) + b'\n' + body)
        run(tool, 'dump', index, 8, good=False)
        index.write_bytes(b'BM_INV_BLOB_INDEX_V1 8 2\n0\n4\n8\n')
        run(tool, 'dump', index, 8, good=False)
        raw, blob = root / 'raw.u32', root / 'data.bv'
        with raw.open('wb') as stream:
            for values in ([1, 7, 90000], list(range(100, 3000)), [5, 65535]):
                array('I', [len(values), *values]).tofile(stream)
        for level in range(0, 11):
            run(binary, '-s', '-threads', 2, '-level', level, '-u32in', raw,
                '-bvout', blob, '-bv-index-out', index)
            assert index.read_bytes().split(b'\n', 1)[0].endswith(f' {level}'.encode())
            boundaries = read_offsets(index, blob.stat().st_size)
            copied = root / 'copy.idx'
            write_offsets(copied, boundaries, level)
            assert read_offsets(copied, blob.stat().st_size) == boundaries
            run(tool, 'check', index, blob)
            run(binary, '-s', '-verify', '-threads', 2, '-u32in', raw,
                '-bvin', blob, '-bv-index-in', index)
        # Valid legacy boundaries must be rejected specifically as an old format.
        retired = root / 'retired.idx'
        retired.write_text(f'BM_INV_BLOB_INDEX_V1 {boundaries[-1]} {len(boundaries)-1}\n' +
                           ''.join(f'{value}\n' for value in boundaries))
        rejected = run(binary, '-s', '-decode', '-threads', 2, '-bvin', blob,
                       '-bv-index-in', retired, good=False)
        assert b'compressed BLOB index header' in rejected.stderr
        print('PASS: compressed indexes, levels 0–10, indexed reads, 64-bit lengths, corruption and legacy rejection')


if __name__ == '__main__':
    main()
