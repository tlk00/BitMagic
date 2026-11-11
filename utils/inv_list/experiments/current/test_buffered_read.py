#!/usr/bin/env python3
"""Check indexed RAM batches, mapping windows, exact ID checks and bounded failures."""
import argparse
from array import array
from pathlib import Path
import re
import subprocess
import sys
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bin', required=True)
args = parser.parse_args()
binary = str(Path(args.bin).resolve())


def run(flags, error=None):
    result = subprocess.run([binary, *map(str, flags)], capture_output=True, text=True, timeout=30)
    text = result.stdout + result.stderr
    assert (result.returncode == 0) == (error is None), text
    if error:
        assert result.returncode == 1, text
        assert error in result.stderr and 'Error:' in result.stderr, text
        assert 'Verification complete.' not in text and 'Decode complete.' not in text, text
    return text


with tempfile.TemporaryDirectory(prefix='bm-buffered-read-') as directory:
    root = Path(directory)
    raw, blob, index = [root / name for name in ('raw.u32', 'data.bv', 'data.idx')]
    with raw.open('wb') as stream:
        for i in range(640):
            row = [n * 11 + i for n in range(1200 + i % 7 * 101)]
            array('I', [len(row), *row]).tofile(stream)
    run(['-s', '-threads', 4, '-u32in', raw, '-bvout', blob, '-bv-index-out', index])
    raw_data, encoded_data = raw.read_bytes(), blob.read_bytes()
    for mapped in [False, True] if sys.platform != 'win32' else [False]:
        flags = ['-read-mmap'] if mapped else []
        for workers, batch in [(1, 7), (4, 512)]:
            common = ['-s', '-threads', workers, '-batch-vectors', batch,
                      '-read-memory-mib', 1, '-bvin', blob, '-bv-index-in', index, *flags]
            text = run([*common, '-verify', '-u32in', raw])
            assert 'Verification complete.' in text
            peak = int(re.search(r'peak_staged_input_bytes=(\d+)', text).group(1))
            assert peak <= 1024 * 1024, peak
            assert 'Decode complete.' in run([*common, '-decode'])
        # Change an ID without changing cardinality/sortedness in first and later tasks.
        raw_offsets = []
        position = 0
        while position < len(raw_data):
            raw_offsets.append(position)
            length = int.from_bytes(raw_data[position:position + 4], sys.byteorder)
            position += 4 + length * 4
        blob_offsets = list(map(int, index.read_text().splitlines()[1:]))
        for failed_vector in [0, 639]:
            changed = bytearray(raw_data)
            position = raw_offsets[failed_vector] + 4
            value = int.from_bytes(changed[position:position + 4], sys.byteorder) + 1
            changed[position:position + 4] = array('I', [value]).tobytes()
            raw.write_bytes(changed)
            for reporting in [[], ['-progress-vectors', 1], ['-progress-vectors', 1, '-silent']]:
                text = run([*common, *reporting, '-verify', '-u32in', raw], 'Vector comparison failed')
                assert f'vector_index={failed_vector}' in text
                assert f'blob_begin={blob_offsets[failed_vector]}' in text
        raw.write_bytes(raw_data)
        # Reject untrusted lengths before allocating their advertised payload.
        raw.write_bytes(array('I', [0x40000000]).tobytes())
        run([*common, '-verify', '-u32in', raw], 'Vector exceeds read batch budget')
        raw.write_bytes(raw_data[:-1])
        run([*common, '-verify', '-u32in', raw], 'Truncated raw input vector')
        raw.write_bytes(raw_data)
        # One extra byte within the final indexed BLOB must fail exact consumption.
        original_index = index.read_text()
        lines = original_index.splitlines()
        header = lines[0].split()
        header[1] = str(int(header[1]) + 1)
        lines[0], lines[-1] = ' '.join(header), str(int(lines[-1]) + 1)
        index.write_text('\n'.join(lines) + '\n')
        blob.write_bytes(encoded_data + b'\0')
        run([*common, '-decode'], 'Indexed BLOB decode or length mismatch')
        index.write_text(original_index)
        blob.write_bytes(encoded_data)
        # A later task must also fail with progress enabled or explicitly silenced.
        last_begin = int(original_index.splitlines()[-2])
        truncated = encoded_data[:last_begin + 1]
        lines = original_index.splitlines()
        header = lines[0].split()
        header[1] = str(len(truncated))
        lines[0], lines[-1] = ' '.join(header), str(len(truncated))
        index.write_text('\n'.join(lines) + '\n')
        blob.write_bytes(truncated)
        for reporting in [[], ['-progress-vectors', 1], ['-progress-vectors', 1, '-silent']]:
            text = run([*common, *reporting, '-decode'], 'Decode worker failed')
            assert 'vector_index=639' in text and f'blob_begin={last_begin}' in text
            assert str(blob) in text
        index.write_text(original_index)
        blob.write_bytes(encoded_data)
        # Oversized indexed payload is rejected before loading/mapping the range.
        blob.write_bytes(b'')
        with blob.open('wb') as stream:
            stream.truncate(1024 * 1024 + 1)
        index.write_text('BM_INV_BLOB_INDEX_V1 1048577 1\n0\n1048577\n')
        run([*common, '-decode'], 'Vector exceeds read batch budget')
        index.write_text(original_index)
        blob.write_bytes(encoded_data)
    run([*common, '-read-memory-mib', 16384, '-decode'])
    for value in ['0', str(2**64 - 1), '-1', 'oops']:
        run(['-read-memory-mib', value], '-read-memory-mib')
    run(['-read-memory-mib', 1, '-bvin', blob, '-decode'], 'Read options require')
print('Buffered/mapped indexed reads, batch boundaries, exact IDs, truncation and memory limits PASS')
