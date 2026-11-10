#!/usr/bin/env python3
"""Summarize verified GAP-XOR collection experiment."""
import pathlib
import re
import sys


def index_count_and_size(path, count):
    with open(path, encoding='ascii') as inp:
        first = inp.readline().split()
        if len(first) != 3 or first[0] != 'BM_INV_BLOB_INDEX_V1':
            raise SystemExit(f'Bad BLOB index: {path}')
        total = int(first[2])
        if count == 0:
            count = total
        if count > total:
            raise SystemExit(f'Reference has {total} vectors; requested {count}')
        offsets = [int(inp.readline()) for _ in range(count + 1)]
        if offsets[0] != 0 or offsets != sorted(offsets):
            raise SystemExit(f'Bad offsets in {path}')
        return offsets[-1]


def elapsed(path):
    return float(path.read_text().split()[0])


def require_pass(out, mode, count):
    markers = {'compress': '', 'verify': 'Verification complete.',
               'decode': 'Decode complete.'}
    for phase in ('compress', 'verify', 'decode'):
        path = out / f'{mode}.{phase}.log'
        if not path.is_file():
            raise SystemExit(f'Missing {mode}.{phase}: {path}')
        content = path.read_text(errors='replace')
        if (not re.search(r'^Total vectors=' + str(count) + r'$', content, re.M) or
                markers[phase] not in content):
            raise SystemExit(f'Missing count or failed {mode}.{phase}: {path}')


def stats(text):
    matches = re.findall(r'GAP XOR ([^\n]+)', text)
    return '\n'.join('  ' + line for line in matches)


def total_in_index(path):
    with open(path, encoding='ascii') as inp:
        parts = inp.readline().split()
    if len(parts) != 3 or parts[0] != 'BM_INV_BLOB_INDEX_V1':
        raise SystemExit(f'Bad BLOB index: {path}')
    return int(parts[2])


def require_full_input(log):
    consumed = re.search(r'^Input consumed bytes=(\d+)$', log, re.M)
    total = re.search(r'^Input file bytes=(\d+)$', log, re.M)
    if not consumed or not total or consumed.group(1) != total.group(1):
        raise SystemExit('Compression did not consume the complete raw input')


if sys.argv[1] == 'check-reference':
    index_count_and_size(pathlib.Path(sys.argv[2]), int(sys.argv[3]))
elif sys.argv[1] == 'report':
    out, requested, ref = pathlib.Path(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    count = requested or total_in_index(out / 'candidate.idx')
    require_pass(out, 'candidate', count)
    if not requested:
        require_full_input((out / 'candidate.compress.log').read_text(errors='replace'))
    candidate_size = index_count_and_size(out / 'candidate.idx', count)
    if ref:
        if not requested and total_in_index(pathlib.Path(ref) / 'auto.idx') != count:
            raise SystemExit('Full candidate vector count differs from reference')
        control_size = index_count_and_size(pathlib.Path(ref) / 'auto.idx', count)
        control_time = None
    else:
        require_pass(out, 'control', count)
        if not requested:
            require_full_input((out / 'control.compress.log').read_text(errors='replace'))
        control_size = index_count_and_size(out / 'control.idx', count)
        control_time = elapsed(out / 'control.compress.time')
    gain = control_size - candidate_size
    print(f'Correctness: PASS; candidate verified and decoded ({count:,} vectors).')
    print(f'Control: {control_size:,} bytes ({"reused reference" if ref else "verified fresh stream"})')
    print(f'Candidate: {candidate_size:,} bytes')
    print(f'Gain: {gain:+,} bytes ({100 * gain / control_size:+.6f}%)')
    print(f'Candidate encode: {elapsed(out / "candidate.compress.time"):.3f}s; verify: '
          f'{elapsed(out / "candidate.verify.time"):.3f}s; decode: '
          f'{elapsed(out / "candidate.decode.time"):.3f}s')
    if control_time is not None:
        print(f'Control encode: {control_time:.3f}s; verify: '
              f'{elapsed(out / "control.verify.time"):.3f}s; decode: '
              f'{elapsed(out / "control.decode.time"):.3f}s')
    print(stats((out / 'candidate.compress.log').read_text(errors='replace')))
else:
    raise SystemExit('usage: report_gap_xor.py check-reference INDEX N | report OUT N REF')
