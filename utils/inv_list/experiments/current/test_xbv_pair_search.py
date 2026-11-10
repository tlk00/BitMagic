#!/usr/bin/env python3
"""Compare deterministic pair-search CSVs/reports against a reference binary."""
from array import array
import argparse
from pathlib import Path
import random
import re
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bin', required=True)
parser.add_argument('--reference-bin', required=True)
args = parser.parse_args()
executables = [str(Path(p).resolve()) for p in (args.reference_bin, args.bin)]


def run(binary, arguments):
    result = subprocess.run([binary, *map(str, arguments)], capture_output=True, text=True)
    assert result.returncode == 0, (arguments, result.stdout, result.stderr)


with tempfile.TemporaryDirectory(prefix='bm-pair-search-') as directory:
    root = Path(directory)
    rng = random.Random(42)
    base = sorted(rng.sample(range(65536), 2000))
    rows = [base, base[:-2], sorted(set(base) | {65534}), base[::2],
            list(range(100, 2500)), list(range(101, 2501)),
            [k * 1009 + 3 for k in range(100)], [1, 65535]]
    raw = root / 'input.u32'
    with raw.open('wb') as stream:
        for row in rows:
            array('I', [len(row), *row]).tofile(stream)
    source = root / 'input.bv'
    run(executables[0], ['-s', '-level', 6, '-u32in', raw, '-bvout', source])
    cases = [(1, 0, 0, False, 0), (4, 0, 0, True, 0),
             (1, 2, 0, True, 0), (4, 1, 0, False, 4),
             (4, 0, 100000000, True, 0)]
    for number, (threads, ratio, gain, inspect, limit) in enumerate(cases):
        paths = []
        for version, binary in enumerate(executables):
            prefix = root / f'{number}-{version}.csv'
            flags = ['-s', '-pair-search', '-bvin', source, '-csv', prefix,
                     '-pair-threads', threads, '-pair-max-ratio', ratio,
                     '-pair-min-gain', gain, '-pair-min-pct', 0, '-pair-limit', limit]
            if inspect:
                flags.append('-pair-inspect')
            run(binary, flags)
            paths.append(prefix)
        suffixes = ['', '.best.csv', '.report.txt'] + (['.measured.csv'] if inspect else [])
        for suffix in suffixes:
            values = [Path(str(p) + suffix).read_text() for p in paths]
            if suffix == '.report.txt':
                values = [re.sub(r'Elapsed seconds=[^\n]*', 'Elapsed seconds=<time>', v) for v in values]
            assert values[0] == values[1], (number, suffix)
print('Pair search: CSV/report parity across inspect, ratio filters, prefix, empty results and 1/4 workers PASS')
