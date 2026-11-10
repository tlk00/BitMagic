#!/usr/bin/env python3
"""Exercise the portable runner with real serial/indexed streams and path spaces."""
import argparse
from array import array
from pathlib import Path
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bin', required=True)
args = parser.parse_args()
binary = str(Path(args.bin).resolve())
runner = Path(__file__).resolve().parents[2] / 'run_collection.sh'


def run(arguments, success=True):
    result = subprocess.run(['bash', str(runner), *map(str, arguments)],
                            capture_output=True, text=True)
    assert (result.returncode == 0) == success, (arguments, result.stdout, result.stderr)


with tempfile.TemporaryDirectory(prefix='bm runner ') as directory:
    root = Path(directory)
    raw = root / 'input lists.u32'
    with raw.open('wb') as stream:
        for row in [list(range(10, 200)), list(range(20, 1000, 3)), [2, 99999]]:
            array('I', [len(row), *row]).tofile(stream)
    base = ['--bin', binary, '--input', raw, '--levels', '6,10', '--max-vectors', '2']
    streams = []
    for workers in [0, 2]:
        output = root / f'result {workers}'
        flags = base + ['--out', output] + (['--threads', workers] if workers else [])
        run(flags)
        for level in [6, 10]:
            destination = output / f'level-{level}'
            assert 'Verification complete.' in (destination / 'verify.log').read_text()
            assert 'Decode complete.' in (destination / 'decode.log').read_text()
            assert 'Total vectors=2' in (destination / 'compress.log').read_text()
            assert (destination / 'collection.bv.idx').exists() == bool(workers)
        streams.append([(output / f'level-{level}/collection.bv').read_bytes() for level in [6, 10]])
        run(flags, False)
        assert not (output / 'FAILURE.txt').exists()  # Refusal must not modify prior results.
    assert streams[0] == streams[1]
    for option, value in [('--levels', '6,6'), ('--levels', '11'),
                          ('--threads', '257'), ('--progress', '0'), ('--max-vectors', '')]:
        output = root / 'invalid'
        run(base + ['--out', output, option, value], False)
        assert not output.exists()
print('Collection runner: serial/indexed round trips, prefixes, spaces, byte parity and refusal checks PASS')
