#!/usr/bin/env python3
"""Focused CLI and serialization regression checks; optional pre-cleanup binary."""
import argparse
from array import array
from pathlib import Path
import random
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bin', required=True)
parser.add_argument('--reference-bin')
args = parser.parse_args()
binary = str(Path(args.bin).resolve())
reference = str(Path(args.reference_bin).resolve()) if args.reference_bin else None


def run(executable, arguments, success=True):
    result = subprocess.run([executable, *map(str, arguments)], capture_output=True, text=True)
    assert (result.returncode == 0) == success, (arguments, result.stdout, result.stderr)
    return result.stdout + result.stderr


with tempfile.TemporaryDirectory(prefix='bm-cli-cleanup-') as directory:
    root = Path(directory)
    raw = root / 'input.u32'
    rng = random.Random(31)
    rows = [sorted(rng.sample(range(65536), 2500)), list(range(17, 131089)),
            [k * 1031 + 7 for k in range(96)], list(range(100, 800))]
    with raw.open('wb') as stream:
        for row in rows:
            array('I', [len(row), *row]).tofile(stream)
    assert 'Compression' in run(binary, ['--help'])
    for option in ['-unknown-option', '-level', '-bvout', '-u32in', '-threads']:
        run(binary, [option], False)
    # Retired experiment switches must fail before touching an output file.
    retired = ['-bwt', '-bwt-full', '-bwt-diag', '-bwt-sample',
               '-periodic', '-periodic-order', '-periodic-budget',
               '-periodic-diag', '-periodic-min-step', '-periodic-step',
               '-periodic-variants', '-periodic-min-length', '-periodic-holes',
               '-periodic-max-chains', '-periodic-sample', '-periodic-formula-chains',
               '-periodic-step1-start-reduce', '-periodic-step1-first-start',
               '-block-select', '-block-select-compose', '-block-select-multistep',
               '-multistep-max-step', '-multistep-min-gain', '-multistep-streak',
               '-multistep-min-stop', '-diag', '--diag', '-d', '--d']
    for option in retired:
        destination = root / 'preserve.bv'
        destination.write_bytes(b'preserve')
        assert 'Unknown option:' in run(binary, ['-u32in', raw, '-bvout', destination, option], False)
        assert destination.read_bytes() == b'preserve'
    for option in ['-svout', '--svout', '-svin', '--svin']:
        destination = root / 'preserve.sv'
        destination.write_bytes(b'preserve')
        assert 'Retired SV option' in run(binary, [option, destination], False)
        assert destination.read_bytes() == b'preserve'
    for level in [4, 6, 7, 8, 9, 10]:
        for workers in [0, 1, 10]:
            output = root / f'{level}-{workers}.bv'
            flags = ['-s', '-level', level, '-u32in', raw]
            threaded = ['-threads', workers, '-batch-vectors', 2] if workers else []
            run(binary, flags + threaded + ['-bvout', output])
            if reference:
                baseline = root / 'baseline.bv'
                run(reference, flags + threaded + ['-bvout', baseline])
                assert output.read_bytes() == baseline.read_bytes(), (level, workers)
            # Mode switches deliberately placed last; serial readers cover both output paths.
            run(binary, ['-s', '-u32in', raw, '-bvin', output, '-verify'])
            run(binary, ['-s', '-u32in', raw, '-bvin', output, '-v'])
            run(binary, ['-s', '-bvin', output, '-decode'])
    # Retained diagnostics/comparison policies still use the old explicit setup.
    for number, extra in enumerate([
            [], ['-gap-xor'], ['-gap-route-stats'], ['-sblock-adaptive'],
            ['-sblock-adaptive', '-gap-xor', '-sblock-wave-adaptive',
             '-sblock-large', '-block-endpoint-select', '-bwt-shared'],
            ['-sblock-correction-select']]):
        output = root / f'policy-{number}.bv'
        index = root / f'policy-{number}.idx'
        flags = ['-s', '-level', 6, '-threads', 2, '-u32in', raw,
                 '-block-select-step2', '-periodic-step1-interleaved', 4,
                 '-bwt-min-gain', 16, '-profile-opcodes', *extra]
        run(binary, flags + ['-bvout', output, '-bv-index-out', index])
        if reference:
            baseline = root / 'baseline.bv'
            run(reference, flags + ['-bvout', baseline])
            assert output.read_bytes() == baseline.read_bytes(), extra
        run(binary, ['-s', '-threads', 2, '-u32in', raw, '-bvin', output,
                     '-bv-index-in', index, '-verify'])
        run(binary, ['-s', '-threads', 2, '-bvin', output,
                     '-bv-index-in', index, '-decode'])
print('CLI rejection, final-position mode switches, 6 levels × serial/1/10 workers, retained policies, indexed round trips and optional byte identity: PASS')
