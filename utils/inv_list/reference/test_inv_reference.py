#!/usr/bin/env python3
"""Fixture checks for fast/full reference tests, logical equality and crash propagation."""
from array import array
import argparse
import csv
import json
from pathlib import Path
import signal
import subprocess
import sys
import tempfile

from gov2_reference import digest

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--bin', required=True, type=Path)
parser.add_argument('--crash-test', action='store_true', help='Deliberately SIGSEGV a fixture process; macOS may display a crash dialog')
parser.add_argument('--generator', type=Path, default=HERE.parent / 'bminv_reference_neon')
args = parser.parse_args()
binary = args.bin.resolve()
generator = args.generator.resolve()


def invoke(flags, code=0):
    """Run the real harness and assert the required shell status."""
    result = subprocess.run(['bash', str(HERE.parent / 'run_inv_test_pl.sh'), *map(str, flags)],
                            capture_output=True, text=True)
    assert result.returncode == code, (result.returncode, result.stdout[-3000:], result.stderr)
    return result


def snapshot(root):
    """Capture every reference file's contents for independent immutability checks."""
    return {str(path.relative_to(root)): digest(path) for path in root.rglob('*') if path.is_file()}


with tempfile.TemporaryDirectory(prefix='bm-inv-reference-') as directory:
    root = Path(directory)
    reference = root / 'reference'
    inputs = reference / 'Gov2Flat'
    inputs.mkdir(parents=True)
    rows = [[1, 2, 3, 7, 1000], list(range(100, 2000)), [k*1009+3 for k in range(100)]]
    for dataset in ('sorted', 'unsorted'):
        with (inputs / ('gov2.' + dataset)).open('wb') as stream:
            for row in (rows if dataset == 'sorted' else rows[::-1]):
                array('I', [len(row), *row]).tofile(stream)
    result = subprocess.run([sys.executable, str(HERE / 'gov2_reference.py'), '--allow-local-output',
                             '--reserve-gib', '0', '--input-dir', str(inputs), '--out', str(reference),
                             '--bin', str(generator), '-j', '2'], capture_output=True, text=True)
    assert result.returncode == 0, (result.stdout[-2000:], result.stderr)
    before = snapshot(reference)
    common = ['--bin', binary, '--reference', reference, '-j', 2, '--reserve-gib', 0]
    for mode in ('fast', 'normal'):
        out = root / mode
        read_flags = ['--read-memory-mib', '1']
        if mode == 'fast' and sys.platform != 'win32':
            read_flags.append('--read-mmap')
        result = invoke([*common, *read_flags, '--out', out, *(['--fast', '--silent', '--progress', '1'] if mode == 'fast' else [])])
        if mode == 'fast':
            assert 'ALL PASS:' in result.stdout and 'verify 1v' not in result.stdout
            assert 'reference-verify-threaded' not in result.stdout
            assert 'Checking reference integrity:' not in result.stdout
            assert any('Progress verify ' in log.read_text() for log in out.glob('*/reference-verify-threaded.log'))
        records = list(csv.DictReader((out / 'summary.tsv').open(), delimiter='\t'))
        assert len(records) == 20 and {row['status'] for row in records} == {'passed'}
        assert (out / 'PASS.txt').exists()
        assert len(list(out.glob('*/*.bv'))) == (0 if mode == 'fast' else 20)
        assert snapshot(reference) == before
    # Fast tests the executable, without consulting checksum/history manifests.
    matrix_path = reference / 'matrix.json'
    original_matrix = matrix_path.read_bytes()
    matrix = json.loads(original_matrix)
    matrix['inputs']['sorted']['sha256'] = '0' * 64
    matrix_path.write_text(json.dumps(matrix))
    manifest_path = next((reference / 'reports').glob('*/metadata.json'))
    original_manifest = manifest_path.read_bytes()
    manifest_path.write_text('not a valid historical manifest')
    output = root / 'fast-with-stale-audit-records'
    result = invoke([*common, '--out', output, '--fast'])
    assert 'Checking reference integrity:' not in result.stdout
    assert 'Checking copied raw input:' not in result.stdout
    assert 'Rechecking reference integrity' not in result.stdout
    assert json.loads((output / 'test.json').read_text())['reference_integrity_audit'] is False
    failed = invoke([*common, '--out', root / 'normal-with-stale-audit-records'], 2)
    assert 'Copied raw input checksum mismatch' in failed.stderr
    matrix_path.write_bytes(original_matrix)
    manifest_path.write_bytes(original_manifest)
    assert snapshot(reference) == before
    invoke([*common, '--out', reference / 'forbidden'], 2)
    alias = root / 'reference-alias'
    alias.symlink_to(reference, target_is_directory=True)
    invoke([*common, '--out', alias / 'forbidden'], 2)
    failures = [('exit17', 'sys.exit(17)', 17)]
    if args.crash_test:
        failures.append(('crash', 'os.kill(os.getpid(), signal.SIGSEGV)', 128+signal.SIGSEGV))
    for name, statement, code in failures:
        wrapper = root / name
        wrapper.write_text('#!' + sys.executable + '\nimport os, sys, signal\n' + statement + '\n')
        wrapper.chmod(0o755)
        out = root / (name + '-output')
        result = invoke(['--bin', wrapper, '--reference', reference, '--out', out, '--fast',
                         '--reserve-gib', 0, '--silent'], code)
        assert 'reference-verify-threaded' in result.stderr and ('SIGSEGV' in result.stderr if name == 'crash' else 'exit 17' in result.stderr)
        assert (out / 'FAILURE.txt').exists() and not (out / 'PASS.txt').exists()
        assert snapshot(reference) == before
    # Simulate changed release encoding: valid level-6 streams for every requested level.
    wrapper = root / 'different-valid-encoding'
    wrapper.write_text('#!' + sys.executable + '\nimport os, sys\na=sys.argv[1:]\n'
                       'if "-bvout" in a: a[a.index("-level")+1]="6"\n'
                       'os.execv(' + repr(str(binary)) + ', [' + repr(str(binary)) + '] + a)\n')
    wrapper.chmod(0o755)
    out = root / 'different-encoding-output'
    invoke(['--bin', wrapper, '--reference', reference, '--out', out, '--reserve-gib', 0, '-j', 2])
    records = list(csv.DictReader((out / 'summary.tsv').open(), delimiter='\t'))
    assert all(row['status'] == 'passed' for row in records)
    assert any(row['serialized_bytes_identical'] == 'False' for row in records)
    assert snapshot(reference) == before
print('20 fast + 20 full + changed-encoding cases, raw-ID equality, exit propagation and read-only references PASS')
if args.crash_test:
    print('Deliberate SIGSEGV propagation: PASS')
