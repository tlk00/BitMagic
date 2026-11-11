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
    for mode in ('fast', 'normal', 'pedantic'):
        out = root / mode
        read_flags = ['--read-memory-mib', '16384' if mode == 'fast' else '1']
        if mode == 'fast' and sys.platform != 'win32':
            read_flags.append('--read-mmap')
        result = invoke([*common, *read_flags, '--out', out, *(['--fast', '--silent', '--progress', '1'] if mode == 'fast' else (['--pedantic'] if mode == 'pedantic' else ['--silent']))])
        if mode == 'fast':
            assert 'ALL PASS:' in result.stdout and 'verify 1v' not in result.stdout
            assert 'reference-verify-threaded' not in result.stdout
            assert 'Checking reference integrity:' not in result.stdout
            assert any('Progress verify ' in log.read_text() for log in out.glob('*/golden-verify-threaded.log'))
        records = list(csv.DictReader((out / 'summary.tsv').open(), delimiter='\t'))
        assert len(records) == 20 and {row['status'] for row in records} == {'passed'}
        assert (out / 'PASS.txt').exists()
        assert not list(out.glob('*/*.bv'))
        if mode != 'fast':
            roundtrip_logs = list(out.glob('*/roundtrip-levels.log'))
            assert len(roundtrip_logs) == 2
            for phase_log in roundtrip_logs:
                assert 'Progress roundtrip completed=3 total=3' in phase_log.read_text()
            if mode == 'normal':
                assert 'Progress roundtrip' not in result.stdout
            assert len(list(out.glob('*/reference-verify-serial.log'))) == (20 if mode == 'pedantic' else 0)
        assert snapshot(reference) == before
    # Normal sorted-only selection must not accidentally validate unsorted.
    out = root / 'normal-sorted-only'
    invoke([*common, '--out', out, '--dataset', 'sorted', '--hash', '--hash-cache-dir', root / 'sorted-cache', '--start-level', 1])
    rows_selected = list(csv.DictReader((out / 'summary.tsv').open(), delimiter='\t'))
    assert len(rows_selected) == 10
    assert {row['dataset'] for row in rows_selected} == {'Gov2.sorted'}
    assert {int(row['level']) for row in rows_selected} == set(range(1, 11))
    metadata = json.loads((out / 'test.json').read_text())
    assert metadata['datasets'] == ['sorted'] and metadata['case_count'] == 10
    assert not list(out.glob('gov2.unsorted*'))

    # Level selection skips lower-level stream I/O in fast mode and counts both datasets.
    for mode, start_level in [('fast', 6), ('fast', 10), ('normal', 10)]:
        out = root / f'{mode}-start-{start_level}'
        flags = ['--fast'] if mode == 'fast' else []
        result = invoke([*common, '--out', out, '--start-level', start_level, *flags])
        records = list(csv.DictReader((out / 'summary.tsv').open(), delimiter='\t'))
        expected_cases = 2 * (11 - start_level)
        assert len(records) == expected_cases
        assert {(row['dataset'], int(row['level'])) for row in records} == {
            (f'Gov2.{dataset}', level) for dataset in ('sorted', 'unsorted')
            for level in range(start_level, 11)}
        assert {row['status'] for row in records} == {'passed'}
        assert f'ALL PASS: {expected_cases} cases' in result.stdout
        assert f'All {expected_cases} ' in (out / 'PASS.txt').read_text()
        metadata = json.loads((out / 'test.json').read_text())
        assert (metadata['start_level'], metadata['end_level'], metadata['case_count']) == (
            start_level, 10, expected_cases)
        assert not list(out.glob('gov2.*.level-01'))
        assert snapshot(reference) == before
    out = root / 'custom-golden'
    invoke([*common, '--out', out, '--fast', '--start-level', 10, '--golden-level', 3])
    assert json.loads((out / 'test.json').read_text())['golden_level'] == 3
    assert len(list(out.glob('*/golden-verify-threaded.log'))) == 2
    assert len(list(out.glob('*/reference-verify-threaded.log'))) == 2
    # A missing prerequisite baseline must fail even outside the selected range.
    missing = reference / 'gov2.sorted.level-06.bv'
    hidden = missing.with_suffix('.hidden')
    missing.rename(hidden)
    try:
        invoke([*common, '--out', root / 'missing-golden', '--fast', '--start-level', 10], 2)
        assert not (root / 'missing-golden/PASS.txt').exists()
    finally:
        hidden.rename(missing)
    lower_stream = reference / 'gov2.sorted.level-01.bv'
    unavailable = lower_stream.with_suffix('.unavailable')
    lower_stream.rename(unavailable)
    try:
        invoke([*common, '--out', root / 'skip-missing-lower-stream', '--fast', '--start-level', 6])
    finally:
        unavailable.rename(lower_stream)
    for option in ['--start-level', '--golden-level']:
        for value in ['0', '11', '-1', 'oops']:
            out = root / ('invalid-' + option[2:] + '-' + value)
            invoke([*common, '--out', out, '--fast', option, value], 2)
            assert not out.exists()
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
    failed = invoke([*common, '--pedantic', '--out', root / 'pedantic-with-stale-audit-records'], 2)
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
        assert 'golden-verify-threaded' in result.stderr and ('SIGSEGV' in result.stderr if name == 'crash' else 'exit 17' in result.stderr)
        assert (out / 'FAILURE.txt').exists() and not (out / 'PASS.txt').exists()
        assert snapshot(reference) == before
    # Real worker failures must survive the child/relay/harness boundary and be logged.
    selected = reference / 'gov2.sorted.level-10.bv'
    selected_index = selected.with_suffix('.idx')
    for failure in ['verify', 'deserialize', 'golden-deserialize', 'golden-mismatch']:
        out = root / ('worker-' + failure)
        golden_blob = reference / 'gov2.sorted.level-06.bv'
        saved_golden = golden_blob.read_bytes()
        saved_raw = (inputs / 'gov2.sorted').read_bytes()
        saved_blob, saved_index = selected.read_bytes(), selected_index.read_bytes()
        try:
            if failure in ('verify', 'golden-mismatch'):
                changed = bytearray(saved_raw)
                changed[4:8] = array('I', [0]).tobytes()
                (inputs / 'gov2.sorted').write_bytes(changed)
                expected_error = 'Vector comparison failed'
                if failure == 'golden-mismatch':
                    # A valid stream with substituted IDs must fail against the golden.
                    changed_raw = root / 'changed.u32'
                    changed_raw.write_bytes(changed)
                    (inputs / 'gov2.sorted').write_bytes(saved_raw)
                    subprocess.run([str(binary), '-s', '-threads', '2', '-u32in', str(changed_raw),
                                    '-bvout', str(selected), '-bv-index-out', str(selected_index)],
                                   check=True, capture_output=True)
                    expected_error = 'Golden vector comparison failed'
            elif failure == 'golden-deserialize':
                golden_blob.write_bytes(bytes([255]) + saved_golden[1:])
                expected_error = 'worker failed'
            else:
                lines = saved_index.decode().splitlines()
                offsets = list(map(int, lines[1:]))
                removed = offsets[1] - 1
                damaged = saved_blob[:1] + saved_blob[offsets[1]:]
                selected.write_bytes(damaged)
                selected_index.write_text(f'BM_INV_BLOB_INDEX_V1 {len(damaged)} {len(offsets)-1}\n0\n' +
                                          '\n'.join(str(value-removed) for value in offsets[1:]) + '\n')
                expected_error = 'worker failed'
            invoke([*common, '--out', out, '--fast', '--start-level', 10, '--silent'], 1)
            phase = 'golden-verify-threaded' if failure in ('verify', 'golden-deserialize') else 'reference-verify-threaded'
            log = (out / ('gov2.sorted.level-10/' + phase + '.log')).read_text()
            assert 'Verification worker failed' in log and expected_error in log
            assert 'vector_index=0' in log and 'blob_begin=0' in log and str(golden_blob if phase == 'golden-verify-threaded' else selected) in log
            assert 'Verification complete.' not in log
            if phase == 'golden-verify-threaded':
                assert not (out / 'gov2.sorted.level-10/reference-verify-threaded.log').exists()
            result = json.loads((out / ('gov2.sorted.level-10/' + phase + '.result.json')).read_text())
            assert result['exit_code'] == 1
            records = list(csv.DictReader((out / 'summary.tsv').open(), delimiter='\t'))
            assert len(records) == 1 and records[0]['status'] == 'failed'
            assert records[0]['exit_code'] == '1'
            assert (out / 'FAILURE.txt').exists() and not (out / 'PASS.txt').exists()
        finally:
            (inputs / 'gov2.sorted').write_bytes(saved_raw)
            selected.write_bytes(saved_blob)
            selected_index.write_bytes(saved_index)
            golden_blob.write_bytes(saved_golden)
        assert snapshot(reference) == before
    # Normal round trips need only the golden stream, even for a missing selected level.
    historical = reference / 'gov2.sorted.level-10.bv'
    saved_historical = historical.with_suffix('.held')
    historical.rename(saved_historical)
    try:
        out = root / 'normal-golden-only'
        invoke([*common, '--out', out, '--start-level', 10])
        assert (out / 'PASS.txt').exists()
        assert not list(out.glob('*/reference-verify-serial.log'))
        invoke([*common, '--pedantic', '--out', root / 'pedantic-missing-history', '--start-level', 10], 1)
    finally:
        saved_historical.rename(historical)
    assert snapshot(reference) == before
    invoke([*common, '--fast', '--pedantic', '--out', root / 'incompatible-modes'], 2)
    # Provider-disabled binaries are rejected before any expensive raw verification.
    unsupported = root / 'no-sha-support'
    unsupported.write_text('#!' + sys.executable + '\nimport sys\n'
                           'if "--help" in sys.argv:\n'
                           ' print("-hash-cache-in FILE; SHA-256 cache support disabled")\n'
                           ' sys.exit(0)\n'
                           'sys.exit(77)\n')
    unsupported.chmod(0o755)
    unavailable_out = root / 'hash-unavailable'
    invoke(['--bin', unsupported, '--reference', reference, '--out', unavailable_out,
            '--fast', '--hash', '--start-level', 10, '--reserve-gib', 0], 2)
    assert not list(unavailable_out.glob('*/golden-verify-threaded.log'))
    # Cached validation is explicitly opt-in and persists next to the references.
    hash_common = [*common, '--fast', '--start-level', 7, '--hash']
    first = root / 'hash-first'
    invoke([*hash_common, '--out', first])
    cache_dir = reference / 'validation-cache'
    assert len(list(cache_dir.glob('*.sha256'))) == 2
    assert len(list(first.glob('*/golden-verify-threaded.log'))) == 2
    assert len(list(first.glob('*/hash-cache-build.log'))) == 2
    second = root / 'hash-reuse'
    invoke([*hash_common, '--out', second])
    assert not list(second.glob('*/golden-verify-threaded.log'))
    assert not list(second.glob('*/hash-cache-build.log'))
    assert len(list(second.glob('*/reference-verify-hash.log'))) == 8
    for command in second.glob('*/reference-verify-hash.command.json'):
        flags = json.loads(command.read_text())
        assert '-hash-cache-in' in flags and '-golden-bvin' not in flags and '-u32in' not in flags
    assert all(digest(reference / name) == checksum for name, checksum in before.items())
    # Independent target directory, and complete mode with the same optional path.
    third = root / 'hash-override'
    invoke([*common, '--out', third, '--start-level', 10, '--hash',
            '--hash-cache-dir', root / 'custom-hash-cache'])
    assert (third / 'PASS.txt').exists()
    # Corruption must stop verification and leave only a failure marker.
    cache = next(cache_dir.glob('*.sha256'))
    data = cache.read_bytes()
    cache.write_bytes(data[:-1])
    failed = root / 'hash-corrupt'
    invoke([*hash_common, '--out', failed], 2)
    assert (failed / 'FAILURE.txt').exists() and not (failed / 'PASS.txt').exists()
    cache.write_bytes(data)
print('Normal round trips, pedantic history, fast/hash scopes, failure propagation and reference immutability PASS')
if args.crash_test:
    print('Deliberate SIGSEGV propagation: PASS')
