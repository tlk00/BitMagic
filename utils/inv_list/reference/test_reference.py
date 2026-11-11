#!/usr/bin/env python3
"""Exercise the real reference matrix and failure/resume rules on tiny collections."""
from array import array
import csv
import json
import hashlib
from pathlib import Path
import random
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
BINARY = HERE.parent / 'bminv_reference_neon'
if len(sys.argv) > 1:
    BINARY = Path(sys.argv[1]).resolve()


def run(arguments, success=True):
    result = subprocess.run([sys.executable, str(HERE / 'gov2_reference.py'), *map(str, arguments)],
                            capture_output=True, text=True)
    assert (result.returncode == 0) == success, (result.stdout[-5000:], result.stderr)
    return result


with tempfile.TemporaryDirectory(prefix='bm-reference-test-') as temporary:
    root = Path(temporary)
    inputs = root / 'inputs'
    inputs.mkdir()
    rng = random.Random(42)
    rows = [sorted(rng.sample(range(65536), 1800)), list(range(100, 3000)),
            [k * 1009 + 3 for k in range(120)], [1, 65535]]
    for dataset in ('sorted', 'unsorted'):
        with (inputs / ('gov2.' + dataset)).open('wb') as stream:
            for row in (rows if dataset == 'sorted' else rows[::-1]):
                array('I', [len(row), *row]).tofile(stream)
    output = root / 'reference'
    common = ['--allow-local-output', '--reserve-gib', '0', '--input-dir', inputs,
              '--out', output, '--bin', BINARY, '-j', '2']
    run(common)
    with (output / 'references.tsv').open() as stream:
        records = list(csv.DictReader(stream, delimiter='\t'))
    assert len(records) == 20
    assert {row['status'] for row in records} == {'validated'}
    assert len({row['integer_count'] for row in records}) == 1
    manifests = list((output / 'reports').glob('*/metadata.json'))
    assert len(manifests) == 20
    before = {path: path.stat().st_mtime_ns for path in manifests}
    again = run(common)
    assert again.stdout.count('Integrity checked; skipping') == 20
    assert before == {path: path.stat().st_mtime_ns for path in manifests}
    # Corruption must be detected before reusing or generating entries.
    target = output / records[0]['stream_file']
    content = target.read_bytes()
    target.write_bytes(content[:-1] + bytes([content[-1] ^ 1]))
    assert 'integrity failure' in run(common, False).stderr
    target.write_bytes(content)
    index = target.with_suffix('.idx')
    original_index = index.read_bytes()
    index.write_bytes(original_index + b'1\n')
    assert 'integrity failure' in run(common, False).stderr
    index.write_bytes(original_index)
    # A missing commit marker must never cause an existing stream to be overwritten.
    manifest = output / 'reports' / target.stem / 'metadata.json'
    metadata = manifest.read_bytes()
    manifest.unlink()
    result = run(common, False)
    assert 'Refusing to overwrite' in result.stderr or 'Uncommitted' in result.stderr or 'Incomplete' in result.stderr
    assert target.read_bytes() == content
    manifest.write_bytes(metadata)
    # Space preflight fails before starting any phase.
    result = run([*common, '--reserve-gib', '100000000'], False)
    assert 'Insufficient free space' in result.stderr
    # Lock acquisition is exclusive.
    (output / '.reference-lock').mkdir()
    result = run(common, False)
    assert '.reference-lock' in result.stderr
    (output / '.reference-lock').rmdir()
    # A failed serial comparison must never publish a reference.
    wrapper = root / 'injected-reader-failure'
    real_binary = str(BINARY.resolve())
    wrapper.write_text('#!' + sys.executable + '\nimport os, sys\n'
                       'if "-verify" in sys.argv and "-threads" not in sys.argv:\n'
                       '    print("Injected serial verification failure"); sys.exit(17)\n'
                       'os.execv(' + repr(real_binary) + ', [' + repr(real_binary) + '] + sys.argv[1:])\n')
    wrapper.chmod(0o755)
    build = json.loads(BINARY.with_suffix('.build.json').read_text())
    build['binary_sha256'] = hashlib.sha256(wrapper.read_bytes()).hexdigest()
    build_log = wrapper.with_suffix('.build.log')
    build_log.write_text('Fixture-only injected serial-reader failure\n')
    build['build_log_sha256'] = hashlib.sha256(build_log.read_bytes()).hexdigest()
    wrapper.with_suffix('.build.json').write_text(json.dumps(build))
    failed = root / 'failed-reference'
    result = run(['--allow-local-output', '--reserve-gib', '0', '--input-dir', inputs,
                  '--out', failed, '--bin', wrapper, '-j', '2', '--levels', '6', '--datasets', 'sorted'], False)
    assert 'verify-serial failed' in result.stderr
    assert not list(failed.glob('*.bv'))
    assert list((failed / '.partial').glob('*/FAILURE.txt'))
    # Different raw content cannot be silently substituted under the same dataset name.
    source = inputs / 'gov2.sorted'
    source.write_bytes(source.read_bytes() + b'\0')
    assert 'Input differs' in run(common, False).stderr
print('20 real codec cases: both readers, decode, accounting, resume, immutability, corruption, space and locking PASS')
