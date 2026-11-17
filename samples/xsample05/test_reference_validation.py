"""Synthetic reader checks: portability, retained output and fail-fast process cleanup."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

SCRIPT = Path(__file__).with_name('reference_matrix.py')


class ValidationTests(unittest.TestCase):
    """Run the real validation CLI against tiny immutable reference fixtures."""

    def test_validation_and_failures(self):
        """A moved bundle validates, while reader errors and signals stop siblings."""
        for fault in ('ok', 'exit 7', 'kill -TERM $$', 'interrupt'):
            with self.subTest(fault=fault), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                reference = root / 'reference'
                (reference / 'inputs').mkdir(parents=True)
                sha = hashlib.sha256(b'x').hexdigest()
                inputs = {k: dict(path='/unavailable/'+k, sha256=sha) for k in ('plain', 'remap')}
                for k in inputs:
                    (reference / 'inputs' / (k+'.sv')).write_bytes(b'x')
                (reference / 'matrix.json').write_text(json.dumps(dict(schema=1, inputs=inputs)))
                for remap in (0, 1):
                    for xor in (0, 1):
                        name = f'L00-remap{remap}-xor{xor}'
                        folder = reference / 'cases' / name
                        folder.mkdir(parents=True)
                        files = ('reference.sv', 'opcodes.txt', 'create.log', 'command.json', 'validation.log')
                        for f in files:
                            (folder / f).write_bytes(b'x')
                        metadata = dict(case=name, level=0, remap=bool(remap), xor=bool(xor), bytes=1,
                                        input_sha256=sha, sha256={f:sha for f in files})
                        (folder / 'metadata.json').write_text(json.dumps(metadata))
                binary = root / 'reader'
                body = 'echo verified\n' if fault == 'ok' else f'case "$2" in *remap0-xor0*) sleep 1; {fault};; *) sleep 20; touch "{root}/survived";; esac\n'
                if fault == 'interrupt':
                    body = 'sleep 20\ntouch "' + str(root/'survived') + '"\n'
                binary.write_text('#!/bin/sh\n'+body)
                binary.chmod(0o755)
                before = {str(f.relative_to(reference)): f.read_bytes() for f in reference.rglob('*') if f.is_file()}
                start = time.monotonic()
                command = [sys.executable, str(SCRIPT), 'validate', '--reference-dir', str(reference),
                           '--bin', str(binary), '--levels', '0', '-j', '2', '--out', str(root/'results')]
                if fault == 'interrupt':
                    child = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    time.sleep(.5)
                    child.terminate()
                    output, errors = child.communicate(timeout=10)
                    result = subprocess.CompletedProcess(command, child.returncode, output, errors)
                else:
                    result = subprocess.run(command, capture_output=True, text=True, timeout=15)
                self.assertEqual(result.returncode == 0, fault == 'ok', result.stdout+result.stderr)
                self.assertLess(time.monotonic()-start, 10)
                self.assertFalse((root/'survived').exists())
                self.assertTrue((root/'results/progress.json').is_file())
                if fault != 'ok':
                    self.assertIn('exit 7' if fault == 'exit 7' else 'SIGTERM', result.stderr)
                else:
                    self.assertEqual(len(list((root/'results').glob('*/validation.log'))), 4)
                self.assertEqual(before, {str(f.relative_to(reference)): f.read_bytes() for f in reference.rglob('*') if f.is_file()})


if __name__ == '__main__':
    unittest.main()
