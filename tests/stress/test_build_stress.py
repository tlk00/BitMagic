"""Isolated builder fault and flag checks using a synthetic Make executable."""
import importlib.util
import json
import os
import re
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location('build_stress', Path(__file__).with_name('build_stress.py'))
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)


class BuildTests(unittest.TestCase):
    """Exercise concurrent output isolation, diagnostics and descendant cleanup."""

    def run_fixture(self, failure, interrupt=False):
        """Run two fake variants; optionally fail one while its sibling is active."""
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            tools = root / 'tools'
            tools.mkdir()
            make = tools / 'make'
            make.write_text('#!' + sys.executable + '\n' + """
import json, os, signal, subprocess, sys, time
from pathlib import Path
values = dict(arg.split('=', 1) for arg in sys.argv[1:] if '=' in arg)
base = Path(values['BUILD_ROOT'])
binary = Path(values['DESTINATION'])
(base / 'arguments.json').write_text(json.dumps(sys.argv[1:]))
(base / 'started').touch()
deadline = time.monotonic() + 5
while len(list(base.parent.glob('*/started'))) < 2:
    if time.monotonic() > deadline: sys.exit(91)
    time.sleep(.01)
print('source.cpp:1: warning: injected warning', flush=True)
if os.environ['INJECT_FAILURE'] == '1':
    if binary.name == 'stress_release' and os.environ.get('INTERRUPT') != '1':
        # Let the sibling record its descendant before failing.
        deadline = time.monotonic() + 5
        while not (base.parent / 'descendant.pid').exists():
            if time.monotonic() > deadline: sys.exit(92)
            time.sleep(.01)
        print('source.cpp:2: error: injected failure', end='', flush=True)
        sys.exit(17)
    child = subprocess.Popen([sys.executable, '-c',
        'import signal,time; signal.signal(signal.SIGTERM,signal.SIG_IGN); time.sleep(60)'])
    (base / 'descendant.pid').write_text(str(child.pid))
    (base.parent / 'descendant.pid').write_text(str(child.pid))
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    time.sleep(60)
else:
    binary.write_text('binary for ' + binary.name)
    binary.chmod(0o755)
""")
            make.chmod(0o755)
            env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ['PATH'],
                       INJECT_FAILURE='1' if failure else '0', INTERRUPT='1' if interrupt else '0')
            code = (
                'import sys; sys.path.insert(0, ' + repr(str(Path(__file__).parent.resolve())) + '); '
                'import build_stress as b; from pathlib import Path; from types import SimpleNamespace; '
                'import signal; signal.signal(signal.SIGTERM,b.interrupted); '
                'sys.exit(b.build(SimpleNamespace(jobs=2,build_dir=Path(sys.argv[1])),'
                'Path(sys.argv[2]),b.variants(True,False)[:2]))')
            command = [sys.executable, '-c', code, str(root / 'outputs'), str(root)]
            if interrupt:
                import time
                proc = subprocess.Popen(command, env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                try:
                    deadline = time.monotonic() + 6
                    while len(list((root / 'outputs').glob('*/descendant.pid'))) < 2 and time.monotonic() < deadline:
                        time.sleep(.01)
                    proc.terminate()
                    stdout, stderr = proc.communicate(timeout=6)
                    result = subprocess.CompletedProcess(command, proc.returncode, stdout, stderr)
                finally:
                    if proc.poll() is None:
                        proc.kill()
                        proc.wait()
            else:
                result = subprocess.run(command, env=env, capture_output=True, text=True, timeout=12)
            sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tagent'))
            from tagent import WarningCounter
            aggregate = root / 'aggregate.log'
            aggregate.write_text(result.stdout)
            counter = WarningCounter()
            counter.observe(aggregate, final=True)
            self.assertEqual(counter.count, 2, result.stdout + result.stderr)
            if failure:
                self.assertEqual(result.returncode, 143 if interrupt else 1, result.stderr)
                self.assertEqual(counter.errors, 0 if interrupt else 1)
                if not interrupt:
                    report = json.loads((root / 'outputs/failure.json').read_text())
                    self.assertEqual((report['variant'], report['exit_code']), ('stress_release', 17))
                self.assertFalse((root / 'stress_release').exists())
                # A killed orphan may remain a zombie until its OS reaper runs.
                for pid_file in (root / 'outputs').glob('*/descendant.pid'):
                    pid = int(pid_file.read_text())
                    state = subprocess.run(['ps', '-o', 'stat=', '-p', str(pid)],
                                           capture_output=True, text=True).stdout.strip()
                    self.assertTrue(not state or state.startswith('Z'), state)
            else:
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(counter.errors, 0)
                for name in ('stress_debug', 'stress_release'):
                    self.assertEqual((root / name).read_text(), 'binary for ' + name)
                    arguments = json.loads((root / 'outputs' / name / 'arguments.json').read_text())
                    self.assertIn('-j1', arguments)
                    self.assertEqual('DEBUG=YES' in arguments, name == 'stress_debug')

    def test_parallel_success_and_warning_accounting(self):
        """Both variants overlap, publish separately and contribute warnings once."""
        self.run_fixture(False)

    def test_failure_stops_descendants_and_preserves_error(self):
        """Exit 17 survives cleanup of a TERM-ignoring compiler descendant."""
        self.run_fixture(True)

    def test_termination_stops_compiler_groups(self):
        """A supervisor SIGTERM also stops isolated compiler groups."""
        self.run_fixture(True, interrupt=True)

    def test_isolated_compiler_options_match_legacy(self):
        """Dry Make expansions retain architecture, SIMD, optimization and debug flags."""
        repo = Path(__file__).resolve().parents[2]
        for target, kernels in (('aarch64-apple-darwin', ['', 'NEON']),
                                ('amd64-linux-gnu', ['', 'SSE2', 'SSE42', 'AVX2', 'AVX512'])):
            for kernel in kernels:
                for debug in (False, True):
                    with self.subTest(target=target, kernel=kernel, debug=debug):
                        command = ['make', '-C', str(repo / 'tests/stress'), '-n', '-B', 'OSNAME=' + target]
                        if kernel:
                            command.append('BMOPTFLAGS=-DBM' + kernel + 'OPT')
                        if debug:
                            command.append('DEBUG=YES')
                        lines = []
                        for extra in ([], ['BUILD_ROOT=/tmp/bm-stress-flag-check', 'DESTINATION=/tmp/bm-stress-flag-check/app']):
                            result = subprocess.run(command + extra + ['all'], capture_output=True, text=True)
                            self.assertEqual(result.returncode, 0, result.stderr)
                            line = next(line for line in result.stdout.splitlines() if ' -c ' in line and 't.cpp' in line)
                            lines.append(re.sub(r'-o \S+', '-o OBJECT', line))
                        self.assertEqual(lines[0], lines[1])

    def test_platform_matrix_preserves_defines(self):
        """The established SIMD/debug variable combinations remain explicit."""
        matrix = dict(builder.variants(False, True))
        self.assertEqual(len(matrix), 11)
        for suffix, define in (('sse2', 'BMSSE2OPT'), ('sse42', 'BMSSE42OPT'),
                               ('avx2', 'BMAVX2OPT'), ('avx512', 'BMAVX512OPT')):
            self.assertEqual(matrix['stress_debug_' + suffix], ['DEBUG=YES', 'BMOPTFLAGS=-D' + define])
            self.assertEqual(matrix['stress_release_' + suffix], ['BMOPTFLAGS=-D' + define])
        self.assertEqual(dict(builder.variants(True, False))['stress_release_neon'],
                         ['BMOPTFLAGS=-DBMNEONOPT'])


if __name__ == '__main__':
    unittest.main()
