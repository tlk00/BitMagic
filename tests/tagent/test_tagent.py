#!/usr/bin/env python3
"""Focused supervisor checks using temporary checkouts and synthetic test jobs."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SPEC = importlib.util.spec_from_file_location('tagent', Path(__file__).with_name('tagent.py'))
ta = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ta)


class SupervisorTests(unittest.TestCase):
    """Exercise orchestration without compiling or running BitMagic workloads."""

    def setUp(self):
        """Create isolated agent state and a synthetic source tree."""
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        with patch.object(sys, 'argv', ['tagent', '--no-sound', '--min-free-gb', '0.001']):
            args = ta.arguments()
        args.state_dir = self.root / 'state'
        self.agent = ta.Agent(args)
        self.agent.repo = self.root
        self.agent.arm = True
        self.agent.commit = 'abc123'
        self.agent.jobs = 2

    def tearDown(self):
        """Stop synthetic processes and remove temporary evidence."""
        self.agent.stop_children()
        self.temp.cleanup()

    def fixture(self, fail=False):
        """Create build scripts that produce tiny executable jobs; fail selects perf failure."""
        for family in ('perf', 'stress'):
            directory = self.root / 'tests' / family
            directory.mkdir(parents=True)
            lines = ['#!/bin/bash', 'set -eu', 'cd "$(dirname "$0")"']
            for name in self.agent.builds(family):
                status = 7 if fail and name == 'perf_release' else 0
                lines += ["cat > " + name + " <<'JOB'", '#!/bin/bash',
                          'touch isolated-output', 'exit ' + str(status), 'JOB', 'chmod +x ' + name]
            (directory / 'build_all.sh').write_text('\n'.join(lines) + '\n')
            if family == 'stress':
                (directory / 'run_all_pl.sh').write_text(
                    '#!/bin/bash\nset -eu\ncd "$(dirname "$0")"\n"$2" -ll\n"$2" -s\n')

    def test_success_and_isolation(self):
        """All expected builds run with perf release-only and isolated stress files."""
        self.fixture()
        with patch.object(self.agent, 'start_test_sound', wraps=self.agent.start_test_sound) as departure:
            self.assertTrue(self.agent.cycle())
            departure.assert_called_once()
        self.assertEqual(self.agent.next_train_at, float('inf'))
        runs = list(self.agent.home.glob('*abc123'))
        self.assertEqual(len(runs), 1)
        self.assertFalse((runs[0] / 'perf_debug').exists())
        self.assertEqual(len(list(runs[0].glob('stress_*/job.*/isolated-output'))), 8)
        self.assertTrue((runs[0] / 'manifest.json').exists())

    def test_test_failure_is_recoverable(self):
        """A nonzero perf result returns failure with logs rather than passing."""
        self.fixture(fail=True)
        self.assertFalse(self.agent.cycle())
        self.assertTrue(list(self.agent.home.glob('*abc123/perf_release/perf.log')))

    def test_build_failure_is_fatal(self):
        """Build failure propagates immediately instead of entering test recovery."""
        self.fixture()
        (self.root / 'tests/perf/build_all.sh').write_text('exit 9\n')
        with self.assertRaises(ta.CommandFailure):
            self.agent.cycle()

    def test_recovery_deadline_without_new_commit(self):
        """Successful fetches with an unchanged tip do not reset the failure window."""
        self.agent.args.retry_minutes = 1
        clock = [0]
        delays = []
        def pause(seconds, label):
            """Advance fake monotonic time by the requested wait."""
            delays.append(seconds)
            clock[0] += seconds
        with patch.object(ta.time, 'monotonic', side_effect=lambda: clock[0]), \
             patch.object(self.agent, 'pause', side_effect=pause), \
             patch.object(self.agent, 'git', return_value='bad'):
            with self.assertRaisesRegex(ta.Failure, 'Recovery window exhausted'):
                self.agent.fetch('bad')
        self.assertEqual(delays, [60])

    def test_network_recovers(self):
        """A transient fetch failure backs off and accepts a later tip."""
        with patch.object(self.agent, 'git', side_effect=[ta.CommandFailure('offline'), '', 'new']), \
             patch.object(self.agent, 'pause') as pause:
            self.assertEqual(self.agent.fetch(), 'new')
            self.assertEqual(pause.call_count, 1)

    def test_timeout_cleans_process_group(self):
        """Timeout kills and reaps the launched process and releases ownership."""
        with self.assertRaisesRegex(ta.Failure, 'timed out'):
            self.agent.command(['bash', '-c', 'sleep 60 & wait'], timeout=0.01)
        self.assertFalse(self.agent.children)

    def test_git_update_and_persistent_success(self):
        """A real local remote advances safely; a completed tip is skipped on restart."""
        remote = self.root / 'remote'
        checkout = self.root / 'checkout'
        def git(directory, *args):
            """Run fixture Git without using the developer checkout."""
            return subprocess.check_output(['git', '-C', str(directory), *args],
                                           stderr=subprocess.DEVNULL, text=True).strip()
        remote.mkdir()
        git(remote, 'init', '-b', 'master')
        git(remote, 'config', 'user.name', 'Fixture')
        git(remote, 'config', 'user.email', 'fixture@example.invalid')
        (remote / 'source').write_text('one')
        git(remote, 'add', '.')
        git(remote, 'commit', '-m', 'first')
        subprocess.check_call(['git', 'clone', '-q', str(remote), str(checkout)])
        (remote / 'source').write_text('two')
        git(remote, 'commit', '-am', 'second')
        self.agent.repo = checkout
        self.agent.args.once = True
        with patch.object(self.agent, 'validate'), patch.object(self.agent, 'cycle', return_value=True) as cycle:
            self.assertEqual(self.agent.run(), 0)
            self.assertEqual(git(checkout, 'rev-parse', 'HEAD'), git(remote, 'rev-parse', 'HEAD'))
            self.assertEqual(self.agent.run(), 0)
            self.assertEqual(cycle.call_count, 1)
        (remote / 'tests/tagent').mkdir(parents=True)
        (remote / 'tests/tagent/version').write_text('new agent')
        git(remote, 'add', '.')
        git(remote, 'commit', '-m', 'agent update')
        self.agent.lock_file = (self.root / 'lock').open('a')
        with patch.object(self.agent, 'validate'), \
             patch.object(ta.os, 'execv', side_effect=RuntimeError('restart sentinel')) as restart:
            with self.assertRaisesRegex(RuntimeError, 'restart sentinel'):
                self.agent.run()
            self.assertEqual(restart.call_args.args[0], '/bin/bash')
            self.assertIn('run_tagent.sh', restart.call_args.args[1][1])


    def test_existing_runner(self):
        """The real stress runner accepts the wrapper interface and records group logs."""
        import shutil
        runner = self.root / 'run_all_pl.sh'
        shutil.copy2(Path(__file__).resolve().parents[1] / 'stress/run_all_pl.sh', runner)
        worker = self.root / 'worker.sh'
        worker.write_text('#!/bin/bash\nexit 0\n')
        worker.chmod(0o755)
        self.agent.command(['bash', str(runner), '--binary', './worker.sh',
                            '--jobs', '2', '--fail-fast', '-ll', '-s'],
                           self.root / 'runner.log')
        self.assertEqual(len(list(self.root.glob('logs/parallel/worker.sh/*.log'))), 2)

    def test_interrupt_stops_children(self):
        """SIGINT exits 130 and prevents a descendant from completing its work."""
        import signal
        import time
        marker = self.root / 'survived'
        ready = self.root / 'ready'
        code = (
            'import sys; sys.path.insert(0, ' + repr(str(Path(__file__).parent.resolve())) + '); '
            'import tagent; '
            'sys.argv=["tagent", "--no-sound", "--state-dir", ' + repr(str(self.root / 'interrupt-state')) + ']; '
            'tagent.Agent.run=lambda self: self.command(["bash", "-c", '
            + repr('touch ' + str(ready) + '; (sleep 2; touch ' + str(marker) + ') & wait') + ']); '
            'sys.exit(tagent.main())')
        proc = subprocess.Popen([sys.executable, '-c', code], stdout=subprocess.DEVNULL)
        try:
            deadline = time.monotonic() + 5
            while not ready.exists() and time.monotonic() < deadline:
                time.sleep(0.01)
            self.assertTrue(ready.exists())
            proc.send_signal(signal.SIGINT)
            self.assertEqual(proc.wait(timeout=10), 130)
            self.assertFalse(marker.exists())
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()

    def test_pedantic_matrix(self):
        """x86 pedantic perf includes debug/release for every SIMD suffix."""
        self.agent.arm = False
        self.agent.args.pedantic = True
        for suffix in ('sse2', 'sse42', 'avx2', '64'):
            for mode in ('debug', 'release'):
                self.assertIn('perf_' + mode + '_' + suffix, self.agent.builds('perf'))


if __name__ == '__main__':
    unittest.main()
