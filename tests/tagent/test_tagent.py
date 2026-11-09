#!/usr/bin/env python3
"""Focused supervisor checks using temporary checkouts and synthetic test jobs."""
import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import time
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
            self.assertEqual([call.args[0] for call in departure.call_args_list],
                             ['perf', 'stress'])
        self.assertEqual(self.agent.next_train_at, float('inf'))
        runs = list(self.agent.home.glob('*abc123'))
        self.assertEqual(len(runs), 1)
        self.assertFalse((runs[0] / 'perf_debug').exists())
        self.assertEqual(len(list(runs[0].glob('stress_*/job.*/isolated-output'))), 8)
        self.assertTrue((runs[0] / 'manifest.json').exists())

    def test_stress_variants_share_history_path_across_cycles(self):
        """Every isolated stress run receives the same persistent history path."""
        self.fixture()
        calls = []
        original = self.agent.command
        def command(argv, *args, **kwargs):
            if '--order-file' in argv:
                calls.append(argv[argv.index('--order-file') + 1])
            return original(argv, *args, **kwargs)
        with patch.object(self.agent, 'command', side_effect=command):
            self.assertTrue(self.agent.cycle())
            self.assertTrue(self.agent.cycle())
        self.assertEqual(len(calls), 2 * len(self.agent.builds('stress')))
        self.assertEqual(set(calls), {str(self.agent.home / 'stress-order.tsv')})

    def test_stress_runner_reuses_explicit_history(self):
        """The real runner honors shared history in a separate script directory."""
        runner = Path(__file__).resolve().parents[1] / 'stress' / 'run_all_pl.sh'
        binary = self.root / 'worker.sh'
        binary.write_text('#!/bin/bash\nif [ "$1" = -s ]; then sleep 1; fi\nexit 0\n')
        binary.chmod(0o755)
        history = self.root / 'shared' / 'order.tsv'
        history.parent.mkdir()
        for index, expected in enumerate((['-ll', '-s'], ['-s', '-ll'])):
            isolated = self.root / ('isolated' + str(index))
            isolated.mkdir()
            copy = isolated / 'run_all_pl.sh'
            copy.write_bytes(runner.read_bytes())
            result = subprocess.run(['bash', str(copy), '--binary', str(binary), '--jobs', '1',
                                     '--order-file', str(history), '-ll', '-s'],
                                    text=True, capture_output=True, timeout=20)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            starts = [line for line in result.stdout.splitlines() if line.startswith('START  ')]
            self.assertEqual([line.split()[1] for line in starts], expected)
            self.assertFalse((isolated / 'logs/parallel/order.tsv').exists())
        self.assertIn('-s ', history.read_text())
        self.assertIn('-ll ', history.read_text())

    def test_test_failure_is_recoverable(self):
        """A nonzero perf result returns failure with logs rather than passing."""
        self.fixture(fail=True)
        self.assertFalse(self.agent.cycle())
        self.assertTrue(list(self.agent.home.glob('*abc123/perf_release/perf.log')))

    def test_build_failure_is_fatal(self):
        """Build failure propagates instead of entering test recovery."""
        self.fixture()
        (self.root / 'tests/perf/build_all.sh').write_text(
            "echo 'a.cpp:1: warning: first'\necho 'a.cpp:2: error: bad'\nprintf 'b.cpp:1: fatal error: missing'\nexit 9\n")
        stress = self.root / 'tests/stress/build_all.sh'
        stress.write_text(stress.read_text() + "echo 'stress.cpp:1: warning: second'\n")
        self.agent.args.sound_theme = 'steampunk'
        with patch.object(self.agent, 'sound') as sound:
            with self.assertRaises(ta.CommandFailure):
                self.agent.cycle()
            sound.assert_any_call('build-final-failed', workers=(2, 2, False))
        self.assertIn('BUILD FAILED | diagnostic errors: 2 | warnings: 2',
                      (self.agent.home / 'agent.log').read_text())

    def test_builds_overlap_and_combine_warnings(self):
        """Both builds reach a barrier together and contribute to the final count."""
        self.fixture()
        barrier = threading.Barrier(2, timeout=5)
        original = self.agent.command
        def command(argv, log=None, *args, **kwargs):
            """Hold both build jobs at the barrier before launching them."""
            if log and log.name in ('build-perf.log', 'build-stress.log'):
                barrier.wait()
            return original(argv, log, *args, **kwargs)
        for family in ('perf', 'stress'):
            script = self.root / 'tests' / family / 'build_all.sh'
            script.write_text(script.read_text() + "echo 'test.cpp:1: warning: synthetic'\n")
        with patch.object(self.agent, 'command', side_effect=command):
            self.assertTrue(self.agent.cycle())
        self.assertIn('BUILD COMPLETE | diagnostic errors: 0 | warnings: 2',
                      (self.agent.home / 'agent.log').read_text())

    def test_single_worker_serializes_builds(self):
        """A one-worker budget never runs both build families at once."""
        self.fixture()
        self.agent.jobs = 1
        active = [0, 0]
        mutex = threading.Lock()
        original = self.agent.command
        def command(argv, log=None, *args, **kwargs):
            """Record overlap while the synthetic builds run."""
            if not log or log.name not in ('build-perf.log', 'build-stress.log'):
                return original(argv, log, *args, **kwargs)
            with mutex:
                active[0] += 1
                active[1] = max(active)
            try:
                time.sleep(.1)
                return original(argv, log, *args, **kwargs)
            finally:
                with mutex:
                    active[0] -= 1
        with patch.object(self.agent, 'command', side_effect=command):
            self.assertTrue(self.agent.cycle())
        self.assertEqual(active[1], 1)

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

    def test_successful_fetch_reports_remote_status(self):
        """A fetch reports whether the remote tip matches the checked-out commit."""
        for tip, event, message in [('abc123', 'repo-unchanged', 'no new commits'),
                                    ('def456', None, 'new commit')]:
            with patch.object(self.agent, 'git', side_effect=['', tip]), \
                 patch.object(self.agent, 'sound') as sound:
                self.assertEqual(self.agent.fetch(), tip)
                if event:
                    sound.assert_called_once_with(event)
                else:
                    sound.assert_not_called()
            self.assertIn(message, (self.agent.home / 'agent.log').read_text())

    def test_remote_check_reports_no_changes(self):
        """A successful unchanged test-phase check prints and signals the result."""
        tip = 'a' * 40
        self.agent.commit = tip
        self.agent.next_remote_poll = 0
        with patch.object(self.agent, 'git', return_value=tip + '\trefs/heads/master'), \
             patch.object(self.agent, 'sound') as sound:
            self.agent.poll_remote()
            sound.assert_called_once_with('repo-unchanged')
        self.assertIn('REMOTE CHECK OK | no new commits', (self.agent.home / 'agent.log').read_text())

    def test_timeout_cleans_process_group(self):
        """Timeout kills and reaps the launched process and releases ownership."""
        with self.assertRaisesRegex(ta.Failure, 'timed out'):
            self.agent.command(['bash', '-c', 'sleep 60 & wait'], timeout=0.01)
        self.assertFalse(self.agent.children)

    def test_git_update_and_persistent_success(self):
        """Manual starts rebuild by default; changed-only and resume are distinct."""
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
        with patch.object(self.agent, 'validate'), patch.object(self.agent, 'sound') as sound, \
             patch.object(self.agent, 'cycle', return_value=True) as cycle:
            self.assertEqual(self.agent.run(), 0)
            sound.assert_any_call('repo-build')
            self.assertEqual(git(checkout, 'rev-parse', 'HEAD'), git(remote, 'rev-parse', 'HEAD'))
            self.assertEqual(self.agent.run(), 0)
            self.assertEqual(cycle.call_count, 2)
            self.agent.args.startup = 'changed'
            self.assertEqual(self.agent.run(), 0)
            self.assertEqual(cycle.call_count, 2)
            self.agent.args.resume_cycle = True
            self.assertEqual(self.agent.run(), 0)
            self.assertEqual(cycle.call_count, 3)
            self.agent.args.resume_cycle = False
            self.agent.args.force = True
            self.assertEqual(self.agent.run(), 0)
            self.assertEqual(cycle.call_count, 4)
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
            self.assertIn('--resume-cycle', restart.call_args.args[1])


    def test_manual_start_builds_once_then_waits(self):
        """The default startup request is consumed and does not rebuild on every poll."""
        self.agent.args.once = False
        self.agent.commit = 'a' * 40
        with patch.object(self.agent, 'validate'), \
             patch.object(self.agent, 'fetch', return_value=self.agent.commit), \
             patch.object(self.agent, 'git', return_value=self.agent.commit), \
             patch.object(self.agent, 'cycle', return_value=True) as cycle, \
             patch.object(self.agent, 'pause', side_effect=[None, StopIteration]):
            with self.assertRaises(StopIteration):
                self.agent.run()
            self.assertEqual(cycle.call_count, 1)

    def test_interrupted_cycle_restarts_even_if_tip_returns_to_previous(self):
        """A partial cycle must be retried even when the watched tip reverts."""
        import json
        self.agent.args.once = True
        commit = 'a' * 40
        (self.agent.home / 'success.json').write_text(json.dumps(
            {'commit': commit, 'pedantic': False, 'avx512': False}))
        with patch.object(self.agent, 'validate'), \
             patch.object(self.agent, 'fetch', return_value=commit), \
             patch.object(self.agent, 'git', return_value=commit), \
             patch.object(self.agent, 'cycle', side_effect=['updated', True]) as cycle:
            self.assertEqual(self.agent.run(), 0)
            self.assertEqual(cycle.call_count, 2)

    def test_failed_initial_fetch_does_not_build(self):
        """Network recovery failure never silently tests a potentially stale checkout."""
        with patch.object(self.agent, 'validate'), \
             patch.object(self.agent, 'fetch', side_effect=ta.Failure('offline')), \
             patch.object(self.agent, 'cycle') as cycle:
            with self.assertRaisesRegex(ta.Failure, 'offline'):
                self.agent.run()
            cycle.assert_not_called()

    def test_update_stops_at_perf_boundary(self):
        """A discovered update waits for all perf jobs and skips old stress builds."""
        self.fixture()
        with patch.object(self.agent, 'update_waiting', return_value=True):
            self.assertEqual(self.agent.cycle(), 'updated')
        run = next(self.agent.home.glob('*abc123'))
        self.assertTrue((run / 'perf_release/perf.log').exists())
        self.assertTrue((run / 'perf_release_neon/perf.log').exists())
        self.assertFalse(list(run.glob('stress_*/runner.log')))
        self.assertIsNone(self.agent.active_run)

    def test_update_skips_queued_perf(self):
        """A pending update lets an active perf job finish but skips queued jobs."""
        self.fixture()
        self.agent.jobs = 1
        command = self.agent.command
        def wrapped(*args, **kwargs):
            """Signal an update immediately after the first perf job finishes."""
            result = command(*args, **kwargs)
            log = args[1] if len(args) > 1 else kwargs.get('log')
            if log and log.name == 'perf.log' and result != 'skipped-update':
                self.agent.pending_tip = 'c' * 40
            return result
        with patch.object(self.agent, 'command', side_effect=wrapped):
            self.assertEqual(self.agent.cycle(), 'updated')
        run = next(self.agent.home.glob('*abc123'))
        self.assertEqual(len(list(run.glob('perf_*/isolated-output'))), 1)
        self.assertTrue(self.agent.perf_skipped)
        self.assertTrue(self.agent.args.force)

    def test_remote_poll_and_recovery(self):
        """A remote tip is remembered and transient failures retain their deadline."""
        tip = 'a' * 40
        self.agent.commit = 'b' * 40
        self.agent.next_remote_poll = 0
        self.agent.stress_stop_file = self.root / 'stop'
        now = [100.0]
        with patch.object(ta.time, 'monotonic', side_effect=lambda: now[0]), \
             patch.object(self.agent, 'git', side_effect=[ta.CommandFailure('offline'),
                 tip + '\trefs/heads/master']) as git:
            self.agent.poll_remote()
            self.assertEqual(self.agent.poll_failure_since, 100)
            now[0] = 159
            self.agent.poll_remote()
            self.assertEqual(git.call_count, 1)
            now[0] = 160
            self.agent.poll_remote()
            self.assertIsNone(self.agent.poll_failure_since)
            self.assertEqual(self.agent.pending_tip, tip)
            self.assertTrue(self.agent.stress_stop_file.exists())
            self.assertTrue(self.agent.update_waiting())

    def test_cleanup_preserves_unowned_and_active(self):
        """Retention removes only expired marked runs, not active or unowned paths."""
        import datetime
        import json
        old = (datetime.datetime.now() - datetime.timedelta(days=8)).strftime('%Y%m%d-%H%M%S-%f')
        fresh = datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f')
        paths = [self.agent.home / (old + '-' + suffix) for suffix in ('aaaaaa', 'bbbbbb', 'cccccc')]
        for path in paths:
            path.mkdir()
        for path in paths[:2]:
            (path / 'manifest.json').write_text(json.dumps({'tagent_run': 1}))
        recent = self.agent.home / (fresh + '-dddddd')
        recent.mkdir()
        (recent / 'manifest.json').write_text(json.dumps({'tagent_run': 1}))
        self.agent.active_run = paths[1]
        self.agent.cleanup_logs()
        self.assertFalse(paths[0].exists())
        self.assertTrue(paths[1].exists())
        self.assertTrue(paths[2].exists())
        self.assertTrue(recent.exists())

    def test_running_build_reports_warning_count(self):
        """A long real subprocess reports accumulated warnings while still running."""
        self.agent.args.sound_theme = 'steampunk'
        self.agent.next_build_report = 0
        log = self.root / 'build-perf.log'
        with patch.object(self.agent, 'sound') as sound:
            self.agent.command(['bash', '-c',
                "printf 'a.cpp:1:2: warning: first\\n'; sleep 11; printf 'b.cpp:3: warning: second\\n'"], log)
            sound.assert_called_once_with('build-warning', workers=(1, 0, 0, 0, 0))
        self.assertEqual(self.agent.build_warnings.observe(log), 2)

    def test_perf_progress_ignores_blank_and_rotates(self):
        """Only completed timing lines count, and reports rotate across active variants."""
        for name in ('perf_release', 'perf_release_avx2'):
            directory = self.root / name
            directory.mkdir()
            (directory / 'perf.log').write_text(
                'BitMagic Library\n\nSIMD code = 8\n\nfirst; 4.924 sec\n\n'
                'InterpolativeCodingTest() OK\nsecond ; 1 min\nunfinished; 2 sec')
        self.agent.perf_logs = {1: self.root / 'perf_release' / 'perf.log',
                                2: self.root / 'perf_release_avx2' / 'perf.log'}
        self.agent.next_perf_report = 0
        with patch.object(ta.time, 'monotonic', return_value=100), \
             patch.object(self.agent, 'sound') as sound:
            self.agent.report_perf_progress()
            sound.assert_called_once_with('perf-progress', workers=('perf_release', 2))
            self.agent.report_perf_progress()
            self.assertEqual(sound.call_count, 1)
            self.agent.next_perf_report = 0
            self.agent.report_perf_progress()
            sound.assert_called_with('perf-progress', workers=('perf_release_avx2', 2))
        output = (self.agent.home / 'agent.log').read_text()
        self.assertIn('perf_release | 2 completed | approximately 1% of test elements | last: second ; 1 min', output)
        self.assertNotIn('unfinished; 2 sec', output)

    def test_build_progress_counts_both_logs_once(self):
        """Concurrent heartbeat reports observe both logs without double counting."""
        logs = [self.root / 'build-perf.log', self.root / 'build-stress.log']
        for path in logs:
            path.write_text('x.cpp:1: warning: synthetic\n')
        self.agent.build_logs = logs
        self.agent.args.sound_theme = 'steampunk'
        self.agent.next_build_report = 0
        with patch.object(self.agent, 'sound') as sound:
            self.agent.build_progress(logs[0])
            self.agent.build_progress(logs[1])
            sound.assert_called_once_with('build-warning', workers=(2, 0, 0, 0, 0))
        self.assertEqual(self.agent.build_warnings.count, 2)

    def test_build_progress_tracks_targets_and_increases_interval(self):
        """Snapshots count completed files and each report adds ten seconds."""
        log = self.root / 'build-perf.log'
        log.write_text('')
        perf = [self.root / 'perf1', self.root / 'perf2']
        stress = [self.root / 'stress1', self.root / 'stress2', self.root / 'stress3']
        for path in (perf[0], stress[0], stress[1]):
            path.write_text('ready')
            path.chmod(0o755)
        self.agent.build_targets = {'perf': perf, 'stress': stress}
        self.agent.args.sound_theme = 'steampunk'
        self.agent.next_build_report = 60
        with patch.object(self.agent, 'sound') as sound, \
             patch.object(ta.time, 'monotonic') as clock:
            clock.return_value = 60
            self.agent.build_progress(log)
            sound.assert_called_with('build-clear', workers=(0, 1, 1, 2, 1))
            self.assertEqual(self.agent.next_build_report, 130)
            clock.return_value = 129
            self.agent.build_progress(log)
            self.assertEqual(sound.call_count, 1)
            clock.return_value = 130
            self.agent.build_progress(log)
            self.assertEqual(self.agent.next_build_report, 210)
        self.assertIn('perf: 1 ready, 1 pending | stress: 2 ready, 1 pending',
                      (self.agent.home / 'agent.log').read_text())

    def test_running_job_polls_real_remote(self):
        """A commit arriving during a long job is detected without changing its checkout."""
        import threading
        import time
        remote = self.root / 'origin'
        checkout = self.root / 'checkout'
        remote.mkdir()
        def git(directory, *args):
            """Execute Git in a temporary repository."""
            return subprocess.check_output(['git', '-C', str(directory), *args],
                                           stderr=subprocess.DEVNULL, text=True).strip()
        git(remote, 'init', '-b', 'master')
        git(remote, 'config', 'user.name', 'Fixture')
        git(remote, 'config', 'user.email', 'fixture@example.invalid')
        (remote / 'source').write_text('one')
        git(remote, 'add', '.')
        git(remote, 'commit', '-m', 'first')
        subprocess.check_call(['git', 'clone', '-q', str(remote), str(checkout)])
        self.agent.repo = checkout
        self.agent.commit = git(checkout, 'rev-parse', 'HEAD')
        self.agent.args.poll_minutes = .001
        self.agent.next_remote_poll = 0
        self.agent.stress_stop_file = self.root / 'stop'
        ready = self.root / 'ready'
        def publish():
            """Advance the remote after the running job has started."""
            deadline = time.monotonic() + 5
            while not ready.exists() and time.monotonic() < deadline:
                time.sleep(.01)
            if ready.exists():
                (remote / 'source').write_text('two')
                git(remote, 'commit', '-am', 'second')
        publisher = threading.Thread(target=publish)
        publisher.start()
        try:
            self.agent.command(['bash', '-c', 'echo START  -ll; touch "$1"; sleep 11; echo PASS   -ll',
                                'job', str(ready)], self.root / 'runner.log')
            publisher.join(timeout=5)
            self.assertEqual(self.agent.pending_tip, git(remote, 'rev-parse', 'HEAD'))
            self.assertEqual(git(checkout, 'rev-parse', 'HEAD'), self.agent.commit)
            self.assertTrue(self.agent.stress_stop_file.exists())
        finally:
            publisher.join(timeout=5)

    def test_runner_drains_active_group_on_update(self):
        """A stop request prevents queued groups starting and lets the active one finish."""
        import shutil
        import time
        runner = self.root / 'run_all_pl.sh'
        shutil.copy2(Path(__file__).resolve().parents[1] / 'stress/run_all_pl.sh', runner)
        worker = self.root / 'worker.sh'
        worker.write_text('#!/bin/bash\necho "$1" >> started\nsleep 1\necho "$1" >> finished\n')
        worker.chmod(0o755)
        stop = self.root / 'stop'
        proc = subprocess.Popen(['bash', str(runner), '--binary', './worker.sh',
                                 '--jobs', '1', '--stop-file', str(stop), '-ll', '-s'],
                                cwd=self.root, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        try:
            deadline = time.monotonic() + 5
            while not (self.root / 'started').exists() and time.monotonic() < deadline:
                time.sleep(.01)
            self.assertTrue((self.root / 'started').exists())
            stop.touch()
            _, stderr = proc.communicate(timeout=10)
            self.assertEqual(proc.returncode, 75, stderr.decode())
            started = (self.root / 'started').read_text().splitlines()
            self.assertEqual(len(started), 1)
            self.assertEqual((self.root / 'finished').read_text().splitlines(), started)
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()

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
