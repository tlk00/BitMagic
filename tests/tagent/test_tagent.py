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
        args.suppress.add('inv_list')
        self.agent = ta.Agent(args)
        self.agent.repo = self.root
        self.agent.arm = True
        self.agent.commit = 'abc123'
        self.agent.jobs = 2

    def tearDown(self):
        """Stop synthetic processes and remove temporary evidence."""
        self.agent.stop_children()
        self.temp.cleanup()

    def test_jemalloc_detection_and_preload(self):
        """Select a loadable library and preserve caller preloads without mutation."""
        library = self.root / 'libjemalloc.so.2'
        library.touch()
        listing = subprocess.CompletedProcess([], 0, 'libjemalloc.so.2 (libc6) => ' + str(library), '')
        probe = subprocess.CompletedProcess([], 0, '', '')
        with patch.object(ta.platform, 'system', return_value='Linux'), \
                patch.object(ta.shutil, 'which', return_value='/sbin/ldconfig'), \
                patch.object(ta.subprocess, 'run', side_effect=[listing, probe]) as run:
            self.agent.detect_jemalloc()
        self.assertEqual(self.agent.jemalloc, str(library))
        self.assertTrue(run.call_args.kwargs['env']['LD_PRELOAD'].startswith(str(library)))
        original = {'LD_PRELOAD': '/existing.so'}
        self.assertEqual(self.agent.test_environment(original)['LD_PRELOAD'], str(library) + ':/existing.so')
        self.assertEqual(original, {'LD_PRELOAD': '/existing.so'})

    def test_jemalloc_loader_rejection_falls_back(self):
        """Loader warnings count as rejection even if the probe exits successfully."""
        library = self.root / 'libjemalloc.so.2'
        library.touch()
        with patch.object(ta.platform, 'system', return_value='Linux'), \
                patch.object(ta.shutil, 'which', return_value='/sbin/ldconfig'), \
                patch.object(ta.subprocess, 'run', side_effect=[
                    subprocess.CompletedProcess([], 0, 'libjemalloc.so.2 => ' + str(library), ''),
                    subprocess.CompletedProcess([], 0, '', 'wrong ELF class')]):
            self.agent.detect_jemalloc()
        self.assertIsNone(self.agent.jemalloc)

    def test_dispatch_reply_waits_for_busy_audio(self):
        """Dispatch joins the queue and runs after the active announcement."""
        entered = threading.Event()
        release = threading.Event()
        order = []
        def active():
            entered.set()
            release.wait(3)
            order.append('prior')
        previous = threading.Thread(target=active)
        previous.start()
        self.assertTrue(entered.wait(2))
        self.agent.audio_thread = previous
        with patch.object(ta.random, 'uniform', return_value=0), \
             patch.object(self.agent, 'play_sound', side_effect=lambda event, _: order.append(event) or True):
            self.agent.queue_joker_reply(self.agent.audio_epoch, self.agent.audio_phase_generation)
            reply = self.agent.audio_thread
            self.assertIsNot(reply, previous)
            self.assertEqual(order, [])
            release.set()
            reply.join(3)
        self.assertEqual(order, ['prior', 'joker-reply'])

    def test_audio_player_has_buffering_grace_period(self):
        """Short clips allow WSL audio startup and buffering time."""
        self.agent.audio_dir.mkdir()
        self.agent.args.tts = 'off'
        with patch.object(ta.soundtrack, 'render_dispatch', return_value=5), \
             patch.object(ta.shutil, 'which', side_effect=lambda x: '/fake/paplay' if x == 'paplay' else None), \
             patch.object(self.agent, 'command', return_value='') as command:
            self.assertTrue(self.agent.play_sound('dispatch', 1))
        self.assertEqual(command.call_args.kwargs['timeout'], 30)

    def test_audio_player_failure_is_logged_and_retried(self):
        """A transient player error does not disable later announcements."""
        self.agent.audio_dir.mkdir()
        self.agent.args.tts = 'off'
        with patch.object(ta.soundtrack, 'render_dispatch', return_value=1), \
             patch.object(ta.shutil, 'which', side_effect=lambda x: '/fake/paplay' if x == 'paplay' else None), \
             patch.object(self.agent, 'command', side_effect=[ta.CommandFailure('device busy'), '']), \
             patch('builtins.open', side_effect=self._deny_console):
            self.agent.play_sound('dispatch', 1)
            self.assertTrue(self.agent.audio_available)
            self.assertTrue(self.agent.play_sound('dispatch', 1))
        self.assertIn('AUDIO player paplay failed for dispatch: device busy', (self.agent.home / 'agent.log').read_text())

    @staticmethod
    def _deny_console(path, *args, **kwargs):
        import io
        if str(path) == '/dev/console':
            raise PermissionError('no console')
        return io.open(path, *args, **kwargs)

    def test_generated_audio_migration_and_conflict(self):
        """Move legacy files while preserving both copies of conflicting names."""
        old = self.agent.home / 'sound-steampunk-dispatch.wav'
        old.write_bytes(b'old')
        demo = self.agent.home / 'demo-steampunk.wav'
        demo.write_bytes(b'demo')
        conflict = self.agent.home / 'sound-steampunk-start.wav'
        conflict.write_bytes(b'legacy')
        self.agent.audio_dir.mkdir()
        (self.agent.audio_dir / conflict.name).write_bytes(b'existing')
        self.agent.prepare_generated_audio()
        self.assertFalse(old.exists())
        self.assertEqual((self.agent.audio_dir / old.name).read_bytes(), b'old')
        self.assertEqual((self.agent.audio_dir / demo.name).read_bytes(), b'demo')
        self.assertEqual(conflict.read_bytes(), b'legacy')
        self.assertEqual((self.agent.audio_dir / conflict.name).read_bytes(), b'existing')

    def test_generated_audio_directory_failure_falls_back(self):
        """Inaccessible audio directories do not stop the agent."""
        with patch.object(Path, 'mkdir', side_effect=PermissionError('denied')):
            self.agent.prepare_generated_audio()
        self.assertEqual(self.agent.audio_dir, self.agent.home)

    def fixture(self, fail=False):
        """Create build scripts that produce tiny executable jobs; fail selects perf failure."""
        poll = patch.object(self.agent, 'poll_remote')
        poll.start()
        self.addCleanup(poll.stop)
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

        inv = self.root / 'utils' / 'inv_list'
        inv.mkdir(parents=True)
        lines = ['#!/bin/bash', 'set -eu', 'cd "$(dirname "$0")"']
        for name in self.agent.inv_builds():
            lines += ['touch ' + name, 'chmod +x ' + name]
        (inv / 'build_all.sh').write_text('\n'.join(lines) + '\n')

        samples = self.root / 'samples'
        samples.mkdir()
        (samples / 'Makefile').write_text('all:\n\t@true\n')

    def test_stress_failure_operator_report(self):
        """Expose the failed group, its assertion, and quoted copy paths."""
        work = self.root / 'stress_debug'
        work.mkdir()
        diagnostic = work / 'group logs' / 'svser.log'
        diagnostic.parent.mkdir()
        diagnostic.write_text("Assertion allocator size failed\n")
        runner = work / 'runner.log'
        runner.write_text('START  -svser\nFAIL   -svser (exit 134, 00:00:13)\n'
                          '       Details: ' + str(diagnostic) + '\n')
        report = self.agent.stress_failure_details(runner)
        self.assertIn('Stress variant: stress_debug', report)
        self.assertIn('Failed group: -svser (exit 134', report)
        self.assertIn('Assertion allocator size failed', report)
        self.assertIn("COPY: cp -- '", report)
        self.assertIn('TAgent-stress_debug-svser.log', report)
        diagnostic.unlink()
        self.assertIn('Cannot read diagnostic log', self.agent.stress_failure_details(runner))
        runner.unlink()
        self.assertIn('Cannot read runner log', self.agent.stress_failure_details(runner))

    def test_external_log_storage(self):
        """Create nested storage and diagnose invalid, denied and full locations."""
        destination = self.root / 'external drive' / 'TAgent-logs'
        ta.prepare_log_directory(destination, 0)
        self.assertEqual(list(destination.iterdir()), [])
        with patch('sys.argv', ['tagent', '--log-dir', str(destination)]):
            self.assertEqual(ta.arguments().state_dir, destination)
        with patch.dict(os.environ, {'LOG_DIR': str(destination)}), patch('sys.argv', ['tagent']):
            self.assertEqual(ta.arguments().state_dir, destination)
        file = self.root / 'file'
        file.write_text('not a directory')
        with self.assertRaisesRegex(ta.Failure, 'Cannot create or write log directory'):
            ta.prepare_log_directory(file, 0)
        with patch.object(ta.shutil, 'disk_usage', return_value=type('Usage', (), {'free': 0})()):
            with self.assertRaisesRegex(ta.Failure, 'requires 5.00 GB'):
                ta.prepare_log_directory(destination, 5)
        with patch.object(ta.tempfile, 'TemporaryDirectory', side_effect=PermissionError('denied')):
            with self.assertRaisesRegex(ta.Failure, 'denied'):
                ta.prepare_log_directory(destination, 0)

    def test_inverted_list_build_script_variants_and_legacy_names(self):
        """Exercise the actual script's x86 and ARM matrices with a fake compiler."""
        root = self.root / 'inv-build'
        root.mkdir()
        script = root / 'build_all.sh'
        script.write_bytes((Path(__file__).resolve().parents[2] / 'utils/inv_list/build_all.sh').read_bytes())
        tools = root / 'tools'
        tools.mkdir()
        (tools / 'uname').write_text('#!/bin/sh\ncase "$1" in -s) echo "$MOCK_OS";; -m) echo "$MOCK_ARCH";; esac\n')
        (tools / 'make').write_text('#!/bin/sh\nprintf "%s\n" "$*" > bminv\n')
        for name in ('uname', 'make'):
            (tools / name).chmod(0o755)
        env = dict(os.environ, PATH=str(tools) + os.pathsep + os.environ['PATH'],
                   MOCK_OS='Linux', MOCK_ARCH='x86_64')
        subprocess.run(['bash', str(script), '--avx512'], env=env, check=True,
                       stdout=subprocess.DEVNULL)
        for suffix in ('sse2', 'sse42', 'avx2', 'avx512'):
            for mode in ('debug', 'release'):
                self.assertTrue((root / ('bminv_' + mode + '_' + suffix)).is_file())
            self.assertEqual((root / ('bminv_' + suffix)).read_bytes(),
                             (root / ('bminv_release_' + suffix)).read_bytes())
        self.assertEqual((root / 'bminv_reg').read_bytes(), (root / 'bminv_release').read_bytes())
        env.update(MOCK_OS='Darwin', MOCK_ARCH='arm64')
        subprocess.run(['bash', str(script)], env=env, check=True, stdout=subprocess.DEVNULL)
        self.assertTrue((root / 'bminv_debug_neon').is_file())
        self.assertEqual((root / 'bminv_neon').read_bytes(), (root / 'bminv_release_neon').read_bytes())
        result = subprocess.run(['bash', str(script), '--avx512'], env=env,
                                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        self.assertNotEqual(result.returncode, 0)

    def test_inverted_list_build_failure_and_matrix(self):
        """Track all variant binaries and stop before tests when the utility fails."""
        self.fixture()
        self.assertEqual(self.agent.inv_builds(), ['bminv_debug', 'bminv_release',
                         'bminv_debug_neon', 'bminv_release_neon'])
        (self.root / 'utils/inv_list/build_all.sh').write_text(
            '#!/bin/bash\necho "inv.cpp:1: warning: sample"\nexit 1\n')
        with self.assertRaises(ta.Failure) as error:
            self.agent.cycle()
        self.assertIn('build-inv_list.log', str(error.exception))
        self.assertEqual(self.agent.build_warnings.count, 1)
        run = next(self.agent.home.glob('*abc123'))
        self.assertIn('sample', (run / 'build-inv_list.log').read_text())
        self.assertFalse(list(run.glob('perf_*/perf.log')))
        self.agent.arm = False
        self.agent.args.avx512 = True
        for suffix in ('sse2', 'sse42', 'avx2', 'avx512'):
            for mode in ('debug', 'release'):
                self.assertIn('bminv_' + mode + '_' + suffix, self.agent.inv_builds())

    def test_samples_failure_and_warnings_stop_cycle(self):
        """Sample diagnostics are counted and a build failure prevents tests."""
        self.fixture()
        (self.root / 'samples' / 'Makefile').write_text(
            'all:\n\t@echo "sample.cpp:3: warning: sample warning"\n\t@false\n')
        with self.assertRaises(ta.Failure) as error:
            self.agent.cycle()
        self.assertIn('build-samples.log', str(error.exception))
        self.assertEqual(self.agent.sample_warnings.count, 1)
        self.assertEqual(self.agent.build_warnings.count, 1)
        run = next(self.agent.home.glob('*abc123'))
        self.assertIn('sample warning', (run / 'build-samples.log').read_text())
        self.assertFalse(list(run.glob('perf_*/perf.log')))

    def test_remote_checked_before_perf_after_build(self):
        """A new tip at the build boundary prevents any perf task from starting."""
        self.fixture()
        self.agent.next_remote_poll = float('inf')
        def check():
            self.assertEqual(self.agent.next_remote_poll, 0)
            self.agent.pending_tip = 'd' * 40
        self.agent.poll_remote.side_effect = check
        self.assertEqual(self.agent.cycle(), 'updated')
        run = next(self.agent.home.glob('*abc123'))
        self.assertFalse(list(run.glob('perf_*/perf.log')))
        self.assertTrue(self.agent.args.force)

    def test_success_and_isolation(self):
        """All expected builds run, including debug perf, with isolated stress files."""
        self.fixture()
        with patch.object(self.agent, 'start_test_sound', wraps=self.agent.start_test_sound) as departure:
            self.assertTrue(self.agent.cycle())
            self.assertEqual([call.args[0] for call in departure.call_args_list],
                             ['perf', 'stress'])
        self.assertEqual(self.agent.next_train_at, float('inf'))
        runs = list(self.agent.home.glob('*abc123'))
        self.assertEqual(len(runs), 1)
        self.assertTrue((runs[0] / 'perf_debug' / 'perf.log').exists())
        self.assertTrue((runs[0] / 'perf_debug_neon' / 'perf.log').exists())
        self.assertEqual(len(list(runs[0].glob('stress_*/job.*/isolated-output'))), 8)
        self.assertTrue((runs[0] / 'manifest.json').exists())

    def test_restart_on_update_interrupts_build_and_preserves_logs(self):
        """An update during a build returns to fetch without treating it as a failure."""
        self.fixture()
        self.agent.args.restart_on_update = True
        self.agent.pending_tip = 'd' * 40
        with patch.object(self.agent, 'command', return_value='skipped-update'):
            self.assertEqual(self.agent.cycle(), 'updated')
        self.assertTrue(list(self.agent.home.glob('*abc123/manifest.json')))
        self.assertFalse((self.agent.home / 'last-test-failure.txt').exists())

    def test_restart_on_update_stops_active_test(self):
        """A running process is terminated after a remote update heartbeat."""
        self.agent.args.restart_on_update = True
        def found_update():
            self.agent.pending_tip = 'd' * 40
        with patch.object(self.agent, 'poll_remote', side_effect=found_update):
            status = self.agent.command(['bash', '-c', 'sleep 30'], self.root / 'perf.log')
        self.assertEqual(status, 'skipped-update')
        self.assertFalse(self.agent.children)

    def test_restart_on_update_skips_new_work(self):
        """No new test or build process starts after a restart is requested."""
        self.agent.args.restart_on_update = True
        self.agent.pending_tip = 'd' * 40
        for name in ('perf.log', 'runner.log', 'build-perf.log'):
            log = self.root / name
            self.assertEqual(self.agent.command(['/definitely/not/a/program'], log), 'skipped-update')

    def test_stress_debug_runs_last_and_fast_in_pedantic_mode(self):
        """Release variants precede fast debug variants without conflicting flags."""
        self.fixture()
        self.agent.args.pedantic = True
        calls = []
        original = self.agent.command
        def command(argv, *args, **kwargs):
            if '--order-file' in argv:
                calls.append((Path(argv[1]).parent.name, list(argv)))
            return original(argv, *args, **kwargs)
        with patch.object(self.agent, 'command', side_effect=command):
            self.assertTrue(self.agent.cycle())
        names = [name for name, _ in calls]
        self.assertTrue(names)
        self.assertEqual(names, sorted(names, key=lambda name: '_debug' in name))
        self.assertTrue(any('_debug' in name for name in names))
        for name, argv in calls:
            self.assertEqual('--fast' in argv, '_debug' in name)
            self.assertEqual('--pedantic' in argv, '_debug' not in name)

    def test_stress_transition_announces_each_success_twice(self):
        """Transitions name actual consecutive runners, repeat, and finish cleanly."""
        self.fixture()
        self.agent.args.no_debug = True
        with patch.object(self.agent, 'sound') as sound:
            self.assertTrue(self.agent.cycle())
        messages = [c.kwargs['workers'] for c in sound.call_args_list if c.args[0] == 'stress-transition']
        self.assertEqual(len(messages), 3)
        self.assertIn('stress run is starting', messages[0])
        self.assertEqual(messages[1].count('stress run completed successfully.'), 2)
        self.assertEqual(messages[1].count('stress run is next.'), 2)
        self.assertIn('All stress runs are complete.', messages[-1])
        self.assertEqual(messages[-1].count('stress run completed successfully.'), 2)

    def test_reference_discovery_rejects_partial_copy(self):
        """A manifest alone must not make a partial copy eligible."""
        import csv
        candidate = self.agent.home / 'gov2-reference-fixture'
        candidate.mkdir()
        rows = []
        for dataset in ('sorted', 'unsorted'):
            for level in range(1, 11):
                rows.append({'dataset': 'Gov2.' + dataset, 'compression_level': level,
                             'stream_file': 'gov2.{}.level-{:02d}.bv'.format(dataset, level),
                             'compressed_bytes': 1, 'index_bytes': 1, 'input_bytes': 1})
        with (candidate / 'references.tsv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=rows[0], delimiter='\t')
            writer.writeheader()
            writer.writerows(rows)
        self.assertFalse(ta.inv_reference.inspect_reference(candidate)[0])
        (candidate / 'Gov2Flat').mkdir()
        for dataset in ('sorted', 'unsorted'):
            (candidate / 'Gov2Flat' / ('gov2.' + dataset)).write_bytes(b'x')
        for row in rows:
            file = candidate / row['stream_file']
            file.write_bytes(b'x')
            file.with_suffix('.idx').write_bytes(b'x')
        self.assertTrue(ta.inv_reference.inspect_reference(candidate)[0])
        self.assertFalse(ta.inv_reference.inspect_reference(candidate, full=True)[0])
        (candidate / 'matrix.json').write_text('{}')
        self.assertTrue(ta.inv_reference.inspect_reference(candidate, full=True)[0])

    def test_inv_fast_selects_one_best_supported_variant(self):
        """Fast mode uses one supported release binary and the CPU budget."""
        self.agent.args.fast = True
        self.agent.inv_reference = self.root / 'reference'
        self.agent.flags = {'sse2', 'sse4_2', 'avx2'}
        self.agent.arm = False
        self.agent.next_remote_poll = float('inf')
        manifest = {'builds': {'inv_list': self.agent.inv_builds()}}
        with patch.object(self.agent, 'command', return_value='') as command, \
             patch.object(self.agent, 'sound'), patch.object(self.agent, 'update_waiting', return_value=False):
            self.assertTrue(self.agent.test_inv_stage(self.root, manifest))
        args = command.call_args.args[0]
        self.assertEqual(Path(args[args.index('--bin') + 1]).name, 'bminv_release_avx2')
        self.assertIn('--fast', args)
        self.assertIn('--silent', args)
        self.assertEqual(args[args.index('--threads') + 1], '2')

    def test_inv_reader_settings_are_forwarded_and_logged(self):
        """Both validation modes receive the read budget and selected backend."""
        self.agent.inv_reference = self.root / 'reference'
        self.agent.flags = {'avx2'}
        manifest = {'builds': {'inv_list': ['bminv_release_avx2']}}
        for fast, mmap, budget in ((True, True, 8192), (False, False, 256)):
            self.agent.args.fast = fast
            self.agent.args.inv_read_mmap = mmap
            self.agent.args.inv_read_memory_mib = budget
            run = self.root / ('fast' if fast else 'normal')
            with patch.object(self.agent, 'command', return_value='') as command, \
                 patch.object(self.agent, 'sound'), patch.object(self.agent, 'update_waiting', return_value=False):
                self.assertTrue(self.agent.test_inv_stage(run, manifest))
            argv = command.call_args.args[0]
            self.assertEqual(argv[argv.index('--read-memory-mib') + 1], str(budget))
            self.assertEqual('--read-mmap' in argv, mmap)
            self.assertEqual('--fast' in argv, fast)
            self.assertEqual(argv[argv.index('--start-level') + 1], '6')
            self.assertEqual(argv[argv.index('--golden-level') + 1], '6')
            import json
            saved = json.loads((run / 'inv_list/bminv_release_avx2/command.json').read_text())
            self.assertEqual(saved, argv)

    def test_inv_read_budget_uses_linux_available_ram(self):
        """Automatic sizing respects available RAM, cgroups, caps and explicit overrides."""
        self.agent.args.inv_read_memory_mib = None
        mib = 1024 * 1024
        with patch.object(ta.platform, 'system', return_value='Linux'), \
             patch.object(ta.memory_monitor.MemoryMonitor, 'sample') as sample:
            for memory, expected in (({'available': 4000 * mib}, 2800),
                                     ({'available': 24000 * mib}, 16800),
                                     ({'available': 16000 * mib, 'cgroup_limit': 8000 * mib,
                                       'cgroup_current': 6000 * mib}, 1400)):
                sample.return_value = memory
                self.assertEqual(self.agent.inv_read_budget(), expected)
            sample.side_effect = OSError('unavailable')
            self.assertEqual(self.agent.inv_read_budget(), 512)
            self.agent.args.inv_read_memory_mib = 1024
            self.assertEqual(self.agent.inv_read_budget(), 1024)
        self.agent.args.inv_read_memory_mib = None
        with patch.object(ta.platform, 'system', return_value='Darwin'):
            self.assertEqual(self.agent.inv_read_budget(), 8192)

    def test_inv_progress_uses_selected_case_count(self):
        """Progress and speech report the selected fast or normal case total."""
        output = self.root / 'results'
        output.mkdir()
        (output / 'summary.tsv').write_text('dataset\tlevel\tphase\tstatus\nGov2.sorted\t6\tverify\tpassed\n')
        for total in (10, 20):
            self.agent.inv_progress_time = 0
            with patch.object(self.agent, 'say') as say, patch.object(self.agent, 'sound') as sound:
                self.agent.report_inv_progress(output, 'bminv_release', total)
            self.assertIn('1 of {} cases completed'.format(total), say.call_args.args[0])
            self.assertIn('1 of {} cases completed'.format(total), sound.call_args.kwargs['workers'])

    def test_inv_reader_cli_validation(self):
        """Reader defaults, overrides and invalid budgets are resolved at startup."""
        with patch.dict(os.environ, {'INV_READ_MEMORY_MIB': '8192', 'INV_READ_MMAP': 'yes'}):
            with patch.object(sys, 'argv', ['tagent']):
                args = ta.arguments()
                self.assertEqual(args.inv_read_memory_mib, 8192)
                self.assertEqual(args.inv_read_mmap, os.name == 'posix')
            with patch.object(sys, 'argv', ['tagent', '--read-memory-mib', '1024', '--no-read-mmap']):
                args = ta.arguments()
                self.assertEqual(args.inv_read_memory_mib, 1024)
                self.assertFalse(args.inv_read_mmap)
            with patch.object(sys, 'argv', ['tagent', '--read-memory-mib', '16384']):
                self.assertEqual(ta.arguments().inv_read_memory_mib, 16384)
            for budget in ('0', str(sys.maxsize), '1.5'):
                with self.subTest(budget=budget), patch.object(sys, 'argv', ['tagent', '--read-memory-mib', budget]):
                    with self.assertRaises(SystemExit):
                        ta.arguments()

    def test_suppress_options_combine(self):
        with patch.object(sys, 'argv', ['tagent', '--suppress', 'perf,inv_list', '--suppress', 'stress']):
            self.assertEqual(ta.arguments().suppress, {'perf', 'inv_list', 'stress'})

    def test_no_debug_skips_execution_but_keeps_builds(self):
        """Pedantic release runs remain enabled while debug binaries are built."""
        self.fixture()
        self.agent.args.no_debug = True
        self.agent.args.pedantic = True
        self.assertTrue(self.agent.cycle())
        run = next(self.agent.home.glob('*abc123'))
        self.assertTrue((self.root / 'tests/perf/perf_debug').exists())
        self.assertTrue((self.root / 'tests/stress/stress_debug').exists())
        self.assertFalse((run / 'perf_debug').exists())
        self.assertFalse((run / 'stress_debug').exists())
        self.assertTrue((run / 'perf_release/perf.log').exists())
        self.assertTrue((run / 'stress_release/runner.log').exists())

    def test_fast_and_pedantic_are_rejected(self):
        """Conflicting workload modes fail before an agent starts."""
        with patch.object(sys, 'argv', ['tagent', '--fast', '--pedantic']):
            with self.assertRaises(SystemExit) as error:
                ta.arguments()
        self.assertEqual(error.exception.code, 2)

    def test_fast_applies_to_all_stress(self):
        """Express tests use fast mode; normal release tests reuse the builds."""
        self.fixture()
        self.agent.args.fast = True
        calls = []
        original = self.agent.command
        def command(argv, *args, **kwargs):
            calls.append((list(argv), args[0] if args else None))
            return original(argv, *args, **kwargs)
        with patch.object(self.agent, 'command', side_effect=command):
            self.assertTrue(self.agent.cycle())
        stress = [(argv, log) for argv, log in calls if '--order-file' in argv]
        self.assertTrue(stress)
        fast_stress = [log.parent.name for _, log in stress if 'normal' not in log.parts]
        self.assertEqual(fast_stress, ['stress_release', 'stress_release_neon'])
        for argv, log in stress:
            normal = 'normal' in log.parts
            debug = '_debug' in log.parent.name
            self.assertEqual('--fast' in argv, not normal or debug)
            self.assertNotIn('--pedantic', argv)
        perf_runs = [(argv, log) for argv, log in calls if Path(argv[0]).name.startswith('perf_release')]
        self.assertTrue(perf_runs)
        fast_perf = {log.parent.name for _, log in perf_runs if 'normal' not in log.parts}
        self.assertEqual(fast_perf, {'perf_release', 'perf_release_neon'})
        for argv, log in perf_runs:
            self.assertEqual('--fast' in argv, 'normal' not in log.parts)
        perf_builds = [argv for argv, _ in calls if len(argv) > 1 and str(argv[1]).endswith('/perf/build_all.sh')]
        self.assertEqual(len(perf_builds), 1)
        self.assertNotIn('--pedantic', perf_builds[0])
        self.assertNotIn('--fast', perf_builds[0])

    def test_fast_variants_by_platform(self):
        """Fast selects baseline and supported SIMD, never debug or AVX-512."""
        names = ['perf_release', 'perf_debug', 'perf_release_neon',
                 'perf_release_sse2', 'perf_release_avx2', 'perf_release_avx512']
        self.agent.arm = True
        self.assertEqual(self.agent.fast_test_variants('perf', names),
                         ['perf_release', 'perf_release_neon'])
        self.agent.arm = False
        self.agent.flags = {'avx2', 'sse2', 'avx512f', 'avx512bw', 'avx512dq', 'avx512vl'}
        self.assertEqual(self.agent.fast_test_variants('perf', names),
                         ['perf_release', 'perf_release_avx2'])
        self.agent.flags = {'sse2'}
        self.assertEqual(self.agent.fast_test_variants('perf', names), ['perf_release'])
        with self.assertRaises(ta.Failure):
            self.agent.fast_test_variants('perf', ['perf_release_avx2'])

    def test_perf_debug_is_always_fast(self):
        """Regular runs use fast for every debug perf binary but not release."""
        self.fixture()
        calls = []
        original = self.agent.command
        def command(argv, *args, **kwargs):
            log = args[0] if args else kwargs.get('log')
            if log and log.name == 'perf.log':
                calls.append(list(argv))
            return original(argv, *args, **kwargs)
        with patch.object(self.agent, 'command', side_effect=command):
            self.assertTrue(self.agent.cycle())
        debug = [argv for argv in calls if '_debug' in Path(argv[0]).name]
        release = [argv for argv in calls if '_release' in Path(argv[0]).name]
        self.assertEqual(len(debug), 2)
        self.assertEqual(len(release), 2)
        self.assertTrue(all(argv[1:] == ['--fast'] for argv in debug))
        self.assertTrue(all(len(argv) == 1 for argv in release))

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

    @unittest.skipUnless(os.environ.get("TAGENT_FAULT_SIGNALS") == "1", "opt-in signal faults can show OS crash dialogs")
    def test_fault_signals_stop_siblings(self):
        """Silent crashes preserve the signal and terminate a sibling job."""
        import signal
        for sig in (signal.SIGSEGV, signal.SIGFPE, signal.SIGABRT):
            with self.subTest(signal=sig):
                self.agent.stopping.clear()
                if hasattr(self.agent, 'first_workload_failure'):
                    del self.agent.first_workload_failure
                marker = self.root / ('survivor-' + str(sig))
                def sibling():
                    try:
                        self.agent.command(['bash', '-c', 'sleep 8; touch "$1"', 'worker', str(marker)],
                                           self.root / 'sibling.log')
                    except ta.Failure:
                        pass
                thread = threading.Thread(target=sibling)
                thread.start()
                time.sleep(.1)
                log = self.root / ('fault-' + str(sig) + '.log')
                with self.assertRaisesRegex(ta.CommandFailure, signal.Signals(sig).name):
                    self.agent.command([sys.executable, '-c',
                        'import os,resource; resource.setrlimit(resource.RLIMIT_CORE,(0,0)); os.kill(os.getpid(), ' + str(sig) + ')'], log)
                thread.join(timeout=6)
                self.assertFalse(thread.is_alive())
                self.assertFalse(marker.exists())
                self.assertEqual(log.read_text(), '')

    def test_main_exits_for_failure_with_pending_update(self):
        """A pending commit never turns a fatal workload result into a restart."""
        self.agent.pending_tip = 'new'
        with patch.object(sys, 'argv', ['tagent']), \
             patch.object(ta, 'Agent', return_value=self.agent), \
             patch.object(self.agent, 'run', side_effect=ta.CommandFailure('perf_release SIGFPE log /tmp/perf.log')):
            self.assertEqual(ta.main(), 1)
        self.assertIn('SIGFPE', (self.agent.home / 'agent.log').read_text())

    def test_stress_fail_fast_kills_unresponsive_sibling(self):
        """The standalone stress runner stops a TERM-ignoring sibling after exit 17."""
        import shutil
        runner = self.root / 'run_all_pl.sh'
        shutil.copy2(Path(__file__).resolve().parents[1] / 'stress/run_all_pl.sh', runner)
        worker = self.root / 'worker'
        worker.write_text('#!' + sys.executable + '\nimport sys,time,signal\n'
                          'if "-ll" in sys.argv: time.sleep(.3); sys.exit(17)\n'
                          'signal.signal(signal.SIGTERM, signal.SIG_IGN)\ntime.sleep(30)\n')
        worker.chmod(0o755)
        result = subprocess.run(['bash', str(runner), '--binary', './worker',
                                 '--jobs', '2', '--fail-fast', '-ll', '-s'],
                                cwd=self.root, capture_output=True, text=True, timeout=8)
        self.assertEqual(result.returncode, 1, result.stdout + result.stderr)
        self.assertIn('-ll exit 17', (self.root / 'first-failure.txt').read_text())
        self.assertNotIn('All 2 test groups passed', result.stdout)

    def test_unwritable_workload_log_is_fatal(self):
        """A missing log parent enters fatal supervision instead of leaving jobs running."""
        log = self.root / 'missing-parent' / 'perf.log'
        with self.assertRaisesRegex(ta.CommandFailure, 'perf.log'):
            self.agent.command([sys.executable, '-c', 'pass'], log)
        self.assertTrue(self.agent.stopping.is_set())
        self.assertFalse(self.agent.children)

    def test_crash_wins_update_boundary(self):
        """A failed binary cannot be treated as update cancellation."""
        proc = subprocess.Popen([sys.executable, '-c', 'raise SystemExit(17)'])
        proc.wait()
        self.agent.pending_tip = 'new'
        with self.assertRaisesRegex(ta.CommandFailure, '17'):
            self.agent.check_update_failure(proc, self.root / 'perf.log')

    def test_ambience_tracks_test_mode(self):
        """Express and pedantic keep their train loop; normal selects the new loop."""
        from unittest.mock import Mock
        self.agent.depot = Mock()
        for fast, pedantic, expected in ((True, False, 'train'),
                                          (False, True, 'train'),
                                          (False, False, 'train-normal')):
            with self.subTest(fast=fast, pedantic=pedantic):
                self.agent.args.fast = fast
                self.agent.args.pedantic = pedantic
                self.agent.start_test_sound('stress')
                self.agent.depot.set_mode.assert_called_with(expected)

    def test_train_stop_preview_includes_conductor(self):
        """The stop preview queues the same kind of spoken transition as a run."""
        self.agent.args.sound_preview = 'train-stop'
        self.agent.args.no_sound = True
        with patch.object(self.agent, 'sound') as sound:
            self.assertEqual(self.agent.run(), 0)
        calls = sound.call_args_list
        self.assertEqual([call.args[0] for call in calls], ['train-stop', 'stress-transition'])
        self.assertIn('SSE two stress run completed successfully', calls[1].kwargs['workers'])
        self.assertIn('AVX two stress run is next', calls[1].kwargs['workers'])

    def test_sound_demo_includes_normal_ambience(self):
        """The show includes a short excerpt without requiring a sound device."""
        self.agent.args.sound_demo = True
        self.agent.args.no_sound = True
        with patch.object(self.agent, 'sound') as sound:
            self.assertEqual(self.agent.run(), 0)
        events = [call.args[0] for call in sound.call_args_list]
        self.assertIn('train-normal-demo', events)
        self.assertIn('train-stop', events)
        self.assertIn('depart-perf', events)
        self.assertIn('fail', events)

    def test_fast_cycle_upgrades_to_normal_without_rebuilding(self):
        """A clear remote check runs normal stages on the existing binaries."""
        self.fixture()
        self.agent.args.fast = True
        observed = []
        def stage(run, manifest):
            observed.append((self.agent.args.fast, run.name))
            return True
        with patch.object(self.agent, 'test_perf_stage', side_effect=stage), \
             patch.object(self.agent, 'test_stress_stage', side_effect=stage), \
             patch.object(self.agent, 'update_waiting', return_value=False), \
             patch.object(self.agent, 'sound') as sound:
            self.assertTrue(self.agent.cycle())
        self.assertEqual([mode for mode, _ in observed], [True, True, False, False])
        self.assertEqual([name for _, name in observed[2:]], ['normal', 'normal'])
        self.assertTrue(self.agent.args.fast)
        events = [call.args[0] for call in sound.call_args_list]
        start = events.index('upgrade-dispatch')
        self.assertEqual(events[start:start + 4],
                         ['upgrade-dispatch', 'upgrade-conductor', 'upgrade-engineer', 'upgrade-depart'])
        self.assertEqual(events.count('train-stop'), 1)

    def test_fast_cycle_skips_normal_when_update_arrives(self):
        """An update at the upgrade boundary takes priority over normal testing."""
        self.fixture()
        self.agent.args.fast = True
        observed = []
        def stage(run, manifest):
            observed.append(self.agent.args.fast)
            return True
        with patch.object(self.agent, 'test_perf_stage', side_effect=stage), \
             patch.object(self.agent, 'test_stress_stage', side_effect=stage), \
             patch.object(self.agent, 'update_waiting', side_effect=[False, True]), \
             patch.object(self.agent, 'sound') as sound:
            self.assertEqual(self.agent.cycle(), 'updated')
        self.assertEqual(observed, [True, True])
        self.assertNotIn('upgrade-dispatch', [call.args[0] for call in sound.call_args_list])

    def test_test_failure_is_fatal(self):
        """A nonzero perf result returns failure with logs rather than passing."""
        self.fixture(fail=True)
        with self.assertRaises(ta.CommandFailure):
            self.agent.cycle()
        self.assertTrue(list(self.agent.home.glob('*abc123/perf_release/perf.log')))
        marker = self.agent.home / 'active-cycle.json'
        self.assertTrue(marker.exists())
        self.assertIn('perf', marker.read_text())

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
            sound.assert_any_call('build-final-failed', workers=(2, self.agent.build_warnings.count, False))
        self.assertIn('BUILD FAILED | diagnostic errors: 2 | warnings:',
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
        """A one-worker budget serializes all four build families."""
        self.fixture()
        self.agent.jobs = 1
        active = [0, 0]
        mutex = threading.Lock()
        original = self.agent.command
        def command(argv, log=None, *args, **kwargs):
            """Record overlap while the synthetic builds run."""
            if not log or log.name not in ('build-perf.log', 'build-stress.log', 'build-inv_list.log', 'build-samples.log'):
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
        with patch.object(self.agent, 'update_waiting', side_effect=[False, True]):
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

    def test_joe_replies_after_first_completed_family(self):
        """Queue Dispatch before one Joe reply and avoid repeats later in the cycle."""
        log = self.root / 'build-perf.log'
        log.write_text('')
        target = self.root / 'perf'
        target.write_text('ready')
        target.chmod(0o755)
        self.agent.build_targets = {'perf': [target], 'stress': [self.root / 'pending']}
        self.agent.next_build_report = 0
        with patch.object(self.agent, 'sound') as sound:
            self.agent.build_progress(log)
            self.assertEqual([x.args[0] for x in sound.call_args_list], ['build-first-clear', 'engineer-build-reply'])
            self.agent.next_build_report = 0
            self.agent.build_progress(log)
            self.agent.finish_build(True)
            self.assertEqual(sum(x.args[0] == 'engineer-build-reply' for x in sound.call_args_list), 2)

    def test_engineer_replies_follow_successful_build_without_repeats(self):
        """Dispatch precedes Joe, who does not joke after a build failure."""
        with patch.object(self.agent, 'sound') as sound:
            for _ in ta.ENGINEER_BUILD_REPLIES:
                self.agent.finish_build(True)
            calls = sound.call_args_list
            self.assertEqual([c.args[0] for c in calls[::2]], ['build-final-ok'] * len(ta.ENGINEER_BUILD_REPLIES))
            replies = [c.kwargs['workers'] for c in calls[1::2]]
            self.assertEqual(set(replies), set(ta.ENGINEER_BUILD_REPLIES))
            sound.reset_mock()
            self.agent.finish_build(False)
            sound.assert_called_once()
            self.assertEqual(sound.call_args.args[0], 'build-final-failed')

    def test_samples_spoken_completion_and_global_warnings(self):
        """Sample completion is explicit; warnings aggregate across all families."""
        logs = []
        for family in ('perf', 'stress', 'inv_list', 'samples'):
            path = self.root / ('build-' + family + '.log')
            path.write_text('file.cpp:1: warning: example\n')
            logs.append(path)
        self.agent.build_logs = logs
        self.agent.build_targets = {'perf': [], 'stress': [], 'inv_list': [], 'samples': []}
        self.agent.next_build_report = 0
        self.agent.samples_complete = True
        with patch.object(self.agent, 'sound') as sound:
            self.agent.build_progress(logs[0])
        report = sound.call_args_list[0].kwargs['workers']
        self.assertEqual(report[0], 4)
        self.assertEqual(self.agent.sample_warnings.count, 1)
        self.assertIn('Sample build is complete.', ta.build_report_text(report[0], report[1:]))
        self.agent.samples_complete = False
        self.agent.next_build_report = 0
        with patch.object(self.agent, 'sound') as sound:
            self.agent.build_progress(logs[0])
        report = sound.call_args.kwargs['workers']
        self.assertIn('Sample build is in progress.', ta.build_report_text(report[0], report[1:]))
        self.assertEqual(report[0], 4)

    def test_spoken_build_completion_and_pending(self):
        """Completed targets use concise speech while in-progress counts remain."""
        text = ta.build_report_text(0, (6, 0, 2, 3, 0, 8, 0))
        self.assertIn('Performance build is complete.', text)
        self.assertIn('Stress builds: 2 ready, 3 pending.', text)
        self.assertIn('Inverted list build is complete.', text)
        self.assertIn('Performance builds: 0 ready, 0 pending.', ta.build_report_text(0, (0, 0, 1, 1)))

    def test_final_debug_uses_successful_release_reference(self):
        """Use matching release counts and issue the mode disclaimer once."""
        debug = self.root / 'perf_debug_avx2' / 'perf.log'
        release = self.root / 'perf_release_avx2' / 'perf.log'
        debug.parent.mkdir()
        release.parent.mkdir()
        debug.write_text('header\na; 1 sec\n')
        release.write_text('a; 1 sec\nb; 2 sec\n')
        self.agent.perf_logs = {1: debug}
        self.agent.perf_pending = {'perf_debug_avx2'}
        self.agent.perf_completed = {'perf_release_avx2': release}
        with patch.object(self.agent, 'sound') as sound:
            self.agent.next_perf_report = 0
            self.agent.report_perf_progress()
            report = sound.call_args.kwargs['workers']
            self.assertEqual(report, ('perf_debug_avx2', 1, 'perf_release_avx2', 2, True, True))
            self.agent.next_perf_report = 0
            debug.write_text('a; 1 sec\nb; 2 sec\nc; 3 sec\n')
            self.agent.report_perf_progress()
            self.assertEqual(sound.call_args.kwargs['workers'], ('perf_debug_avx2', 3, 'perf_release_avx2', 2, False, False))
            with patch.object(self.agent, 'synthesize_announcement', return_value=1) as speak:
                self.agent.render_perf_progress(self.root / 'out.wav', sound.call_args.kwargs['workers'])
                self.assertIn('delayed', speak.call_args.args[1])

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
            log = (self.root / 'interrupt-state' / 'agent.log').read_text()
            self.assertIn('SIGNAL SIGINT (2) received', log)
            self.assertIn('sender unavailable', log)
            self.assertIn('pid=', log)
            self.assertIn('ppid=', log)
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait()

    def test_termination_signals_are_reported(self):
        """SIGTERM and SIGHUP retain their identity and conventional exit codes."""
        import signal
        for sig in (signal.SIGTERM, signal.SIGHUP):
            state = self.root / sig.name
            ready = self.root / (sig.name + '.ready')
            code = (
                'import sys,time; sys.path.insert(0, ' + repr(str(Path(__file__).parent.resolve())) + '); '
                'import tagent; '
                'sys.argv=["tagent", "--no-sound", "--state-dir", ' + repr(str(state)) + ']; '
                'tagent.Agent.run=lambda self: (open(' + repr(str(ready)) + ', "w").close(), time.sleep(30)); '
                'sys.exit(tagent.main())')
            proc = subprocess.Popen([sys.executable, '-c', code], stdout=subprocess.DEVNULL)
            try:
                deadline = time.monotonic() + 5
                while not ready.exists() and time.monotonic() < deadline:
                    time.sleep(.01)
                self.assertTrue(ready.exists())
                proc.send_signal(sig)
                self.assertEqual(proc.wait(timeout=10), 128 + sig)
                self.assertIn('SIGNAL ' + sig.name, (state / 'agent.log').read_text())
            finally:
                if proc.poll() is None:
                    proc.kill()
                    proc.wait()

    def test_concurrent_log_entries_remain_intact(self):
        """Concurrent workers cannot interleave lines in the shared agent log."""
        def write(worker):
            for number in range(30):
                self.agent.say('worker={} entry={}'.format(worker, number))
        threads = [threading.Thread(target=write, args=(n,)) for n in range(8)]
        with patch('builtins.print'):
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
        lines = (self.agent.home / 'agent.log').read_text().splitlines()
        self.assertEqual(len(lines), 240)
        self.assertEqual(len(set(line[20:] for line in lines)), 240)

    def test_pedantic_matrix(self):
        """x86 pedantic perf includes debug/release for every SIMD suffix."""
        self.agent.arm = False
        self.agent.args.pedantic = True
        for suffix in ('sse2', 'sse42', 'avx2', '64'):
            for mode in ('debug', 'release'):
                self.assertIn('perf_' + mode + '_' + suffix, self.agent.builds('perf'))


if __name__ == '__main__':
    unittest.main()
