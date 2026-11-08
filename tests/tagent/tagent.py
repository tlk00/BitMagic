#!/usr/bin/env python3
"""Remote BitMagic supervisor; standard library only, Python 3.9+ and POSIX.

Owns subprocess groups and one checkout lock. Never discards repository changes.
Build/test outputs live in persistent per-cycle directories for investigation.
"""
import argparse
import concurrent.futures
import datetime
import fcntl
import json
import math
import os
from pathlib import Path
import platform
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
import soundtrack
import wave


class Failure(RuntimeError):
    """An operational failure, or a failed build/test with its recorded command."""


class CommandFailure(Failure):
    """A child command returned nonzero; test phases may wait for a fix."""


def arguments():
    """Return validated CLI configuration, defaulting to launcher environment."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in [('poll-minutes', 5), ('retry-minutes', 60),
                          ('reserve-cpus', 2), ('jobs', 0), ('min-free-gb', 5)]:
        parser.add_argument('--' + name, type=float, default=float(
            os.environ.get(name.upper().replace('-', '_'), default)))
    for name, default in [('remote', 'origin'), ('branch', 'master')]:
        parser.add_argument('--' + name, default=os.environ.get(name.upper(), default))
    parser.add_argument('--state-dir', type=Path)
    for name in ['pedantic', 'avx512', 'once', 'check', 'force']:
        parser.add_argument('--' + name, action='store_true')
    parser.add_argument('--no-sound', action='store_true',
                        default=os.environ.get('SOUND', 'yes') == 'no')
    parser.add_argument('--sound-theme', choices=soundtrack.THEMES,
                        default=os.environ.get('SOUND_THEME', 'dieselpunk'))
    parser.add_argument('--sound-volume', type=float,
                        default=float(os.environ.get('SOUND_VOLUME', '0.35')))
    parser.add_argument('--sound-demo', action='store_true',
                        help='Play a complete theme demonstration and exit without Git/builds')
    parser.add_argument('--sound-preview', choices=soundtrack.EVENTS,
                        help='Play one event locally and exit without Git/build validation')
    args = parser.parse_args()
    if args.sound_theme not in soundtrack.THEMES:
        parser.error('Unknown SOUND_THEME')
    if not math.isfinite(args.sound_volume) or not 0 <= args.sound_volume <= 1:
        parser.error('sound-volume must be between 0 and 1')
    for name in ['poll_minutes', 'retry_minutes', 'min_free_gb']:
        if not math.isfinite(getattr(args, name)) or getattr(args, name) <= 0:
            parser.error(name + ' must be positive and finite')
    for name in ['jobs', 'reserve_cpus']:
        value = getattr(args, name)
        if not math.isfinite(value) or value < 0 or value != int(value):
            parser.error(name + ' must be a nonnegative integer')
        setattr(args, name, int(value))
    if args.remote.startswith('-') or args.branch.startswith('-'):
        parser.error('remote/branch must not start with a dash')
    return args


def available_cpus():
    """Return logical CPU budget, respecting affinity and cgroup quotas if exposed."""
    count = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    try:
        quota, period = Path('/sys/fs/cgroup/cpu.max').read_text().split()
        if quota != 'max':
            count = min(count, max(1, int(quota) // int(period)))
    except (OSError, ValueError):
        try:
            root = Path('/sys/fs/cgroup/cpu')
            quota = int((root / 'cpu.cfs_quota_us').read_text())
            period = int((root / 'cpu.cfs_period_us').read_text())
            if quota > 0:
                count = min(count, max(1, quota // period))
        except (OSError, ValueError):
            pass
    return count


class Agent:
    """Owns configuration, process lifetime, logs and update/recovery state."""

    def __init__(self, args):
        """Initialize local paths and tracking; args is validated CLI configuration."""
        self.args = args
        self.repo = Path(__file__).resolve().parents[2]
        self.home = (args.state_dir or self.repo / '.tagent').resolve()
        self.home.mkdir(parents=True, exist_ok=True)
        self.children = set()
        self.mutex = threading.RLock()
        self.stopping = threading.Event()
        self.jobs = args.jobs or max(1, available_cpus() - args.reserve_cpus)
        self.last_fetch = 'never'
        self.attempts = 0
        self.commit = 'unknown'
        self.last_sound = 0
        self.next_train_at = math.inf
        self.audio_thread = None
        self.audio_mutex = threading.Lock()
        self.audio_warned = False
        self.perf_workers = set()
        self.lock_file = None
        self.flags = set()
        self.arm = platform.system() == 'Darwin' and platform.machine() in ('arm64', 'aarch64')

    def say(self, message):
        """Print and persist a timestamped lifecycle message."""
        line = time.strftime('%Y-%m-%d %H:%M:%S ') + message
        print('\r' + line, flush=True)
        with (self.home / 'agent.log').open('a') as stream:
            stream.write(line + '\n')

    def resource_check(self):
        """Raise Failure if configured storage reserve is unavailable."""
        if shutil.disk_usage(self.home).free < self.args.min_free_gb * 1024**3:
            raise Failure('Free disk below configured reserve: ' + str(self.home))
        if shutil.disk_usage(self.repo).free < self.args.min_free_gb * 1024**3:
            raise Failure('Free disk below configured reserve: ' + str(self.repo))

    def signal_group(self, proc, sig):
        """Signal an owned group; tolerate vanished groups after their leader exits.

        Some macOS sandboxes report EPERM rather than ESRCH for an empty group.
        A permission failure with a live leader remains an operational error.
        """
        try:
            os.killpg(proc.pid, sig)
        except ProcessLookupError:
            pass
        except PermissionError:
            if proc.poll() is None:
                raise

    def stop_children(self):
        """Terminate all owned process groups, escalating after three seconds."""
        with self.mutex:
            processes = list(self.children)
        for sig in (signal.SIGTERM, signal.SIGKILL):
            for proc in processes:
                self.signal_group(proc, sig)
            if sig == signal.SIGTERM and processes:
                time.sleep(3)
        for proc in processes:
            proc.wait()

    def command(self, argv, log=None, cwd=None, timeout=None, fail_fast=False):
        """Run argv in an owned group; return stdout or raise with log diagnostics.

        log optionally redirects combined output; cwd defaults to the checkout.
        timeout bounds operational commands. Every exit cleans the process group.
        """
        output = log.open('w') if log else subprocess.PIPE
        proc = None
        try:
            with self.mutex:
                if self.stopping.is_set():
                    raise Failure('Job cancelled')
                env = dict(os.environ, GIT_TERMINAL_PROMPT='0',
                           GIT_SSH_COMMAND='ssh -oBatchMode=yes -oConnectTimeout=20')
                proc = subprocess.Popen(argv, cwd=cwd or self.repo, stdout=output,
                                        stderr=subprocess.STDOUT, text=True,
                                        start_new_session=True, env=env)
                self.children.add(proc)
                if log and log.name == 'perf.log':
                    self.perf_workers.add(proc)
            started = time.monotonic()
            while True:
                try:
                    result, _ = proc.communicate(timeout=min(10, timeout) if timeout else 10)
                    break
                except subprocess.TimeoutExpired:
                    self.resource_check()
                    if fail_fast and log:
                        with log.open('rb') as stream:
                            stream.seek(max(0, log.stat().st_size - 16384))
                            tail = stream.read().decode(errors='replace')
                        if any(line.startswith('FAIL   ') for line in tail.splitlines()):
                            raise CommandFailure('Stress group failed | ' + str(log))
                    if timeout and time.monotonic() - started >= timeout:
                        raise CommandFailure('Command timed out: ' + shlex.join(map(str, argv)))
                    if self.stopping.is_set():
                        raise Failure('Job cancelled')
                    self.say('RUN ' + Path(str(argv[0])).name + ' | logs ' + str(log))
                    if log and log.name in ('perf.log', 'runner.log'):
                        workers = len(self.perf_workers)
                        if log.name == 'runner.log':
                            lines = log.read_text(errors='replace').splitlines()
                            workers = sum(line.startswith('START  ') for line in lines) - sum(
                                line.startswith(('PASS   ', 'FAIL   ')) for line in lines)
                        if workers > 0:
                            self.sound('work', workers=workers)
            if proc.returncode:
                detail = str(log) if log else (result or '')[-2000:]
                raise CommandFailure('{} exited {} | {}'.format(shlex.join(map(str, argv)), proc.returncode, detail))
            return result or ''
        finally:
            if proc:
                self.signal_group(proc, signal.SIGKILL)
                proc.wait()
                if proc.stdout:
                    proc.stdout.close()
                with self.mutex:
                    self.children.discard(proc)
                    self.perf_workers.discard(proc)
            if log:
                output.close()

    def git(self, *args):
        """Run a bounded noninteractive Git command and return stripped stdout."""
        return self.command(['git', *args], timeout=120).strip()

    def validate(self):
        """Validate platform, tools, checkout, lock and CPU instruction capabilities."""
        for tool in ('git', 'make', os.environ.get('CXX', 'c++'), 'bash'):
            if not shutil.which(tool):
                raise Failure('Required executable missing: ' + tool)
            self.say(self.command([tool, '--version'], timeout=20).splitlines()[0])
        if not self.arm and not (platform.system() == 'Linux' and platform.machine() in ('x86_64', 'amd64')):
            raise Failure('Unsupported platform: ' + platform.platform())
        # Git common directory makes the lock shared even across linked worktrees.
        common = Path(self.git('rev-parse', '--git-common-dir'))
        if not common.is_absolute():
            common = self.repo / common
        self.lock_file = (common / 'tagent.lock').open('a')
        try:
            fcntl.flock(self.lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise Failure('Another TAgent owns this repository') from error
        if self.git('status', '--porcelain', '--untracked-files=no'):
            raise Failure('Tracked local changes present; use a clean checkout')
        if self.git('branch', '--show-current') != self.args.branch:
            raise Failure('Checkout must be on configured branch ' + self.args.branch)
        self.git('check-ref-format', 'refs/heads/' + self.args.branch)
        if not self.arm:
            self.flags = set(Path('/proc/cpuinfo').read_text().lower().replace('.', '_').split())
        self.resource_check()
        self.say('Validated Python {} | {} workers | {}'.format(platform.python_version(), self.jobs, self.repo))

    def sound(self, event, force=False, workers=1):
        """Schedule local themed audio without delaying jobs; workers sets cadence.

        Steam trains pass every 300 seconds within the active test cycle; diesel
        pulses retain their 15-second interval. force bypasses timing for demos.
        Important events finish the previous clip; audio failures remain nonfatal.
        """
        if self.args.no_sound or self.args.sound_volume == 0:
            return
        with self.audio_mutex:
            now = time.monotonic()
            if event == 'work' and not force:
                due = self.next_train_at if self.args.sound_theme == 'steampunk' else self.last_sound + 15
                if now < due or (self.audio_thread and self.audio_thread.is_alive()):
                    return
            if self.audio_thread and self.audio_thread.is_alive():
                self.audio_thread.join(timeout=55)
                if self.audio_thread.is_alive():
                    return
            self.last_sound = time.monotonic()
            if event == 'work' and self.args.sound_theme == 'steampunk':
                self.next_train_at = self.last_sound + 300
            self.audio_thread = threading.Thread(target=self.play_sound,
                                                args=(event, workers), daemon=True)
            self.audio_thread.start()

    def start_test_sound(self):
        """Announce the whole perf/stress test cycle and start its train clock."""
        self.next_train_at = time.monotonic() + 300
        if self.args.sound_theme == 'steampunk':
            self.sound('depart')

    def play_sound(self, event, workers):
        """Render and play one event through local audio, then console fallback.

        Called by the audio thread; event and workers are captured on scheduling.
        All playback processes belong to the agent's normal cleanup registry.
        """
        try:
            suffix = '-' + str(workers) if event == 'work' else ''
            path = self.home / ('sound-' + self.args.sound_theme + '-' + event + suffix + '.wav')
            duration = soundtrack.render(path, self.args.sound_theme, event,
                              self.args.sound_volume, workers, self.jobs)
            for player in ('afplay', 'pw-play', 'paplay', 'aplay'):
                if shutil.which(player):
                    try:
                        self.command([player, str(path)], timeout=duration + 2)
                        return
                    except (OSError, Failure):
                        pass
            with open('/dev/console', 'w') as console:
                console.write('\a' * (3 if event == 'fail' else 1))
        except (OSError, ValueError, subprocess.SubprocessError):
            if not self.audio_warned:
                self.audio_warned = True
                self.say('Local audio unavailable; continuing. Generated WAV files: ' + str(self.home))

    def pause(self, seconds, label):
        """Wait interruptibly with resource checks and visible fetch heartbeat."""
        deadline = time.monotonic() + seconds
        while time.monotonic() < deadline:
            self.resource_check()
            line = '{} | fetch {} | attempts {} | commit {} | next {:.0f}s | free {:.1f} GB'.format(
                label, self.last_fetch, self.attempts, self.commit[:12],
                deadline - time.monotonic(), shutil.disk_usage(self.home).free / 1024**3)
            if sys.stdout.isatty():
                print('\r\033[K' + line, end='', flush=True)
            else:
                self.say(line)
            if self.stopping.wait(min(15, max(0, deadline - time.monotonic()))):
                raise Failure('Stopped')

    def fetch(self, failed_commit=None):
        """Fetch with bounded recovery; optionally wait for a fix to failed_commit."""
        deadline = time.monotonic() + self.args.retry_minutes * 60
        delay = 60
        while True:
            self.attempts += 1
            try:
                self.git('fetch', '--no-tags', self.args.remote,
                         '+refs/heads/{}:refs/remotes/{}/{}'.format(
                             self.args.branch, self.args.remote, self.args.branch))
                self.last_fetch = time.strftime('%H:%M:%S')
                tip = self.git('rev-parse', self.args.remote + '/' + self.args.branch)
                if tip != failed_commit:
                    return tip
                self.say('Waiting for a fix to ' + failed_commit[:12])
            except CommandFailure as error:
                self.say('Fetch failed: ' + str(error))
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise Failure('Recovery window exhausted; operator intervention required')
            self.sound('fail')
            self.pause(min(delay, remaining), 'RECOVERY')
            delay = min(600, delay * 2)

    def supported(self, name):
        """Return whether this build suffix can execute on the current host."""
        if name.endswith('_avx512'):
            return {'avx512f', 'avx512bw', 'avx512dq', 'avx512vl'} <= self.flags
        for suffix, feature in [('avx2', 'avx2'),
                                ('sse42', 'sse4_2'), ('sse2', 'sse2')]:
            if name.endswith('_' + suffix):
                return feature in self.flags
        return True

    def builds(self, family):
        """Return explicit expected build names for perf or stress."""
        modes = ['debug', 'release']
        names = [family + '_' + mode for mode in modes]
        if self.arm:
            names += [family + '_' + mode + '_neon' for mode in modes]
        else:
            if family == 'stress':
                names += ['stress_release_gcc']
            for suffix in ['sse2', 'sse42', 'avx2'] + (['avx512'] if self.args.avx512 else []):
                for mode in modes if family == 'stress' or self.args.pedantic else ['release']:
                    names.append(family + '_' + mode + '_' + suffix)
            if family == 'perf':
                names += ['perf_release_64']
                if self.args.pedantic:
                    names += ['perf_debug_64']
        return names

    def cycle(self):
        """Build exact current commit; return False on test failure, True on success."""
        run = self.home / (datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f') + '-' + self.commit[:12])
        run.mkdir()
        manifest = {'commit': self.commit, 'platform': platform.platform(),
                    'jobs': self.jobs, 'pedantic': self.args.pedantic, 'builds': {}}
        manifest['builds'] = {family: self.builds(family) for family in ('perf', 'stress')}
        (run / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        for family in ('perf', 'stress'):
            directory = self.repo / 'tests' / family
            names = self.builds(family)
            manifest['builds'][family] = names
            for name in names:
                (directory / name).unlink(missing_ok=True)
            argv = ['bash', str(directory / 'build_all.sh')]
            if self.args.avx512:
                argv.append('--avx512')
            if family == 'perf' and self.args.pedantic:
                argv.append('--pedantic')
            self.sound('building')
            self.say('BUILD ' + family + ' | ' + str(run))
            self.command(argv, run / ('build-' + family + '.log'))
            for name in names:
                if not os.access(directory / name, os.X_OK):
                    raise Failure('Expected build missing: ' + str(directory / name))
        (run / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        self.sound('build')
        self.start_test_sound()
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as pool:
                futures = []
                for name in manifest['builds']['perf']:
                    if (not self.args.pedantic and '_debug' in name) or not self.supported(name):
                        self.say('SKIP perf ' + name)
                        continue
                    work = run / name
                    work.mkdir()
                    futures.append(pool.submit(self.command,
                        [str(self.repo / 'tests/perf' / name)], work / 'perf.log', work))
                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()
                    except BaseException:
                        self.stopping.set()
                        for pending in futures:
                            pending.cancel()
                        self.stop_children()
                        raise
            for name in manifest['builds']['stress']:
                if not self.supported(name):
                    self.say('SKIP unsupported ' + name)
                    continue
                work = run / name
                work.mkdir()
                # Runner selects groups; wrapper gives each invocation its own cwd.
                runner = work / 'run_all_pl.sh'
                shutil.copy2(self.repo / 'tests/stress/run_all_pl.sh', runner)
                wrapper = work / 'worker.sh'
                wrapper.write_text('#!/bin/bash\nset -eu\ncd "$(mktemp -d ' +
                    shlex.quote(str(work / 'job.XXXXXX')) + ')"\nexec ' +
                    shlex.quote(str(self.repo / 'tests/stress' / name)) + ' "$@"\n')
                wrapper.chmod(0o755)
                argv = ['bash', str(runner), '--binary', './worker.sh',
                        '--jobs', str(self.jobs), '--fail-fast']
                if self.args.pedantic:
                    argv.append('--pedantic')
                self.say('STRESS ' + name + ' | ' + str(work))
                self.command(argv, work / 'runner.log', work, fail_fast=True)
            return True
        except CommandFailure as error:
            self.say('TEST FAILURE ' + str(error))
            self.stopping.clear()
            self.sound('fail')
            return False
        finally:
            self.next_train_at = math.inf
            self.stopping.clear()

    def run(self):
        """Run update/build/test/wait loop; return exit status, or exec on self-update."""
        if self.args.sound_demo:
            self.say('SOUND DEMO | ' + self.args.sound_theme +
                     ' | volume ' + str(self.args.sound_volume))
            sequence = [('start', 1), ('building', 1), ('build', 1),
                        ('work', 1), ('work', self.jobs),
                        ('success', 1), ('idle', 1), ('fail', 1)]
            if self.args.sound_theme == 'steampunk':
                sequence.insert(3, ('depart', 1))
            clips = []
            for event, workers in sequence:
                self.say('DEMO ' + event + (' | workers ' + str(workers) if event == 'work' else ''))
                self.sound(event, force=True, workers=workers)
                if self.audio_thread:
                    self.audio_thread.join(timeout=55)
                suffix = '-' + str(workers) if event == 'work' else ''
                clip = self.home / ('sound-' + self.args.sound_theme + '-' + event + suffix + '.wav')
                if not self.args.no_sound and self.args.sound_volume > 0 and clip.exists():
                    with wave.open(str(clip), 'rb') as stream:
                        clips.append(stream.readframes(stream.getnframes()))
            if clips:
                with wave.open(str(self.home / ('demo-' + self.args.sound_theme + '.wav')), 'wb') as stream:
                    stream.setparams((2, 2, soundtrack.RATE, 0, 'NONE', 'not compressed'))
                    stream.writeframes((b'\0' * (soundtrack.RATE // 3 * 4)).join(clips))
            self.say('Sound demo complete | WAV files: ' + str(self.home))
            return 0
        if self.args.sound_preview:
            self.sound(self.args.sound_preview, force=True, workers=self.jobs)
            return 0
        self.validate()
        if self.args.check:
            return 0
        self.sound('start')
        failed = None
        state_path = self.home / 'success.json'
        while True:
            tip = self.fetch(failed)
            self.commit = self.git('rev-parse', 'HEAD')
            if tip != self.commit:
                if self.git('status', '--porcelain', '--untracked-files=no'):
                    raise Failure('Tracked modifications prevent update')
                agent_changed = bool(self.git('diff', '--name-only', self.commit, tip, '--', 'tests/tagent'))
                self.git('merge', '--ff-only', tip)
                self.commit = tip
                if agent_changed:
                    self.say('TAgent updated; restarting launcher. If restart fails, rerun run_tagent.sh.')
                    self.lock_file.close()
                    os.execv('/bin/bash', ['bash', str(self.repo / 'tests/tagent/run_tagent.sh'), *sys.argv[1:]])
            key = {'commit': self.commit, 'pedantic': self.args.pedantic, 'avx512': self.args.avx512}
            previous = json.loads(state_path.read_text()) if state_path.exists() else None
            if self.args.force or previous != key or failed:
                if not self.cycle():
                    if self.args.once:
                        return 3
                    failed = self.commit
                    continue
                temp = state_path.with_suffix('.tmp')
                temp.write_text(json.dumps(key))
                temp.replace(state_path)
                self.sound('success')
                self.say('SUCCESS ' + self.commit)
                self.args.force = False
                failed = None
            if self.args.once:
                return 0
            self.sound('idle')
            self.pause(self.args.poll_minutes * 60, 'WAIT')


def main():
    """Start agent, handle signals, and guarantee child cleanup and lock release."""
    agent = Agent(arguments())
    def interrupted(signum, frame):
        """Convert termination signals to an interrupt handled by the owner."""
        agent.stopping.set()
        agent.stop_children()
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    try:
        return agent.run()
    except KeyboardInterrupt:
        agent.say('Interrupted; stopping all child jobs')
        return 130
    except (Failure, OSError, ValueError) as error:
        agent.say('ERROR ' + str(error) + ' | logs ' + str(agent.home))
        agent.sound('fail')
        return 3 if 'Recovery window exhausted' in str(error) else 1
    finally:
        if agent.audio_thread and not agent.stopping.is_set():
            agent.audio_thread.join(timeout=55)
        agent.stopping.set()
        agent.stop_children()
        if agent.lock_file:
            agent.lock_file.close()


if __name__ == '__main__':
    sys.exit(main())
