#!/usr/bin/env python3
"""Remote BitMagic supervisor; standard library only, Python 3.9+ and POSIX.

Owns subprocess groups and one checkout lock. Never discards repository changes.
Build/test outputs live in persistent per-cycle directories for investigation.
"""
import argparse
import concurrent.futures
import datetime
import depot
import fcntl
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import soundtrack
import wave


class Failure(RuntimeError):
    """An operational failure, or a failed build/test with its recorded command."""


class CommandFailure(Failure):
    """A child command returned nonzero; test phases may wait for a fix."""


class WarningCounter:
    """Incrementally count complete GCC/Clang-style warning diagnostic lines."""

    def __init__(self):
        """Start an empty counter for a single build cycle."""
        self.count = 0
        self.errors = 0
        self.offsets = {}
        self.partial = {}

    def count_line(self, line):
        """Count a single complete diagnostic line after removing ANSI colors."""
        line = re.sub(rb'\x1b\[[0-9;]*m', b'', line)
        if re.search(rb'(?:^|:\s*)warning:', line, re.IGNORECASE):
            self.count += 1
        if re.search(rb'(?:^|:\s*)(?:fatal\s+)?error:', line, re.IGNORECASE):
            self.errors += 1

    def observe(self, path, final=False):
        """Read appended bytes from path and return the total across build logs.

        Count case-insensitive warning:/error:/fatal error: diagnostic lines.
        Summary lines and incomplete lines are ignored unless final=True, which
        includes the last unterminated line. ANSI colors are stripped.
        """
        with path.open('rb') as stream:
            stream.seek(self.offsets.get(path, 0))
            while True:
                chunk = stream.read(65536)
                if not chunk:
                    break
                lines = (self.partial.get(path, b'') + chunk).split(b'\n')
                self.partial[path] = lines.pop()
                for line in lines:
                    self.count_line(line)
            self.offsets[path] = stream.tell()
            if final and self.partial.get(path):
                self.count_line(self.partial.pop(path))
        return self.count


def build_report_text(count, progress=None):
    """Return a railway-themed status, accurately scoped to warnings seen so far."""
    status = ''
    if progress is not None:
        perf_ready, perf_pending, stress_ready, stress_pending = progress
        status = (' Performance builds: {} ready, {} pending. '
                  'Stress builds: {} ready, {} pending.'.format(
                      perf_ready, perf_pending, stress_ready, stress_pending))
    if count == 0:
        weather = random.choice((
            'Tracks are clear, no snowdrifts reported.',
            'The mountain pass is clear of snow and avalanches.',
            'Clear signals ahead and fair weather on the line.',
            'The snowplough reports clear tracks all the way to the junction.',
            'No ice on the rails and no delays at the crossing.',
            'The line is open, and the points are set for a clear run.',
        ))
        return 'Dispatch to BitMagic. Build in progress. No compiler warnings so far.' + status + ' ' + weather
    word = 'warning' if count == 1 else 'warnings'
    return ('Dispatch to BitMagic. Caution on the line. {} compiler {} so far.'.format(count, word)
            + status + ' Please check the build log.')


def arguments():
    """Return validated CLI configuration, defaulting to launcher environment."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name, default in [('poll-minutes', 5), ('retry-minutes', 60),
                          ('reserve-cpus', 2), ('jobs', 0), ('min-free-gb', 5),
                          ('log-retention-days', 7), ('build-report-seconds', 60)]:
        parser.add_argument('--' + name, type=float, default=float(
            os.environ.get(name.upper().replace('-', '_'), default)))
    for name, default in [('remote', 'origin'), ('branch', 'master')]:
        parser.add_argument('--' + name, default=os.environ.get(name.upper(), default))
    parser.add_argument('--state-dir', type=Path)
    parser.add_argument('--startup', choices=('build', 'changed'),
                        default=os.environ.get('STARTUP', 'build'),
                        help='build on every manual start (default), or only when changed')
    parser.add_argument('--resume-cycle', action='store_true', help=argparse.SUPPRESS)
    for name in ['pedantic', 'avx512', 'once', 'check', 'force']:
        parser.add_argument('--' + name, action='store_true')
    parser.add_argument('--no-sound', action='store_true',
                        default=os.environ.get('SOUND', 'yes') == 'no')
    parser.add_argument('--no-depot', action='store_true',
                        default=os.environ.get('DEPOT_SOUND', 'yes') == 'no',
                        help='disable steampunk depot and train ambience')
    parser.add_argument('--sound-theme', choices=soundtrack.THEMES,
                        default=os.environ.get('SOUND_THEME', 'steampunk'))
    parser.add_argument('--sound-volume', type=float,
                        default=float(os.environ.get('SOUND_VOLUME', '0.35')))
    parser.add_argument('--tts', choices=('auto', 'off'),
                        default=os.environ.get('TTS_MODE', 'auto'),
                        help='use local eSpeak NG when available, or use recorded voice')
    parser.add_argument('--tts-voice', default=os.environ.get('TTS_VOICE', 'en-us+m3'))
    parser.add_argument('--dispatcher-voice', default=os.environ.get('DISPATCHER_VOICE', 'en-us+f3,en-us+f4,en-us+m3'),
                        help='comma-separated eSpeak NG voices for build reports')
    parser.add_argument('--sound-demo', action='store_true',
                        help='Play a complete theme demonstration and exit without Git/builds')
    parser.add_argument('--sound-preview', choices=soundtrack.EVENTS + ('depot', 'train-ambient', 'dispatch', 'depart-perf', 'depart-stress', 'build-clear', 'build-warning', 'build-final-ok', 'build-final-failed'),
                        help='Play one event locally and exit without Git/build validation')
    args = parser.parse_args()
    if args.tts not in ('auto', 'off'):
        parser.error('TTS_MODE must be auto or off')
    if args.startup not in ('build', 'changed'):
        parser.error('STARTUP must be build or changed')
    if args.sound_theme not in soundtrack.THEMES:
        parser.error('Unknown SOUND_THEME')
    if args.sound_preview in ('depot', 'train-ambient', 'dispatch', 'depart-perf', 'depart-stress', 'build-clear', 'build-warning', 'build-final-ok', 'build-final-failed') and args.sound_theme != 'steampunk':
        parser.error('voice preview requires --sound-theme steampunk')
    if not math.isfinite(args.sound_volume) or not 0 <= args.sound_volume <= 1:
        parser.error('sound-volume must be between 0 and 1')
    for name in ['poll_minutes', 'retry_minutes', 'min_free_gb', 'log_retention_days', 'build_report_seconds']:
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
        self.home = (args.state_dir or self.repo / 'TAgent-logs').resolve()
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
        self.audio_active = None
        self.audio_epoch = 0
        self.audio_mutex = threading.Lock()
        self.depot = None
        self.audio_warned = False
        self.audio_available = True
        self.console_notice = False
        self.tts_warned = False
        self.build_warnings = WarningCounter()
        self.build_lock = threading.Lock()
        self.build_logs = []
        self.build_targets = {}
        self.build_report_interval = args.build_report_seconds
        self.next_build_report = math.inf
        self.tts_engine = (os.environ.get('TAGENT_TTS_ENGINE') if
                           os.environ.get('TAGENT_TTS_CHECKED') == '1' else
                           shutil.which('espeak-ng'))
        self.perf_workers = set()
        self.perf_skipped = False
        self.lock_file = None
        self.poll_lock = threading.Lock()
        self.cleanup_lock = threading.Lock()
        self.next_remote_poll = math.inf
        self.pending_tip = None
        self.poll_failure_since = None
        self.poll_retry_delay = 60
        self.active_run = None
        self.stress_stop_file = None
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

    def cleanup_logs(self):
        """Remove expired owned cycle folders; never follow symlinks or remove an active run."""
        if not self.cleanup_lock.acquire(blocking=False):
            return
        try:
            cutoff = time.time() - self.args.log_retention_days * 86400
            for path in self.home.iterdir():
                if path.is_symlink() or not path.is_dir() or path == self.active_run:
                    continue
                if not re.fullmatch(r'\d{8}-\d{6}-\d{6}-[0-9a-f]{6,12}', path.name):
                    continue
                try:
                    manifest = json.loads((path / 'manifest.json').read_text())
                    if manifest.get('tagent_run') != 1:
                        continue
                    created = datetime.datetime.strptime(path.name[:22], '%Y%m%d-%H%M%S-%f').timestamp()
                    if created < cutoff:
                        shutil.rmtree(path)
                        self.say('Removed expired cycle logs: ' + str(path))
                except (OSError, ValueError, TypeError) as error:
                    self.say('Log cleanup skipped {}: {}'.format(path, error))
        finally:
            self.cleanup_lock.release()

    def poll_remote(self):
        """Check remote without changing checkout; enforce bounded network recovery."""
        if time.monotonic() < self.next_remote_poll or not self.poll_lock.acquire(False):
            return
        try:
            now = time.monotonic()
            if now < self.next_remote_poll:
                return
            self.cleanup_logs()
            self.attempts += 1
            try:
                output = self.git('ls-remote', '--exit-code', self.args.remote,
                                  'refs/heads/' + self.args.branch)
                lines = output.splitlines()
                if len(lines) != 1 or len(lines[0].split()) != 2:
                    raise CommandFailure('Unexpected remote branch response')
                tip, ref = lines[0].split()
                if ref != 'refs/heads/' + self.args.branch or not re.fullmatch(r'[0-9a-fA-F]{40,64}', tip):
                    raise CommandFailure('Invalid remote branch response')
                self.last_fetch = time.strftime('%H:%M:%S') + ' (remote check)'
                self.poll_failure_since = None
                self.poll_retry_delay = 60
                self.next_remote_poll = time.monotonic() + self.args.poll_minutes * 60
                if tip != self.commit:
                    self.pending_tip = tip
                    if self.stress_stop_file is not None:
                        self.stress_stop_file.touch()
                    self.say('New commit ' + tip[:12] + '; waiting for current tests to finish')
                else:
                    self.pending_tip = None
            except CommandFailure as error:
                now = time.monotonic()
                if self.poll_failure_since is None:
                    self.poll_failure_since = now
                    self.sound('fail')
                self.say('Remote check failed: ' + str(error))
                remaining = self.args.retry_minutes * 60 - (now - self.poll_failure_since)
                if remaining <= 0:
                    raise Failure('Recovery window exhausted during remote checks')
                self.next_remote_poll = now + min(self.poll_retry_delay, remaining)
                self.poll_retry_delay = min(600, self.poll_retry_delay * 2)
        finally:
            self.poll_lock.release()

    def update_waiting(self):
        """Check due remote updates at a safe test boundary; return whether to rebuild."""
        self.poll_remote()
        # Do not declare success or start another build while the latest remote
        # check is in recovery. Continue the same bounded backoff window.
        while self.poll_failure_since is not None:
            self.pause(max(0, self.next_remote_poll - time.monotonic()), 'RECOVERY')
            self.poll_remote()
        if self.pending_tip:
            self.say('Current tests finished; fetching update and restarting cycle')
            return True
        return False

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

    def command(self, argv, log=None, cwd=None, timeout=None, fail_fast=False, allowed_exit_codes=()):
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
                if log and log.name == 'perf.log' and self.pending_tip:
                    self.perf_skipped = True
                    return 'skipped-update'
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
                    if log and log.name in ('build-perf.log', 'build-stress.log'):
                        self.build_progress(log)
                    else:
                        self.say('RUN ' + Path(str(argv[0])).name + ' | logs ' + str(log))
                    if log and log.name in ('perf.log', 'runner.log'):
                        self.poll_remote()
                        workers = len(self.perf_workers)
                        if log.name == 'runner.log':
                            lines = log.read_text(errors='replace').splitlines()
                            workers = sum(line.startswith('START  ') for line in lines) - sum(
                                line.startswith(('PASS   ', 'FAIL   ')) for line in lines)
                        if workers > 0:
                            self.sound('work', workers=workers)
            if proc.returncode in allowed_exit_codes:
                return str(proc.returncode)
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

    def report_audio(self):
        """Report optional TTS/audio availability without requiring either resource."""
        if self.args.no_sound or self.args.sound_volume == 0:
            self.say('Sound disabled')
            return
        players = [p for p in ('afplay', 'pw-play', 'paplay', 'aplay') if shutil.which(p)]
        if not players:
            if os.access('/dev/console', os.W_OK):
                self.say('No audio player found; will try the local console bell')
            else:
                self.audio_available = False
                self.say('No audio player found and console unavailable; continuing without sound')
        if self.args.sound_theme == 'steampunk' and self.args.tts == 'auto':
            if self.tts_engine:
                self.say('TTS: eSpeak NG available (' + self.tts_engine + ')')
            else:
                self.tts_warned = True
                self.say('TTS: eSpeak NG unavailable; using prerecorded voice')

    def wait_for_audio(self, timeout=300):
        """Wait in short interruptible steps for a bounded in-flight audio clip.

        Return whether the audio thread has finished. This includes synthesis
        time as well as playback, which matters for the longer station departure.
        """
        deadline = time.monotonic() + timeout
        thread = self.audio_thread
        while thread and thread.is_alive() and not self.stopping.is_set():
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            thread.join(timeout=min(1, remaining))
        return not thread or not thread.is_alive()

    def sound(self, event, force=False, workers=1):
        """Schedule local themed audio without delaying jobs; workers sets cadence.

        Steam trains pass every 300 seconds within the active test cycle.
        force bypasses timing for demos.
        Events queue behind the current clip without blocking tests. A failure
        drops queued chatter and sounds after the currently playing clip.
        """
        if self.args.no_sound or self.args.sound_volume == 0:
            return
        if not self.audio_available and not (self.args.sound_demo or self.args.sound_preview):
            return
        with self.audio_mutex:
            now = time.monotonic()
            if event == 'work' and not force:
                due = self.next_train_at
                if now < due or (self.audio_thread and self.audio_thread.is_alive()):
                    return
            if event == 'fail':
                # Discard queued non-emergency chatter; finish the one clip
                # already playing before the distress signal.
                self.audio_epoch += 1
                previous = self.audio_active if self.audio_active and self.audio_active.is_alive() else None
            else:
                previous = self.audio_thread if self.audio_thread and self.audio_thread.is_alive() else None
            epoch = self.audio_epoch
            self.last_sound = time.monotonic()
            if event == 'work' and self.args.sound_theme == 'steampunk':
                self.next_train_at = self.last_sound + 300
            selected = event
            if event == 'work' and not force and self.args.sound_theme == 'steampunk':
                selected = random.choice(('work', 'dispatch'))
                self.say('STEAMPUNK SOUND ' + ('passing train' if selected == 'work' else 'dispatch'))
            self.audio_thread = threading.Thread(target=self.play_after,
                                                args=(previous, selected, workers, epoch), daemon=True)
            self.audio_thread.start()

    def play_after(self, previous, event, workers, epoch):
        """Play queued audio after its predecessor without blocking test work.

        epoch invalidates queued clips when a failure signal takes precedence.
        Already playing audio finishes; pending chatter cannot delay the alarm.
        """
        while previous and previous.is_alive() and not self.stopping.is_set():
            if epoch != self.audio_epoch:
                return
            previous.join(timeout=.25)
        if self.stopping.is_set() or epoch != self.audio_epoch:
            return
        with self.audio_mutex:
            if epoch != self.audio_epoch:
                return
            self.audio_active = threading.current_thread()
        try:
            if not self.audio_available and not (self.args.sound_demo or self.args.sound_preview):
                return
            if self.depot:
                self.depot.mute()
            if self.stopping.is_set() or epoch != self.audio_epoch:
                return
            self.play_sound(event, workers)
        finally:
            if self.depot:
                self.depot.unmute()
            with self.audio_mutex:
                if self.audio_active is threading.current_thread():
                    self.audio_active = None

    def start_test_sound(self, phase='perf'):
        """Announce a test phase without resetting the shared five-minute clock."""
        if self.depot:
            self.depot.set_mode('train')
        if self.next_train_at == math.inf:
            self.next_train_at = time.monotonic() + 300
        if self.args.sound_theme == 'steampunk':
            self.sound('depart-' + phase)

    def build_progress(self, log):
        """Report aggregate warnings and completed target counts on a growing interval."""
        with self.build_lock:
            for path in self.build_logs or [log]:
                if path.exists():
                    self.build_warnings.observe(path)
            count = self.build_warnings.count
            now = time.monotonic()
            if now >= self.next_build_report:
                self.build_report_interval += 10
                self.next_build_report = now + self.build_report_interval
                progress = []
                for family in ('perf', 'stress'):
                    targets = self.build_targets.get(family, [])
                    ready = sum(path.is_file() and os.access(path, os.X_OK) for path in targets)
                    progress.extend((ready, len(targets) - ready))
                self.say('BUILD STATUS: {} warning(s) | perf: {} ready, {} pending | '
                         'stress: {} ready, {} pending'.format(count, *progress))
                if self.args.sound_theme == 'steampunk':
                    self.sound('build-clear' if count == 0 else 'build-warning',
                               workers=(count, *progress))

    def dispatcher_voice(self):
        """Choose a configured eSpeak voice for this announcement."""
        voices = [voice.strip() for voice in self.args.dispatcher_voice.split(',') if voice.strip()]
        return random.choice(voices or ['en-us+f3'])

    def render_build_report(self, path, report):
        """Synthesize a status snapshot or use the prerecorded warning fallback."""
        if isinstance(report, tuple):
            count, *progress = report
        else:
            count, progress = report, None
        if self.args.tts == 'auto' and self.tts_engine:
            try:
                with tempfile.TemporaryDirectory(prefix='tts-', dir=self.home) as temp:
                    speech = Path(temp) / 'speech.wav'
                    self.command([self.tts_engine, '-v', self.dispatcher_voice(),
                                  '-s', '145', '-p', '57' if count else '50',
                                  '-w', str(speech), build_report_text(count, progress)], timeout=15)
                    return soundtrack.mix_dynamic_dispatch(speech, path, self.args.sound_volume,
                                                           1, self.jobs, background='build')
            except (OSError, ValueError, EOFError, wave.Error, Failure) as error:
                if not self.tts_warned:
                    self.tts_warned = True
                    self.say('TTS unavailable or failed ({}); using prerecorded voice'.format(error))
        asset = Path(soundtrack.__file__).resolve().parent / 'assets' / (
            'build-warning.wav' if count else 'build-clear.wav')
        return soundtrack.mix_dynamic_dispatch(asset, path, self.args.sound_volume,
                                               1, self.jobs, background='build')

    def render_build_summary(self, path, result):
        """Announce final success and exact diagnostic counts, with offline fallback."""
        errors, warnings, success = result
        if self.args.tts == 'auto' and self.tts_engine:
            try:
                with tempfile.TemporaryDirectory(prefix='tts-', dir=self.home) as temp:
                    speech = Path(temp) / 'speech.wav'
                    text = ('Dispatch to BitMagic. Build {}. {} compiler errors. '
                            '{} compiler warnings.'.format('complete' if success else 'failed',
                                                           errors, warnings))
                    self.command([self.tts_engine, '-v', self.dispatcher_voice(),
                                  '-s', '145', '-p', '50' if success else '57',
                                  '-w', str(speech), text], timeout=15)
                    return soundtrack.mix_dynamic_dispatch(speech, path, self.args.sound_volume,
                                                           1, self.jobs, background='build')
            except (OSError, ValueError, EOFError, wave.Error, Failure) as error:
                if not self.tts_warned:
                    self.tts_warned = True
                    self.say('TTS unavailable or failed ({}); using prerecorded voice'.format(error))
        return soundtrack.render_build_summary_fallback(path, errors, warnings, success,
                                                        self.args.sound_volume)

    def finish_build(self, success):
        """Print and announce final recognized error/warning line counts."""
        counts = (self.build_warnings.errors, self.build_warnings.count, success)
        self.say('BUILD {} | diagnostic errors: {} | warnings: {}'.format(
            'COMPLETE' if success else 'FAILED', counts[0], counts[1]))
        if self.args.sound_theme == 'steampunk':
            self.sound('build-final-ok' if success else 'build-final-failed', workers=counts)

    def render_dispatch(self, path, workers):
        """Create a live crew-count announcement, falling back to the bundled WAV.

        path is the output file; workers is the last observed active test-process
        count. A missing engine, bad voice, timeout or invalid WAV never stops tests.
        """
        if self.args.tts == 'auto' and not self.tts_engine and not self.tts_warned:
            self.tts_warned = True
            self.say('TTS: eSpeak NG unavailable; using prerecorded voice')
        if self.args.tts == 'auto' and self.tts_engine:
            try:
                with tempfile.TemporaryDirectory(prefix='tts-', dir=self.home) as temp:
                    speech = Path(temp) / 'speech.wav'
                    count = max(1, int(workers))
                    crew = 'crew' if count == 1 else 'crews'
                    text = 'BitMagic freight rolling. {} test {} at work.'.format(count, crew)
                    self.command([self.tts_engine, '-v', self.args.tts_voice,
                                  '-s', '145', '-p', '40', '-w', str(speech), text], timeout=15)
                    return soundtrack.mix_dynamic_dispatch(speech, path, self.args.sound_volume,
                                                           count, self.jobs)
            except (OSError, ValueError, EOFError, wave.Error, Failure) as error:
                if not self.tts_warned:
                    self.tts_warned = True
                    self.say('TTS unavailable or failed ({}); using prerecorded voice'.format(error))
        return soundtrack.render_dispatch(path, self.args.sound_volume)

    def play_sound(self, event, workers):
        """Render and play one event through local audio, then console fallback.

        Called by the audio thread; event and workers are captured on scheduling.
        All playback processes belong to the agent's normal cleanup registry.
        """
        rendered = False
        try:
            suffix = '-' + str(workers) if event == 'work' else ''
            path = self.home / ('sound-' + self.args.sound_theme + '-' + event + suffix + '.wav')
            if event in ('depot', 'train-ambient'):
                duration = depot.render_preview(path, self.args.sound_volume,
                                                'depot' if event == 'depot' else 'train')
            elif event in ('depart-perf', 'depart-stress'):
                duration = soundtrack.render(path, self.args.sound_theme, 'depart',
                                             self.args.sound_volume, workers, self.jobs,
                                             announcement=event.split('-')[1])
            elif event in ('build-final-ok', 'build-final-failed'):
                duration = self.render_build_summary(path, workers)
            elif event in ('build-clear', 'build-warning'):
                duration = self.render_build_report(path, workers)
            elif event == 'dispatch':
                duration = self.render_dispatch(path, workers)
            else:
                duration = soundtrack.render(path, self.args.sound_theme, event,
                                             self.args.sound_volume, workers, self.jobs)
            rendered = True
            # Explicit demos/previews still export WAVs on a headless machine.
            if not self.audio_available:
                return
            for player in ('afplay', 'pw-play', 'paplay', 'aplay'):
                if shutil.which(player):
                    try:
                        self.command([player, str(path)], timeout=duration + 2)
                        return
                    except (OSError, Failure):
                        pass
            if not self.console_notice:
                self.console_notice = True
                self.say('No working audio player; trying the local console bell')
            with open('/dev/console', 'w') as console:
                console.write('\a' * (3 if event == 'fail' else 1))
        except (OSError, ValueError, EOFError, wave.Error, subprocess.SubprocessError) as error:
            if rendered:
                self.audio_available = False
            if not self.audio_warned:
                self.audio_warned = True
                if rendered:
                    self.say('Local audio unavailable; continuing without sound. Generated WAV files: ' + str(self.home))
                else:
                    self.say('Sound generation failed: {}; continuing tests'.format(error))

    def pause(self, seconds, label):
        """Wait interruptibly with resource checks and visible fetch heartbeat."""
        deadline = time.monotonic() + seconds
        while time.monotonic() < deadline:
            self.resource_check()
            self.cleanup_logs()
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
        if self.depot:
            self.depot.set_mode('depot')
        run = self.home / (datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f') + '-' + self.commit[:12])
        run.mkdir()
        self.active_run = run
        self.build_warnings = WarningCounter()
        self.build_report_interval = self.args.build_report_seconds
        self.next_build_report = time.monotonic() + self.args.build_report_seconds
        manifest = {'tagent_run': 1, 'commit': self.commit, 'platform': platform.platform(),
                    'jobs': self.jobs, 'pedantic': self.args.pedantic, 'builds': {}}
        manifest['builds'] = {family: self.builds(family) for family in ('perf', 'stress')}
        (run / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        self.build_logs = [run / ('build-' + family + '.log') for family in ('perf', 'stress')]
        build_tasks = []
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
            build_tasks.append((family, directory, names, argv, run / ('build-' + family + '.log')))
        self.build_targets = {family: [directory / name for name in names]
                              for family, directory, names, _, _ in build_tasks}
        self.sound('building')

        def run_build(task):
            """Run one build family and retain its separate diagnostic log."""
            family, _, _, argv, build_log = task
            self.say('BUILD ' + family + ' | ' + str(run))
            try:
                self.command(argv, build_log)
            finally:
                with self.build_lock:
                    if build_log.exists():
                        self.build_warnings.observe(build_log, final=True)

        failures = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(2, self.jobs)) as pool:
            futures = [pool.submit(run_build, task) for task in build_tasks]
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Failure as error:
                    failures.append(error)
        for path in self.build_logs:
            if path.exists():
                self.build_warnings.observe(path, final=True)
        if failures:
            self.finish_build(False)
            raise failures[0]
        for _, directory, names, _, _ in build_tasks:
            for name in names:
                if not os.access(directory / name, os.X_OK):
                    self.finish_build(False)
                    raise Failure('Expected build missing: ' + str(directory / name))
        (run / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        self.finish_build(True)
        self.sound('build')
        self.pending_tip = None
        self.perf_skipped = False
        self.poll_failure_since = None
        self.poll_retry_delay = 60
        self.next_remote_poll = time.monotonic() + self.args.poll_minutes * 60
        self.start_test_sound('perf')
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
            if self.update_waiting() or self.perf_skipped:
                self.args.force = True
                return 'updated'
            self.start_test_sound('stress')
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
                self.stress_stop_file = work / 'stop-for-update'
                argv = ['bash', str(runner), '--binary', './worker.sh',
                        '--jobs', str(self.jobs), '--fail-fast',
                        '--stop-file', str(self.stress_stop_file)]
                if self.args.pedantic:
                    argv.append('--pedantic')
                self.say('STRESS ' + name + ' | ' + str(work))
                status = self.command(argv, work / 'runner.log', work, fail_fast=True,
                                      allowed_exit_codes=(75,))
                self.stress_stop_file = None
                if status == '75':
                    # A branch may have moved back while the workers drained.
                    # Rebuild even then, since this stress suite is incomplete.
                    self.args.force = True
                    self.update_waiting()
                    return 'updated'
                if self.update_waiting():
                    return 'updated'
            return True
        except CommandFailure as error:
            self.say('TEST FAILURE ' + str(error))
            self.stopping.clear()
            self.sound('fail')
            return False
        finally:
            if self.depot:
                self.depot.set_mode('depot')
            self.active_run = None
            self.stress_stop_file = None
            self.next_remote_poll = math.inf
            self.next_train_at = math.inf
            self.stopping.clear()

    def run(self):
        """Run update/build/test/wait loop; return exit status, or exec on self-update."""
        self.report_audio()
        if self.args.sound_demo:
            self.say('SOUND DEMO | ' + self.args.sound_theme +
                     ' | volume ' + str(self.args.sound_volume))
            sequence = [('start', 1), ('building', 1), ('build', 1),
                        ('work', 1), ('work', self.jobs),
                        ('success', 1), ('idle', 1), ('fail', 1)]
            if self.args.sound_theme == 'steampunk':
                sequence.insert(3, ('depart-perf', 1))
                sequence.insert(4, ('depart-stress', 1))
                sequence.insert(7, ('dispatch', 1))
                sequence.extend([('build-clear', (0, 2, 4, 3, 7)),
                                 ('build-warning', (3, 5, 1, 7, 3)),
                                 ('build-final-ok', (0, 3, True)),
                                 ('build-final-failed', (2, 3, False))])
            clips = []
            for event, workers in sequence:
                self.say('DEMO ' + event + (' | workers ' + str(workers) if event == 'work' else ''))
                self.sound(event, force=True, workers=workers)
                if self.audio_thread:
                    self.wait_for_audio()
                suffix = '-' + str(workers) if event == 'work' else ''
                clip = self.home / ('sound-' + self.args.sound_theme + '-' + event + suffix + '.wav')
                if not self.args.no_sound and self.args.sound_volume > 0 and clip.exists():
                    with wave.open(str(clip), 'rb') as stream:
                        clips.append(stream.readframes(stream.getnframes()))
            if clips:
                with wave.open(str(self.home / ('demo-' + self.args.sound_theme + '.wav')), 'wb') as stream:
                    stream.setparams((1, 2, soundtrack.RATE, 0, 'NONE', 'not compressed'))
                    stream.writeframes((b'\0' * (soundtrack.RATE // 3 * 2)).join(clips))
            self.say('Sound demo complete | WAV files: ' + str(self.home))
            return 0
        if self.args.sound_preview:
            preview_count = ((0, 2, 4, 3, 7) if self.args.sound_preview == 'build-clear' else
                             (3, 5, 1, 7, 3) if self.args.sound_preview == 'build-warning' else
                             (0, 3, True) if self.args.sound_preview == 'build-final-ok' else
                             (2, 3, False) if self.args.sound_preview == 'build-final-failed' else self.jobs)
            self.sound(self.args.sound_preview, force=True, workers=preview_count)
            return 0
        self.validate()
        if self.args.check:
            return 0
        if (self.args.sound_theme == 'steampunk' and not self.args.no_sound and
                not self.args.no_depot and self.args.sound_volume > 0 and self.audio_available):
            self.depot = depot.Depot(self)
            self.depot.start()
        self.cleanup_logs()
        self.sound('start')
        failed = None
        state_path = self.home / 'success.json'
        startup_pending = True
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
                    restart_args = list(sys.argv[1:])
                    if '--resume-cycle' not in restart_args:
                        restart_args.append('--resume-cycle')
                    if self.depot:
                        self.depot.close()
                    self.stop_children()
                    os.execv('/bin/bash', ['bash', str(self.repo / 'tests/tagent/run_tagent.sh'), *restart_args])
            key = {'commit': self.commit, 'pedantic': self.args.pedantic, 'avx512': self.args.avx512}
            previous = json.loads(state_path.read_text()) if state_path.exists() else None
            start_cycle = startup_pending and (self.args.startup == 'build' or self.args.resume_cycle)
            startup_pending = False
            if start_cycle or self.args.force or previous != key or failed:
                result = self.cycle()
                if result == 'updated':
                    self.args.force = True
                    failed = None
                    continue
                if not result:
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
            agent.wait_for_audio()
        agent.stopping.set()
        if agent.depot:
            agent.depot.close()
        agent.stop_children()
        if agent.lock_file:
            agent.lock_file.close()


if __name__ == '__main__':
    sys.exit(main())
