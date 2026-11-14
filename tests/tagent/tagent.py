#!/usr/bin/env python3
"""Remote BitMagic supervisor; standard library only, Python 3.9+ and POSIX.

Owns subprocess groups and one checkout lock. Never discards repository changes.
Build/test outputs live in persistent per-cycle directories for investigation.
"""
import argparse
import array
import concurrent.futures
import datetime
import depot
import fcntl
import fireman
import hashlib
import inv_cache
import json
import csv
import inv_reference
import joker
import math
import memory_monitor
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
import urllib.request
import wave


# Initial estimate until a successful run teaches the target for a variant.
PERF_ESTIMATED_TESTS = 299
FOREMAN_VOICE_FILES = {
    'en_AU-librivox-medium.onnx': 'f34849bce64d803d0fe90087766d5d6f2a8215a12edad43e486d77072df7dd41',
    'en_AU-librivox-medium.onnx.json': 'e0d8d381bf3eedeb7483cfa4ab5333d45284d30822bbc93f39fd10dc2935c1c1',
}


def perf_estimated_percent(count, target=PERF_ESTIMATED_TESTS):
    """Estimate completed test elements, keeping active runs below 100 percent."""
    return min(99, max(0, round(count * 100 / target)))


def perf_log_progress(path):
    """Count complete timed records, ignoring headers, blank and partial lines."""
    pattern = re.compile(r';\s*[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?\s+(?:ms|sec|min)\s*$')
    count, last = 0, None
    with path.open('r', errors='replace') as stream:
        for line in stream:
            if line.endswith('\n') and line.split(';', 1)[0].strip().upper() != 'TOTAL' and pattern.search(line):
                count += 1
                last = line.strip()
    return count, last


STRESS_TRANSITION_PHRASES = (
    'The crew is changing engines for the next leg.',
    'We are taking water before the next departure.',
    'The fireman is bunkering coal for the next leg.',
    'The crew is servicing the locomotive before we continue.',
)


def build_status_label(argv, log):
    """Derive display configuration from the actual executable or variant log directory.

    inv_list invokes a shell wrapper with --bin; stress wrappers identify their
    variant through the log's parent directory. Unknown suffixes remain visible.
    Return an empty label when no executable build configuration is available.
    """
    names = []
    if '--bin' in argv:
        position = argv.index('--bin') + 1
        if position < len(argv):
            names.append(Path(str(argv[position])).name)
    names.extend((Path(str(argv[0])).name, log.parent.name))
    display = {'sse42': 'SSE4.2', 'sse41': 'SSE4.1', '64': '64-bit'}
    for name in names:
        name = name.removesuffix('.exe')
        family, separator, configuration = name.partition('_')
        if separator and family in ('bminv', 'perf', 'stress', 'stress64'):
            parts = configuration.split('_')
            return ' '.join(part if part in ('debug', 'release') else
                            display.get(part, part.upper()) for part in parts)
    return ''


def spoken_stress_name(name):
    """Convert a stress executable name into a pronounceable variant label."""
    parts = name.split('_')[1:]
    names = {'sse2': 'S S E two', 'sse42': 'S S E four point two',
             'avx2': 'A V X two', 'avx512': 'A V X five twelve',
             'neon': 'neon', '64': 'sixty four'}
    return ' '.join(names.get(part, part) for part in parts) or 'standard'


ENGINEER_BUILD_REPLIES = (
    'Engineer Joe to Dispatch. Copy that. I am getting steamed up waiting to hit the rails. Over.',
    'Engineer Joe to Dispatch. One build complete? Excellent. This locomotive has places to choo and people to see. Over.',
    'Engineer Joe to Dispatch. Roger. I have been sitting on this track so long I am losing my train of thought. Over.',
    'Engineer Joe to Dispatch. Acknowledged. Let us get this freight moving before I grow roots at the depot. Over.',
    'Engineer Joe to Dispatch. Copy. I am ready to make tracks; even the timetable is getting impatient. Over.',
    'Engineer Joe to Dispatch. All aboard? Good. I have run out of excuses to stay on this siding. Over.',
    'Engineer Joe to Dispatch. Assembly is nearly complete. The crew has begun disassembly to find out why. Over.',
    'Engineer Joe to Dispatch. The assembly crew says everything is working. The disassembly crew says that is not what the code does. Over.',
    'Engineer Joe to Dispatch. We assembled the locomotive twice. The compiler insists on doing it one instruction at a time. Over.',
    'Engineer Joe to Dispatch. The crew found some loose instructions. They are not sure whether to assemble them or put them back in the register. Over.',
    'Engineer Joe to Dispatch. The mechanic asked for an assembly drawing. I gave him the disassembly. We are no longer on speaking terms. Over.',
    'Engineer Joe to Dispatch. The engine passed inspection, but the optimizer removed the caboose. Apparently nobody was using it. Over.',
)

FOREMAN_BUILD_REPLIES = (
    'Depot foreman Angus to Dispatch. The build is finished. Where the hell is that engineer who was so keen to get moving? My crew has been working their backsides off, and Joe has gone walkabout. Send him back to his locomotive. Over.',
    'Depot foreman Angus to Dispatch. The engine is assembled, the crews have worked like mad, and our impatient engineer is nowhere to be found. If Joe is at the tea room, tell him his train is ready and so is my temper. Over.',
    'Depot foreman Angus to Dispatch. We have finished the bloody build. Joe could not wait to leave, could he? Well, we have busted our backsides in this shed, and now the man has vanished. Find that engineer and put him on his engine. Over.',
    'Depot foreman Angus to Dispatch. All build crews report finished. We have been sweating over this locomotive while Joe has apparently gone on an extended inspection of the refreshment room. Tell him to get his backside back here. Over.',
)


DISPATCH_JOKER_REPLIES = (
    'Dispatch to Engineer Joe. Thank you. Please keep the radio clear for railway business.',
    'Dispatch to Engineer Joe. That was noted. Radio silence, please.',
    'Dispatch to Engineer Joe. This frequency is for train operations, not stand-up comedy.',
    'Dispatch to Engineer Joe. Maintain radio discipline. The timetable is suffering enough.',
    'Dispatch to Engineer Joe. Enough. One more joke and you can explain it to the yardmaster.',
    'Dispatch to Engineer Joe. Copy that. Please leave the comedy on the platform.',
    'Dispatch to Engineer Joe. Message received. The signal box requests a quieter channel.',
    'Dispatch to Engineer Joe. Roger. Let us keep this frequency clear for the working timetable.',
    'Dispatch to Engineer Joe. Your remark has been logged. Railway traffic takes priority.',
    'Dispatch to Engineer Joe. Acknowledged. Please return the microphone to railway service.',
    'Dispatch to Engineer Joe. The yardmaster heard that. I suggest we both return to work.',
    'Dispatch to Engineer Joe. Copy. The signal is clear; let us keep the frequency the same way.',
    'Dispatch to Engineer Joe. Understood. Further commentary can wait until the next siding.',
    'Dispatch to Engineer Joe. The train has right of way. Your material does not.',
    'Dispatch to Engineer Joe. Please reserve the airwaves for signals, reports, and genuine emergencies.',
    'Dispatch to Engineer Joe. I will add that to the incident book under unnecessary transmissions.',
    'Dispatch to Engineer Joe. Acknowledged. We are operating a railway, not an open microphone night.',
    'Dispatch to Engineer Joe. Very well. Let us try a stretch of uninterrupted silence.',
    'Dispatch to Engineer Joe. Your audience is currently occupied with the timetable.',
    'Dispatch to Engineer Joe. Noted. Please keep the next transmission strictly operational.',
    'Dispatch to Engineer Joe. I have notified the yardmaster that morale is apparently under control.',
    'Dispatch to Engineer Joe. The line is busy. Save the next one for the crew room.',
    'Dispatch to Engineer Joe. That concludes the entertainment portion of this transmission.',
    'Dispatch to Engineer Joe. Negative on further punchlines. Maintain radio discipline.',
    'Dispatch to Engineer Joe. Please give the other trains a chance to use the air.',
)


class Failure(RuntimeError):
    """An operational failure, or a failed build/test with its recorded command."""


class ReceivedSignal(KeyboardInterrupt):
    """Termination signal captured for logging outside the signal handler."""

    def __init__(self, signum, commit, cycle, children, location):
        """Remember signal and active work before stack unwinding changes it."""
        self.signum = signum
        self.commit = commit
        self.cycle = cycle
        self.children = children
        self.location = location


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
        perf_ready, perf_pending, stress_ready, stress_pending = progress[:4]
        def target_status(name, ready, pending):
            if ready and pending == 0:
                return ' {} build is complete.'.format(name)
            return ' {} builds: {} ready, {} pending.'.format(name, ready, pending)
        status = target_status('Performance', perf_ready, perf_pending)
        status += target_status('Stress', stress_ready, stress_pending)
        if len(progress) > 4:
            if len(progress) > 7 and progress[7]:
                status += ' Sample build is complete.'
            else:
                status += ' Sample build is in progress.'
        if len(progress) > 6:
            status += target_status('Inverted list', progress[5], progress[6])
    if count == 0:
        depot_status = random.choice((
            'The workshop crew reports smooth running in the depot.',
            'The fitters report all is well in the engine shed.',
            'Assembly is proceeding smoothly at the depot.',
            'The repair crew reports no trouble at the workbench.',
            'The locomotive is coming together nicely in the shed.',
            'The depot crew reports a clear inspection so far.',
        ))
        return 'Depot foreman to Dispatch. Build in progress. No compiler warnings so far.' + status + ' ' + depot_status
    word = 'warning' if count == 1 else 'warnings'
    return ('Depot foreman to Dispatch. Attention in the depot. {} compiler {} so far.'.format(count, word)
            + status + ' Please check the build log.')


def arguments():
    """Return validated CLI configuration, defaulting to launcher environment."""
    parser = argparse.ArgumentParser(
        prog='run_tagent.sh',
        description='Build and test BitMagic on new commits, with optional railway sounds.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=('Examples:\n'
                '  ./run_tagent.sh --sound-volume 0.20\n'
                '  ./run_tagent.sh --log-dir /mnt/d/TAgent-logs\n'
                '  ./run_tagent.sh --sound-preview joker\n'
                '  ./run_tagent.sh --once --jobs 4\n'
                '  ./run_tagent.sh --fast --suppress perf,stress --log-dir /mnt/d/TAgent-logs\n'
                '  ./run_tagent.sh --log-dir /mnt/d/TAgent-logs --sound-volume 1.0 --force --log-retention-days 2 --poll-minutes 15 --fast\n\n'
                'Preview events: ' + ', '.join(soundtrack.EVENTS + ('joker', 'joker-reply', 'interrupt', 'station-perf', 'station-stress', 'all', 'perf-progress',
                'repo-restart', 'repo-new', 'repo-build', 'repo-unchanged', 'depot',
                'train-ambient', 'train-normal-ambient', 'train-stop', 'station-manchuria', 'station-visit', 'toronto', 'compilation-depot', 'cache-progress', 'foreman-build-reply', 'fireman', 'inv-progress', 'dispatch', 'depart-perf', 'depart-stress', 'build-clear',
                'build-warning', 'build-final-ok', 'build-final-failed')) + '\n\n'
                'Defaults reflect launcher settings or environment variables; command-line options override them.'))
    numeric_options = [
        ('poll-minutes', 5, 'minutes between remote checks'),
        ('retry-minutes', 60, 'recovery window in minutes for network errors or waiting for a fix'),
        ('reserve-cpus', 2, 'CPUs reserved when calculating the automatic worker budget'),
        ('jobs', 0, 'maximum parallel workers; 0 chooses available CPUs minus reserved CPUs'),
        ('min-free-gb', 5, 'minimum free GiB required on both repository and log filesystems'),
        ('log-retention-days', 7, 'remove eligible old cycle log directories after this many days'),
        ('build-report-seconds', 60, 'initial build status interval; each subsequent interval grows by 10 seconds'),
    ]
    for name, default, description in numeric_options:
        value = float(os.environ.get(name.upper().replace('-', '_'), default))
        parser.add_argument('--' + name, type=float, default=value,
                            help=description + ' (default: %(default)s)')
    for name, default in [('remote', 'origin'), ('branch', 'master')]:
        parser.add_argument('--' + name, default=os.environ.get(name.upper(), default),
                            help='Git ' + name + ' to monitor (default: %(default)s)')
    parser.add_argument('--log-dir', '--state-dir', dest='state_dir', type=Path,
                        default=os.environ.get('LOG_DIR') or None,
                        help='directory for logs and persistent state (default: LOG_DIR or repo/TAgent-logs)')
    parser.add_argument('--startup', choices=('build', 'changed'),
                        default=os.environ.get('STARTUP', 'build'),
                        help='build on every manual start, or only when changed (default: %(default)s)')
    parser.add_argument('--path-to-invref', type=Path, help='Gov2 reference directory; default: discover gov2-reference-* under log directory')
    parser.add_argument('--inv-cache-dir', type=Path, default=os.environ.get('INV_CACHE_DIR') or None,
                        help='optional existing reference cache directory (default: /mnt/c/TAgent-cache on Linux if present)')
    parser.add_argument('--no-inv-cache', action='store_true', help='read inv_list reference directly without a local copy')
    parser.add_argument('--read-memory-mib', '--inv-read-memory-mib', dest='inv_read_memory_mib',
                        type=int, default=os.environ.get('INV_READ_MEMORY_MIB') or None,
                        help='inv_list read-batch budget in MiB; default: 70%% of Linux available RAM; macOS: 8192')
    reader = parser.add_mutually_exclusive_group()
    reader.add_argument('--read-mmap', '--inv-read-mmap', dest='inv_read_mmap', action='store_true',
                        help='use Unix read-only BLOB mapping windows for inv_list (default on Unix, including WSL)')
    reader.add_argument('--no-read-mmap', '--no-inv-read-mmap', dest='inv_read_mmap', action='store_false',
                        help='use buffered inv_list reads instead of Unix mapping windows')
    parser.set_defaults(inv_read_mmap=os.name == 'posix' and os.environ.get('INV_READ_MMAP', 'yes') != 'no')
    parser.add_argument('--suppress', action='append', default=[], metavar='STAGES', help='skip test stages: perf, inv_list, stress; comma-separated or repeat option')
    parser.add_argument('--resume-cycle', action='store_true', help=argparse.SUPPRESS)
    flags = {
        'pedantic': 'run expanded perf and release stress coverage; incompatible with --fast',
        'fast': 'run all supported release perf variants and scalar/AVX2 or NEON release stress in fast mode, then normal matrix; incompatible with --pedantic',
        'no-debug': 'skip debug perf and stress execution; debug builds are still compiled (default: run debug tests)',
        'avx512': 'build and test AVX-512 variants when supported',
        'once': 'run the selected startup cycle once and exit; exit immediately on test failure',
        'check': 'validate local prerequisites only; do not fetch, build, or test',
        'force': 'force a build and test cycle even if the commit was already tested',
        'restart-on-update': 'stop active builds/tests and restart launcher when a new commit is detected',
    }
    for name, description in flags.items():
        parser.add_argument('--' + name, action='store_true', help=description)
    parser.add_argument('--no-sound', action='store_true',
                        default=os.environ.get('SOUND', 'yes') == 'no',
                        help='disable all local sound and speech')
    parser.add_argument('--no-depot', action='store_true',
                        default=os.environ.get('DEPOT_SOUND', 'yes') == 'no',
                        help='disable steampunk depot and train ambience')
    parser.add_argument('--sound-theme', choices=soundtrack.THEMES,
                        default=os.environ.get('SOUND_THEME', 'steampunk'),
                        help='sound theme (default: %(default)s)')
    parser.add_argument('--sound-volume', type=float,
                        default=float(os.environ.get('SOUND_VOLUME', '0.35')),
                        help='overall TAgent sound level from 0.0 (silent) to 1.0 (maximum); '
                             'independent of system volume (default: %(default)s)')
    parser.add_argument('--tts', choices=('auto', 'off'),
                        default=os.environ.get('TTS_MODE', 'auto'),
                        help='auto: use Piper then eSpeak NG when available; off: use prerecorded audio (default: %(default)s)')
    parser.add_argument('--piper-python', default=os.environ.get('PIPER_PYTHON', '~/piper-demo/.venv/bin/python'),
                        help='Python executable for Piper (default: %(default)s)')
    parser.add_argument('--piper-model', default=os.environ.get('PIPER_MODEL') or None,
                        help='Piper dispatcher model (default: LOG_DIR/piper-voices/en_US-lessac-medium.onnx)')
    parser.add_argument('--piper-conductor-model', default=os.environ.get('PIPER_CONDUCTOR_MODEL') or None,
                        help='Piper conductor model (default: LOG_DIR/piper-voices/en_GB-southern_english_female-low.onnx)')
    parser.add_argument('--piper-engineer-model', default=os.environ.get('PIPER_ENGINEER_MODEL') or None,
                        help='Piper engineer model (default: LOG_DIR/piper-voices/en_US-danny-low.onnx)')
    parser.add_argument('--piper-fireman-model', default=os.environ.get('PIPER_FIREMAN_MODEL') or None,
                        help='Piper fireman voice (default: LOG_DIR/piper-voices/en_US-ryan-medium.onnx)')
    parser.add_argument('--piper-foreman-model', default=os.environ.get('PIPER_FOREMAN_MODEL') or None,
                        help='Piper foreman Angus model (default: LOG_DIR/piper-voices/en_AU-librivox-medium.onnx, speaker 4)')
    parser.add_argument('--no-joker', action='store_true', help='disable engineer jokes')
    parser.add_argument('--no-piper-install', action='store_true',
                        default=os.environ.get('PIPER_INSTALL', 'yes') == 'no',
                        help='do not create a user virtual environment or install Piper')
    parser.add_argument('--no-piper-download', action='store_true',
                        default=os.environ.get('PIPER_DOWNLOAD', 'yes') == 'no',
                        help='do not download a missing Piper voice')
    parser.add_argument('--tts-voice', default=os.environ.get('TTS_VOICE', 'en-us+m3'),
                        help='eSpeak voice for crew-count reports (default: %(default)s)')
    parser.add_argument('--dispatcher-voice', default=os.environ.get('DISPATCHER_VOICE', 'en-us+f3,en-us+f4,en-us+m3'),
                        help='comma-separated eSpeak voices for build and remote reports (default: %(default)s)')
    parser.add_argument('--sound-demo', action='store_true',
                        help='Play a complete theme demonstration and exit without Git/builds')
    parser.add_argument('--sound-preview', nargs='?', const='all', metavar='EVENT', choices=soundtrack.EVENTS + ('joker', 'joker-reply', 'interrupt', 'station-perf', 'station-stress', 'all', 'perf-progress', 'repo-restart', 'repo-new', 'repo-build', 'repo-unchanged', 'depot', 'train-ambient', 'train-normal-ambient', 'train-stop', 'station-manchuria', 'station-visit', 'toronto', 'compilation-depot', 'cache-progress', 'foreman-build-reply', 'fireman', 'inv-progress', 'dispatch', 'depart-perf', 'depart-stress', 'build-clear', 'build-warning', 'build-final-ok', 'build-final-failed'),
                        help='play EVENT, or all events if omitted, without Git/build validation; examples: joker, perf-progress, repo-restart, depot, train-ambient')
    args = parser.parse_args()
    args.suppress = {part.strip() for group in args.suppress for part in group.split(',') if part.strip()}
    if args.suppress - {'perf', 'inv_list', 'stress'}:
        parser.error('--suppress accepts perf, inv_list, stress')
    if args.inv_read_memory_mib is not None and not 1 <= args.inv_read_memory_mib <= sys.maxsize // (1024 * 1024):
        parser.error('--read-memory-mib must be positive and fit host address limits')
    if args.inv_read_mmap and os.name != 'posix':
        parser.error('--read-mmap requires Unix')
    if args.fast and args.pedantic:
        parser.error('--fast and --pedantic are mutually exclusive')
    args.full_sound_preview = args.sound_preview == 'all'
    if args.full_sound_preview:
        args.sound_demo = True
        args.sound_preview = None
    if args.tts not in ('auto', 'off'):
        parser.error('TTS_MODE must be auto or off')
    if args.startup not in ('build', 'changed'):
        parser.error('STARTUP must be build or changed')
    if args.sound_theme not in soundtrack.THEMES:
        parser.error('Unknown SOUND_THEME')
    if args.sound_preview in ('depot', 'train-ambient', 'train-normal-ambient', 'train-stop', 'station-visit', 'toronto', 'compilation-depot', 'fireman', 'inv-progress', 'dispatch', 'depart-perf', 'depart-stress', 'build-clear', 'build-warning', 'build-final-ok', 'build-final-failed') and args.sound_theme != 'steampunk':
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


def prepare_log_directory(path, min_free_gb):
    """Create and probe log storage; reject unwritable or undersized destinations.

    path is the resolved destination; min_free_gb is the required free-space
    reserve. Probe creation, writing and removal using the current user's rights.
    """
    try:
        path.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(prefix='.tagent-check-', dir=path) as probe:
            test = Path(probe) / 'write-test'
            with test.open('wb') as stream:
                stream.write(b'TAgent storage check\n')
                stream.flush()
                os.fsync(stream.fileno())
            test.unlink()
        free_gb = shutil.disk_usage(path).free / 1024**3
        if free_gb < min_free_gb:
            raise Failure('Log directory {} has {:.2f} GB free; requires {:.2f} GB'.format(
                path, free_gb, min_free_gb))
    except OSError as error:
        raise Failure('Cannot create or write log directory {}: {}'.format(path, error)) from error


class Agent:
    """Owns configuration, process lifetime, logs and update/recovery state."""

    def __init__(self, args):
        """Initialize local paths and tracking; args is validated CLI configuration."""
        self.args = args
        self.repo = Path(__file__).resolve().parents[2]
        self.home = Path(args.state_dir or self.repo / 'TAgent-logs').expanduser().resolve()
        prepare_log_directory(self.home, args.min_free_gb)
        self.audio_dir = self.home / 'generated-audio'
        self.temp_dir = self.home / 'tmp'
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        self.voices_dir = self.home / 'piper-voices'
        self.voices_dir.mkdir(parents=True, exist_ok=True)
        self.inv_reference = None
        self.inv_progress_time = 0
        self.inv_active = None
        self.memory_monitor = None
        self.memory_phase = 'startup'
        self.test_started_at = None
        self.jemalloc = None
        self.children = set()
        self.mutex = threading.RLock()
        self.log_mutex = threading.RLock()
        self.stopping = threading.Event()
        self.jobs = args.jobs or max(1, available_cpus() - args.reserve_cpus)
        self.last_fetch = 'never'
        self.attempts = 0
        self.commit = 'unknown'
        self.last_sound = 0
        self.last_station_visit = -math.inf
        self.next_intermediate_station = math.inf
        self.station_visit_lock = threading.RLock()
        self.next_train_at = math.inf
        self.audio_thread = None
        self.audio_active = None
        self.audio_epoch = 0
        self.audio_phase_generation = 0
        self.audio_build_generation = 0
        self.audio_closing = False
        self.audio_pending = {}
        self.audio_failure_thread = None
        self.joker_reply_thread = None
        self.dispatch_reply_deck = []
        self.audio_mutex = threading.Lock()
        self.depot = None
        self.audio_warned = False
        self.audio_available = True
        self.console_notice = False
        self.tts_warned = False
        self.build_warnings = WarningCounter()
        self.sample_warnings = WarningCounter()
        self.samples_complete = False
        self.engineer_completed_families = set()
        self.build_lock = threading.Lock()
        self.build_logs = []
        self.build_targets = {}
        self.build_report_interval = args.build_report_seconds
        self.next_build_report = math.inf
        self.tts_engine = (os.environ.get('TAGENT_TTS_ENGINE') if
                           os.environ.get('TAGENT_TTS_CHECKED') == '1' else
                           shutil.which('espeak-ng'))
        self.piper_python = Path(args.piper_python).expanduser()
        self.piper_model = Path(args.piper_model or self.voices_dir / 'en_US-lessac-medium.onnx').expanduser()
        self.piper_conductor_model = Path(args.piper_conductor_model or self.voices_dir /
                                           'en_GB-southern_english_female-low.onnx').expanduser()
        self.piper_available = (self.piper_python.is_file() and os.access(self.piper_python, os.X_OK)
                                and self.piper_model.is_file() and self.piper_model.stat().st_size > 0
                                and Path(str(self.piper_model) + '.json').is_file()
                                and Path(str(self.piper_model) + '.json').stat().st_size > 0)
        self.piper_failed = False
        self.piper_conductor_available = (self.piper_python.is_file() and
            os.access(self.piper_python, os.X_OK) and self.piper_conductor_model.is_file() and
            self.piper_conductor_model.stat().st_size > 0 and
            Path(str(self.piper_conductor_model) + '.json').is_file() and
            Path(str(self.piper_conductor_model) + '.json').stat().st_size > 0)
        self.piper_conductor_failed = False
        self.joker_model = Path(args.piper_engineer_model or self.voices_dir / 'en_US-danny-low.onnx').expanduser()
        self.joker_available = (self.piper_python.is_file() and
            os.access(self.piper_python, os.X_OK) and self.joker_model.is_file() and
            self.joker_model.stat().st_size > 0 and
            Path(str(self.joker_model) + '.json').is_file() and
            Path(str(self.joker_model) + '.json').stat().st_size > 0)
        self.joker_failed = False
        self.engineer_engine = None
        self.engineer_build_replies = []
        self.joker_clock = joker.Schedule()
        self.fireman_clock = fireman.Schedule()
        self.fireman_model = Path(args.piper_fireman_model or self.voices_dir / 'en_US-ryan-medium.onnx').expanduser()
        self.fireman_available = self.fireman_model.is_file() and Path(str(self.fireman_model) + '.json').is_file() and self.piper_python.is_file()
        self.fireman_failed = False
        self.foreman_model = Path(args.piper_foreman_model or self.voices_dir / 'en_AU-librivox-medium.onnx').expanduser()
        self.foreman_available = (self.piper_python.is_file() and os.access(self.piper_python, os.X_OK) and
            self.foreman_model.is_file() and self.foreman_model.stat().st_size > 0 and
            Path(str(self.foreman_model) + '.json').is_file() and Path(str(self.foreman_model) + '.json').stat().st_size > 0)
        self.foreman_failed = False
        self.foreman_engine = None
        self.foreman_build_replies = []
        self.foreman_download_thread = None
        self.conductor_engine = None
        self.conductor_download_thread = None
        self.conductor_download_lock = threading.Lock()
        self.piper_download_thread = None
        self.perf_workers = set()
        self.perf_logs = {}
        self.perf_pending = set()
        self.perf_completed = {}
        self.perf_last_disclaimer = False
        self.perf_last_announced = False
        self.perf_report_lock = threading.Lock()
        self.next_perf_report = math.inf
        self.perf_report_index = 0
        self.station_report_lock = threading.Lock()
        self.station_phase = None
        self.next_station_report = math.inf
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

    def say(self, message, compact=False):
        """Print and persist a timestamped lifecycle message."""
        line = time.strftime('%Y-%m-%d %H:%M:%S ') + message
        with self.log_mutex:
            console_line = line[5:] if compact else line
            print('\r' + console_line, flush=True)
            with (self.home / 'agent.log').open('a') as stream:
                stream.write(line + '\n')

    def prepare_generated_audio(self):
        """Move old generated WAVs without overwriting files or blocking tests."""
        try:
            self.audio_dir.mkdir(parents=True, exist_ok=True)
        except OSError as error:
            self.say('Generated audio directory unavailable; using log root: ' + str(error))
            self.audio_dir = self.home
            return
        for pattern in ('sound-steampunk-*.wav', 'demo-steampunk.wav'):
            for source in self.home.glob(pattern):
                destination = self.audio_dir / source.name
                try:
                    # Exclusive creation prevents replacing an existing destination.
                    with source.open('rb') as src:
                        with destination.open('xb') as dst:
                            try:
                                shutil.copyfileobj(src, dst)
                            except OSError:
                                dst.close()
                                destination.unlink(missing_ok=True)
                                raise
                    source.unlink()
                except FileExistsError:
                    self.say('Generated audio migration conflict; leaving ' + str(source))
                except OSError as error:
                    self.say('Generated audio migration failed for {}: {}'.format(source, error))

    def process_context(self):
        """Describe this process and terminal without claiming a signal sender."""
        try:
            tty = os.ttyname(0)
        except OSError:
            tty = 'none'
        return 'pid={} ppid={} pgid={} sid={} tty={}'.format(
            os.getpid(), os.getppid(), os.getpgrp(), os.getsid(0), tty)

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
                    if self.stress_stop_file is not None and not self.args.restart_on_update:
                        self.stress_stop_file.touch()
                    action = 'stopping active work for restart' if self.args.restart_on_update else 'waiting for current tests to finish'
                    self.say('New commit ' + tip[:12] + '; ' + action)
                else:
                    self.pending_tip = None
                    self.say('REMOTE CHECK OK | no new commits | ' + tip[:12])
                    self.sound('repo-unchanged')
                if tip != self.commit:
                    self.sound('repo-new')
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

    def check_update_failure(self, proc, log):
        """Preserve an observed workload failure before cancelling for an update."""
        with self.mutex:
            failure = getattr(self, 'first_workload_failure', None)
            processes = list(self.children)
        if failure:
            raise failure
        for child in processes:
            child_log = getattr(child, 'tagent_log', None)
            if child_log and child_log.name == 'runner.log' and (child_log.parent / 'first-failure.txt').exists():
                raise CommandFailure(self.stress_failure_details(child_log))
            status = child.poll()
            if child_log and status and not (status == 75 and child_log.name == 'runner.log'):
                raise CommandFailure('{} | exit {} | observed before update cancellation | current log {}'.format(child.args, status, child_log))
        status = proc.poll()
        if log and log.name == 'runner.log' and (log.parent / 'first-failure.txt').exists():
            raise CommandFailure(self.stress_failure_details(log))
        if status and not (status == 75 and log and log.name == 'runner.log'):
            reason = signal.Signals(-status).name if status < 0 else 'exit ' + str(status)
            raise CommandFailure('{} | {} | log {}'.format(proc.args, reason, log))

    def command(self, argv, log=None, cwd=None, timeout=None, fail_fast=False, allowed_exit_codes=(), test_allocator=False):
        """Run argv in an owned group; return stdout or raise with log diagnostics.

        log optionally redirects combined output; cwd defaults to the checkout.
        timeout bounds operational commands. Every exit cleans the process group.
        test_allocator enables the detected jemalloc only for test subprocesses.
        """
        output = None
        proc = None
        try:
            output = log.open('w') if log else subprocess.PIPE
            with self.mutex:
                if self.stopping.is_set():
                    raise Failure('Job cancelled')
                if log and self.pending_tip and (log.name == 'perf.log' or
                        (self.args.restart_on_update and log.name in ('runner.log', 'inv-runner.log', 'build-perf.log', 'build-stress.log', 'build-samples.log', 'build-inv_list.log'))):
                    self.perf_skipped = True
                    return 'skipped-update'
                env = dict(os.environ, TMPDIR=str(self.temp_dir), TMP=str(self.temp_dir), TEMP=str(self.temp_dir),
                           GIT_TERMINAL_PROMPT='0',
                           GIT_SSH_COMMAND='ssh -oBatchMode=yes -oConnectTimeout=20')
                if test_allocator:
                    env = self.test_environment(env)
                proc = subprocess.Popen(argv, cwd=cwd or self.repo, stdout=output,
                                        stderr=subprocess.STDOUT, text=True,
                                        start_new_session=True, env=env)
                proc.tagent_log = log
                self.children.add(proc)
                if log and log.name == 'perf.log':
                    self.perf_workers.add(proc)
                    self.perf_logs[proc] = log
            started = time.monotonic()
            next_job_message = started
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
                        if (log.parent / 'first-failure.txt').exists() or any(line.startswith('FAIL   ') for line in tail.splitlines()):
                            raise CommandFailure('Stress group failed\n' + self.stress_failure_details(log))
                    if timeout and time.monotonic() - started >= timeout:
                        raise CommandFailure('Command timed out: ' + shlex.join(map(str, argv)))
                    if self.stopping.is_set():
                        raise Failure('Job cancelled')
                    if log and log.name in ('build-perf.log', 'build-stress.log', 'build-samples.log', 'build-inv_list.log'):
                        self.build_progress(log)
                        if self.args.restart_on_update:
                            self.poll_remote()
                            if self.pending_tip:
                                self.check_update_failure(proc, log)
                                return 'skipped-update'
                    else:
                        now = time.monotonic()
                        is_test = log and log.name in ('perf.log', 'runner.log', 'inv-runner.log')
                        if not is_test or self.test_started_at is None or now >= next_job_message:
                            if log:
                                mode = 'Fast!' if self.args.fast else 'SloMo' if self.args.pedantic else 'Norm'
                                build = build_status_label(argv, log)
                                self.say(mode + (' ' + build if build else '') + ' | ' + str(log), compact=True)
                            else:
                                self.say('RUN ' + Path(str(argv[0])).name)
                            if is_test and self.test_started_at is not None:
                                elapsed = max(0, now - self.test_started_at)
                                interval = 10 + 290 * min(elapsed / 3600, 1)
                                next_job_message = now + interval
                    if log and log.name in ('perf.log', 'runner.log', 'inv-runner.log'):
                        self.report_intermediate_station()
                        self.report_fireman()
                    if log and log.name == 'inv-runner.log':
                        if self.inv_active:
                            self.report_inv_progress(*self.inv_active)
                        self.poll_remote()
                        if self.args.restart_on_update and self.pending_tip:
                            self.check_update_failure(proc, log)
                            return 'skipped-update'
                    if log and log.name in ('perf.log', 'runner.log'):
                        self.report_station_status()
                        self.report_joker()
                        self.poll_remote()
                        if self.args.restart_on_update and self.pending_tip:
                            self.check_update_failure(proc, log)
                            self.perf_skipped = True
                            return 'skipped-update'
                        if log.name == 'perf.log':
                            self.report_perf_progress()
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
                detail = (self.stress_failure_details(log) if log and log.name == 'runner.log'
                          else str(log) if log else (result or '')[-2000:])
                reason = signal.Signals(-proc.returncode).name if proc.returncode < 0 else 'exit ' + str(proc.returncode)
                raise CommandFailure('{} | {} | variant {} | {}'.format(shlex.join(map(str, argv)), reason, log.parent.name if log else 'operational', detail))
            return result or ''
        except (CommandFailure, OSError) as error:
            if log:
                if not isinstance(error, CommandFailure):
                    error = CommandFailure('{} | variant {} | log {} | {}'.format(argv, log.parent.name, log, error))
                with self.mutex:
                    if not hasattr(self, 'first_workload_failure'):
                        self.first_workload_failure = error
                self.stopping.set()
                self.stop_children()
                raise self.first_workload_failure
            raise
        finally:
            if proc:
                self.signal_group(proc, signal.SIGKILL)
                proc.wait()
                if proc.stdout:
                    proc.stdout.close()
                with self.mutex:
                    self.children.discard(proc)
                    self.perf_workers.discard(proc)
                    self.perf_logs.pop(proc, None)
            if log and output is not None:
                output.close()

    def stress_failure_details(self, log):
        """Describe failed parallel groups and provide operator-ready copy commands.

        log is the stress runner log. Missing or incomplete details preserve the
        runner path so diagnostics never hide the original process failure.
        """
        details = ['Stress variant: ' + log.parent.name, 'Runner log: ' + str(log)]
        try:
            record = log.parent / 'first-failure.txt'
            if record.exists():
                details.append('First failure: ' + record.read_text(errors='replace').strip())
            lines = log.read_text(errors='replace').splitlines()
            for index, line in enumerate(lines):
                if not line.startswith('FAIL   '):
                    continue
                details.append('Failed group: ' + line[len('FAIL   '):].strip())
                for following in lines[index + 1:index + 3]:
                    if following.strip().startswith('Details: '):
                        path = Path(following.strip()[len('Details: '):])
                        details.append('Diagnostic log: ' + str(path))
                        destination = 'TAgent-' + log.parent.name + '-' + path.name
                        details.append('COPY: cp -- ' + shlex.quote(str(path)) + ' ' +
                                       shlex.quote(destination))
                        try:
                            with path.open('rb') as stream:
                                stream.seek(max(0, path.stat().st_size - 8192))
                                tail = stream.read().decode(errors='replace').splitlines()[-12:]
                            details.extend('  ' + entry for entry in tail)
                        except OSError as error:
                            details.append('Cannot read diagnostic log: ' + str(error))
                        break
        except OSError as error:
            details.append('Cannot read runner log: ' + str(error))
        details.append('COPY runner: cp -- ' + shlex.quote(str(log)) + ' ' +
                       shlex.quote('TAgent-' + log.parent.name + '-runner.log'))
        return '\n'.join(details)

    def git(self, *args):
        """Run a bounded noninteractive Git command and return stripped stdout."""
        return self.command(['git', *args], timeout=120).strip()

    def test_environment(self, env):
        """Return a copied test environment, preserving existing preload libraries."""
        env = dict(env)
        if self.jemalloc:
            previous = env.get('LD_PRELOAD', '')
            env['LD_PRELOAD'] = self.jemalloc + (':' + previous if previous else '')
        return env

    def detect_jemalloc(self):
        """Probe Linux's loader cache and validate preload; fall back without failing."""
        self.jemalloc = None
        if platform.system() != 'Linux':
            self.say('ALLOCATOR: jemalloc preload is Linux-only; using system allocator')
            return
        ldconfig = shutil.which('ldconfig')
        if not ldconfig:
            ldconfig = next((p for p in ('/sbin/ldconfig', '/usr/sbin/ldconfig')
                             if os.access(p, os.X_OK)), None)
        if not ldconfig:
            self.say('ALLOCATOR: ldconfig unavailable; using system allocator')
            return
        try:
            listing = subprocess.run([ldconfig, '-p'], capture_output=True, text=True,
                                     timeout=10, check=True).stdout
            for line in listing.splitlines():
                fields = line.split()
                if not fields or fields[0] != 'libjemalloc.so.2' or '=>' not in fields:
                    continue
                library = fields[-1]
                if not Path(library).is_file() or not os.access(library, os.R_OK):
                    continue
                # A separate process avoids loading a second allocator into Python.
                self.jemalloc = library
                probe = subprocess.run(['/bin/true'], env=self.test_environment(os.environ),
                                       capture_output=True, text=True, timeout=10)
                if probe.returncode == 0 and not probe.stderr.strip():
                    self.say('ALLOCATOR: jemalloc enabled for inv_list/perf/stress | ' + library)
                    return
                self.jemalloc = None
            self.say('ALLOCATOR: usable libjemalloc.so.2 not found; using system allocator')
        except (OSError, subprocess.SubprocessError) as error:
            self.jemalloc = None
            self.say('ALLOCATOR: jemalloc check failed; using system allocator | ' + str(error))

    def discover_inv_reference(self):
        """Select one complete reference collection or fail with diagnostics."""
        if 'inv_list' in self.args.suppress:
            self.say('INV REFERENCE | validation suppressed')
            return
        candidates = ([self.args.path_to_invref] if self.args.path_to_invref else
                      sorted(self.home.glob('gov2-reference-*')))
        matches = []
        for candidate in candidates:
            required_levels = range(6, 11) if self.args.fast else range(1, 11) if self.args.pedantic else (6,)
            required_datasets = None if self.args.fast or self.args.pedantic else ('sorted',)
            valid, detail = inv_reference.inspect_reference(candidate, full=self.args.pedantic,
                                                             required_levels=required_levels,
                                                             required_datasets=required_datasets)
            self.say('INV REFERENCE | {} | {}'.format(candidate, detail))
            if valid:
                matches.append(candidate.resolve())
        if len(matches) != 1:
            raise Failure('Expected exactly one complete inv_list reference; use --path-to-invref or --suppress inv_list')
        if self.jobs > 256:
            raise Failure('inv_list runner supports at most 256 workers; set --jobs 256 or lower')
        self.inv_reference = matches[0]
        if self.home == self.inv_reference or self.home.is_relative_to(self.inv_reference):
            raise Failure('inv_list reference cannot contain the TAgent log directory')
        if self.args.no_inv_cache:
            self.say('INV CACHE | disabled; reading ' + str(self.inv_reference))
            return
        cache = self.args.inv_cache_dir
        if cache is None and platform.system() == 'Linux':
            candidate = Path('/mnt/c/TAgent-cache')
            if candidate.is_dir():
                cache = candidate
        if cache is None:
            self.say('INV CACHE | no cache directory configured; reading ' + str(self.inv_reference))
            return
        original = self.inv_reference
        next_announcement = 0.0
        def cache_progress(phase, completed, total, filename):
            """Announce occasional cache progress without delaying the copy."""
            nonlocal next_announcement
            if self.args.check:
                return
            now = time.monotonic()
            if phase == 'ready':
                with self.audio_mutex:
                    self.audio_phase_generation += 1
                self.sound('cache-complete', workers=(
                    'Depot foreman to Dispatch. Reference supplies have arrived at the depot. '
                    'The crew has checked the cargo and the train is ready for assembly. Over.'))
                return
            if now < next_announcement:
                return
            next_announcement = now + 300
            if phase == 'verifying':
                message = ('Depot foreman to Dispatch. The crew is checking the reference cargo. '
                           '{} of {} gibibytes have been copied. Over.'.format(
                               round(completed / 1024**3), round(total / 1024**3)))
            else:
                message = ('Depot foreman to Dispatch. Reference supplies are being loaded at the depot. '
                           '{} of {} gibibytes copied. Over.'.format(
                               round(completed / 1024**3), round(total / 1024**3)))
            self.sound('cache-progress', workers=message)
        try:
            self.inv_reference = inv_cache.copy_reference(original, Path(cache).expanduser(),
                int(self.args.min_free_gb * 1024**3), self.say, cache_progress)
            self.say('INV REFERENCE | using cached copy ' + str(self.inv_reference))
        except (OSError, ValueError) as error:
            with self.audio_mutex:
                self.audio_phase_generation += 1
            if not self.args.check:
                self.sound('cache-complete', workers=(
                    'Depot foreman to Dispatch. Reference loading was interrupted. '
                    'The crew will use the original supplies. Over.'))
            self.inv_reference = original
            self.say('INV CACHE | unavailable: {} | falling back to {}'.format(error, original))

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
        self.detect_jemalloc()
        if (not self.args.check and self.args.sound_theme == 'steampunk' and
                not self.args.no_sound and not self.args.no_depot and
                self.args.sound_volume > 0 and self.audio_available):
            self.depot = depot.Depot(self)
            self.depot.start()
        self.discover_inv_reference()
        self.memory_monitor = memory_monitor.MemoryMonitor(self)
        self.memory_monitor.start()
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
                self.say('No audio player found and console unavailable; playback will be retried')
        if self.args.sound_theme == 'steampunk' and self.args.tts == 'auto':
            if self.piper_available:
                self.say('TTS: Piper available (' + str(self.piper_model) + ')')
            elif self.tts_engine:
                self.say('TTS: eSpeak NG available (' + self.tts_engine + ')')
            else:
                self.tts_warned = True
                self.say('TTS: eSpeak NG unavailable; using prerecorded voice')
            if self.piper_conductor_available:
                self.say('TTS: Conductor Piper voice files present (' + str(self.piper_conductor_model) + ')')

    def start_piper_download(self):
        """Prepare Piper and its voice asynchronously without delaying test work."""
        if (self.args.tts != 'auto' or self.args.no_sound or self.args.sound_volume == 0
                or self.args.check):
            return
        python_exists = self.piper_python.is_file() and os.access(self.piper_python, os.X_OK)
        if not python_exists:
            root = self.piper_python.parent.parent
            safe_new_venv = (self.piper_python.name in ('python', 'python3') and
                             self.piper_python.parent.name == 'bin' and root.name == '.venv' and
                             not root.exists())
            if self.args.no_piper_install or not safe_new_venv:
                self.say('Piper environment unavailable; continuing with TTS fallback')
                return
        if not self.piper_available and self.args.no_piper_download:
            return
        voice = self.piper_model.stem
        if self.piper_model.suffix != '.onnx' or not re.fullmatch(r'[A-Za-z]{2,3}_[A-Za-z]{2}-[A-Za-z0-9_]+-(?:x_low|low|medium|high)', voice):
            self.say('Piper voice path is not a downloadable voice name; using fallback')
            return
        self.say('Checking Piper environment and voice {} in background'.format(voice))
        self.piper_available = False
        self.piper_download_thread = threading.Thread(target=self.download_piper_voice,
                                                      args=(voice,), daemon=True)
        self.piper_download_thread.start()

    def ensure_piper_engine(self):
        """Check Piper and install it only into an explicitly user-owned venv."""
        root = self.piper_python.parent.parent
        if not self.piper_python.is_file():
            if self.args.no_piper_install or root.exists():
                raise Failure('Piper virtual environment unavailable')
            root.parent.mkdir(parents=True, exist_ok=True)
            self.say('Creating Piper virtual environment: ' + str(root))
            self.command([sys.executable, '-m', 'venv', str(root)], timeout=120)
        probe = [str(self.piper_python), '-c',
                 'import importlib.util,sys; sys.exit(0 if importlib.util.find_spec("piper") else 1)']
        try:
            self.command(probe, timeout=20)
            return
        except CommandFailure:
            if self.args.no_piper_install or not (root / 'pyvenv.cfg').is_file():
                raise Failure('Piper is not installed in a writable virtual environment')
        self.say('Installing Piper in virtual environment: ' + str(root))
        self.command([str(self.piper_python), '-m', 'pip', 'install',
                      '--disable-pip-version-check', '--no-cache-dir', 'piper-tts'], timeout=300)
        self.command(probe, timeout=20)

    def start_conductor_download(self):
        """Prepare the conductor's separate Piper voice without delaying tests."""
        if (self.args.tts != 'auto' or self.args.no_sound or self.args.sound_volume == 0
                or self.args.check or self.piper_conductor_available or self.args.no_piper_download):
            return
        with self.conductor_download_lock:
            if self.conductor_download_thread and self.conductor_download_thread.is_alive():
                return
            voice = self.piper_conductor_model.stem
            if (self.piper_conductor_model.suffix != '.onnx' or not re.fullmatch(
                    r'[A-Za-z]{2,3}_[A-Za-z]{2}-[A-Za-z0-9_]+-(?:x_low|low|medium|high)', voice)):
                self.say('Conductor Piper voice path is not a downloadable voice name; using fallback')
                return
            def prepare():
                """Wait for shared engine setup before provisioning the conductor voice."""
                if self.piper_download_thread:
                    while self.piper_download_thread.is_alive() and not self.stopping.is_set():
                        self.piper_download_thread.join(timeout=.25)
                if not self.stopping.is_set():
                    self.download_piper_voice(voice, self.piper_conductor_model, conductor=True)
            self.say('TTS: Preparing conductor Piper voice (' + str(self.piper_conductor_model) + ')')
            self.conductor_download_thread = threading.Thread(target=prepare, daemon=True)
            self.conductor_download_thread.start()

    def download_piper_voice(self, voice, target_model=None, conductor=False, engineer=False, fireman_voice=False):
        """Prepare Piper and stage a missing voice; failures preserve TTS fallback."""
        target_model = target_model or self.piper_model
        try:
            self.ensure_piper_engine()
            config_path = Path(str(target_model) + '.json')
            if (target_model.is_file() and target_model.stat().st_size > 0 and
                    config_path.is_file() and config_path.stat().st_size > 0):
                if fireman_voice:
                    self.fireman_available = True
                    self.fireman_failed = False
                elif engineer:
                    self.joker_available = True
                    self.joker_failed = False
                elif conductor:
                    self.piper_conductor_available = True
                    self.piper_conductor_failed = False
                else:
                    self.piper_available = True
                    self.piper_failed = False
                self.say('TTS: Piper ready (' + str(target_model) + ')')
                return
            if self.args.no_piper_download:
                return
            self.say('Piper voice missing; downloading {}'.format(voice))
            target_model.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(prefix='.tagent-piper-', dir=target_model.parent) as temp:
                directory = Path(temp)
                model = directory / (voice + '.onnx')
                config = Path(str(model) + '.json')
                self.command([str(self.piper_python), '-m', 'piper.download_voices',
                              '--data-dir', str(directory), voice], timeout=120)
                if not model.is_file() or not model.stat().st_size or not config.is_file():
                    raise ValueError('download did not produce both voice files')
                json.loads(config.read_text())
                probe = directory / 'probe.wav'
                self.command([str(self.piper_python), '-m', 'piper', '-m', str(model),
                              '-f', str(probe), '--', 'Voice ready.'], timeout=60)
                with wave.open(str(probe), 'rb') as stream:
                    if stream.getnchannels() != 1 or stream.getsampwidth() != 2 or not stream.getnframes():
                        raise ValueError('Piper produced invalid test audio')
                for source, target in ((model, target_model),
                                       (config, Path(str(target_model) + '.json'))):
                    if not target.is_file() or not target.stat().st_size:
                        os.replace(source, target)
                ready = target_model.is_file() and Path(str(target_model) + '.json').is_file()
                if fireman_voice:
                    self.fireman_available = ready
                    self.fireman_failed = False
                elif engineer:
                    self.joker_available = ready
                    self.joker_failed = False
                elif conductor:
                    self.piper_conductor_available = ready
                    self.piper_conductor_failed = False
                else:
                    self.piper_available = ready
                    self.piper_failed = False
                if ready:
                    self.say('TTS: Piper voice ready (' + str(target_model) + ')')
        except (OSError, ValueError, EOFError, wave.Error, Failure) as error:
            if not self.stopping.is_set():
                self.say('Piper setup or voice download failed ({}); continuing with TTS fallback'.format(error))

    def prepare_foreman_voice(self):
        """Provision and validate Angus asynchronously; leave other voices available on failure."""
        try:
            self.ensure_piper_engine()
            model = self.foreman_model
            config = Path(str(model) + '.json')
            if not (model.is_file() and model.stat().st_size and config.is_file() and config.stat().st_size):
                if self.args.no_piper_download or self.args.piper_foreman_model:
                    return
                self.say('TTS: Downloading Australian foreman voice Angus')
                model.parent.mkdir(parents=True, exist_ok=True)
                base = 'https://huggingface.co/DataCraftsmanAustralia/piper-en_AU-librivox-medium/resolve/main/'
                with tempfile.TemporaryDirectory(prefix='.tagent-foreman-', dir=model.parent) as temp:
                    directory = Path(temp)
                    for name, expected in FOREMAN_VOICE_FILES.items():
                        digest = hashlib.sha256()
                        with urllib.request.urlopen(base + name, timeout=30) as response, (directory / name).open('wb') as output:
                            while not self.stopping.is_set():
                                chunk = response.read(1024 * 1024)
                                if not chunk:
                                    break
                                digest.update(chunk)
                                output.write(chunk)
                        if self.stopping.is_set():
                            return
                        if digest.hexdigest() != expected:
                            raise ValueError('foreman voice checksum mismatch: ' + name)
                    for name, target in ((model.name, model), (config.name, config)):
                        os.replace(directory / name, target)
            settings = json.loads(config.read_text())
            if settings.get('num_speakers', 0) <= 4:
                raise ValueError('foreman model does not contain speaker 4')
            with tempfile.TemporaryDirectory(prefix='foreman-probe-', dir=self.temp_dir) as temp:
                probe = Path(temp) / 'probe.wav'
                self.command([str(self.piper_python), '-m', 'piper', '-m', str(model), '-s', '4',
                              '-f', str(probe), '--', 'Depot foreman ready.'], timeout=60)
                with wave.open(str(probe), 'rb') as audio:
                    if audio.getnchannels() != 1 or audio.getsampwidth() != 2 or not audio.getnframes():
                        raise ValueError('foreman Piper voice produced invalid audio')
            self.foreman_available = True
            self.foreman_failed = False
            self.say('TTS: Foreman Angus Piper ready (' + str(model) + ')')
        except (OSError, ValueError, EOFError, wave.Error, Failure) as error:
            self.foreman_available = False
            if not self.stopping.is_set():
                self.say('TTS: Foreman Angus unavailable ({}); using existing voice fallback'.format(error))

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

    def finish_audio(self, timeout=300):
        """Reject new chatter, then drain the current serialized audio queue."""
        with self.audio_mutex:
            self.audio_closing = True
        return self.wait_for_audio(timeout)

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
            if self.audio_closing:
                return
            self.audio_pending = {key: thread for key, thread in self.audio_pending.items() if thread.is_alive()}
            if event == 'fail' and self.audio_failure_thread and self.audio_failure_thread.is_alive():
                return
            coalesced = ('cache-progress', 'station-perf', 'station-stress', 'station-intermediate', 'perf-progress', 'inv-progress', 'joker', 'fireman', 'work',
                         'dispatch', 'repo-unchanged', 'build-clear', 'build-warning')
            phase_events = ('cache-progress', 'station-perf', 'station-stress', 'station-intermediate', 'perf-progress', 'inv-progress', 'joker', 'fireman', 'work', 'dispatch')
            pending_key = (event, self.audio_phase_generation if event in phase_events else None)
            if event in coalesced and pending_key in self.audio_pending:
                return
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
            phase = self.audio_phase_generation if selected in phase_events else None
            build = self.audio_build_generation if selected in ('build-clear', 'build-warning', 'compilation-depot-start') else None
            self.audio_thread = threading.Thread(target=self.play_after,
                                                args=(previous, selected, workers, epoch, phase, build), daemon=True)
            if event in coalesced:
                self.audio_pending[pending_key] = self.audio_thread
            if event == 'fail':
                self.audio_failure_thread = self.audio_thread
            self.audio_thread.start()

    def play_after(self, previous, event, workers, epoch, phase=None, build=None):
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
            if (epoch != self.audio_epoch or (phase is not None and phase != self.audio_phase_generation) or
                    (build is not None and build != self.audio_build_generation)):
                return
            self.audio_active = threading.current_thread()
        try:
            if not self.audio_available and not (self.args.sound_demo or self.args.sound_preview):
                return
            if self.depot:
                self.depot.mute(.8) if event == 'fireman' else self.depot.mute()
            with self.audio_mutex:
                obsolete = ((phase is not None and phase != self.audio_phase_generation) or
                            (build is not None and build != self.audio_build_generation))
            if self.stopping.is_set() or epoch != self.audio_epoch or obsolete:
                return
            if event == 'interrupt':
                self.play_sound('interrupt-engineer', 1)
                if not self.stopping.is_set() and epoch == self.audio_epoch:
                    self.play_sound('interrupt-dispatch', 1)
                return
            played = self.play_sound(event, workers)
            if event == 'fail' and played and self.args.tts == 'auto':
                for distress_event in ('distress-engineer', 'distress-dispatch'):
                    if self.stopping.wait(.7) or epoch != self.audio_epoch:
                        break
                    self.play_sound(distress_event, 1)
            if event == 'joker' and played and epoch == self.audio_epoch and random.random() >= .25:
                self.joker_reply_thread = threading.Thread(
                    target=self.queue_joker_reply, args=(epoch, phase), daemon=True)
                self.joker_reply_thread.start()
        finally:
            if self.depot:
                self.depot.unmute()
            with self.audio_mutex:
                if self.audio_active is threading.current_thread():
                    self.audio_active = None

    def queue_joker_reply(self, epoch, phase=None):
        """Schedule dispatcher chatter after a pause, unless interrupted or superseded."""
        delay = random.uniform(3, 7)
        deadline = time.monotonic() + delay
        while not self.stopping.is_set() and epoch == self.audio_epoch:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            self.stopping.wait(min(.25, remaining))
        if self.stopping.is_set() or epoch != self.audio_epoch:
            return
        # An important announcement already queued takes precedence over chatter.
        with self.audio_mutex:
            if (self.audio_closing or epoch != self.audio_epoch or
                    (phase is not None and phase != self.audio_phase_generation)):
                return
            previous = self.audio_thread if self.audio_thread and self.audio_thread.is_alive() else None
            if not self.dispatch_reply_deck:
                self.dispatch_reply_deck = list(DISPATCH_JOKER_REPLIES)
                random.shuffle(self.dispatch_reply_deck)
            message = self.dispatch_reply_deck.pop()
            self.say('DISPATCH: ' + message)
            self.audio_thread = threading.Thread(target=self.play_after,
                args=(previous, 'joker-reply', message, epoch, phase), daemon=True)
            self.audio_thread.start()

    def start_test_sound(self, phase='perf'):
        """Announce a test phase without resetting the shared five-minute clock."""
        self.joker_clock.start()
        with self.audio_mutex:
            self.audio_phase_generation += 1
        with self.station_report_lock:
            self.station_phase = phase
            self.next_station_report = time.monotonic() + random.uniform(600, 1200)
        if self.depot:
            self.depot.set_mode('train-normal' if not self.args.fast and not self.args.pedantic else 'train')
        if phase == 'perf':
            self.start_conductor_download()
            self.next_perf_report = time.monotonic() + 600
            self.perf_report_index = 0
        if self.next_train_at == math.inf:
            self.next_train_at = time.monotonic() + 300
        if self.args.sound_theme == 'steampunk':
            self.sound('depart-' + phase)
            self.sound('station-' + phase)

    def announce_interruption(self):
        """Stop existing audio, then play a bounded interrupt exchange if available."""
        if self.args.no_sound or self.args.sound_volume == 0 or self.args.tts != 'auto' or not self.audio_available:
            return
        # The signal handler has set stopping; terminate children before any speech.
        self.stop_children()
        if self.depot:
            self.depot.close()
            self.depot = None
        for thread in (self.audio_thread, self.joker_reply_thread):
            if thread and thread.is_alive():
                thread.join(timeout=5)
                if thread.is_alive():
                    self.say('Interrupt announcement skipped: previous audio did not stop')
                    return
        with self.audio_mutex:
            self.audio_epoch += 1
            self.audio_active = None
            self.audio_thread = None
            self.audio_pending.clear()
            self.audio_failure_thread = None
            self.audio_closing = False
        self.stopping.clear()
        try:
            self.sound('interrupt', force=True)
            self.finish_audio(timeout=90)
        finally:
            self.stopping.set()

    def report_station_status(self):
        """Repeat conductor destination reports every 10–20 minutes during tests."""
        if self.args.no_sound or self.args.sound_volume == 0 or self.args.tts != 'auto':
            return
        if not self.audio_available or self.args.sound_theme != 'steampunk':
            return
        with self.station_report_lock:
            now = time.monotonic()
            if self.station_phase is None or now < self.next_station_report:
                return
            if self.audio_thread and self.audio_thread.is_alive():
                return
            self.next_station_report = now + random.uniform(600, 1200)
            self.sound('station-' + self.station_phase)

    def report_fireman(self):
        """Offer a measured, infrequent boiler report during running tests."""
        if self.args.no_sound or self.args.sound_volume <= 0 or self.args.tts != 'auto' or not self.audio_available:
            return
        if self.audio_thread and self.audio_thread.is_alive():
            return
        if time.monotonic() < self.fireman_clock.deadline:
            return
        try:
            sample = self.memory_monitor.sample() if self.memory_monitor else {}
        except (OSError, ValueError, subprocess.SubprocessError):
            sample = {}
        message = self.fireman_clock.due(sample)
        if message:
            self.say('FIREMAN: ' + message)
            self.sound('fireman', workers=message)

    def report_joker(self):
        """Offer occasional engineer chatter only while test jobs are running."""
        if self.args.no_joker or self.args.no_sound or self.args.tts != 'auto':
            return
        if self.audio_thread and self.audio_thread.is_alive():
            return
        message = self.joker_clock.due()
        if message:
            self.say('ENGINEER: ' + message)
            self.sound('joker', workers=message)

    def perf_target_mode(self, name):
        """Distinguish cycle configuration and the actual fast/normal execution."""
        cycle = 'pedantic' if self.args.pedantic else 'fast' if self.args.fast else 'normal'
        execution = 'fast' if self.args.fast or '_debug' in name else 'normal'
        return cycle + '-' + execution

    def perf_target_path(self, name):
        """Return a semantic cache filename for a supported executable name."""
        if not re.fullmatch(r'perf_[A-Za-z0-9_]+', name):
            raise ValueError('Invalid perf variant name: ' + name)
        return self.home / 'perf-test-targets.cache' / (name + '--' + self.perf_target_mode(name) + '.json')

    def perf_target(self, name):
        """Read the last successful count, or use the initial estimate."""
        try:
            data = json.loads(self.perf_target_path(name).read_text())
            count = data['completed_tests']
            if type(count) is int and count > 0:
                return count
        except (OSError, ValueError, TypeError, KeyError):
            pass
        return PERF_ESTIMATED_TESTS

    def save_perf_target(self, name, path):
        """Atomically retain nonzero complete results after successful execution."""
        try:
            count, _ = perf_log_progress(path)
            if not count:
                return
            target = self.perf_target_path(name)
            target.parent.mkdir(parents=True, exist_ok=True)
            data = {'variant': name, 'mode': self.perf_target_mode(name), 'completed_tests': count}
            with tempfile.NamedTemporaryFile(mode='w', dir=target.parent, prefix='.target-', delete=False) as stream:
                temporary = Path(stream.name)
                try:
                    json.dump(data, stream)
                    stream.write('\n')
                except BaseException:
                    temporary.unlink(missing_ok=True)
                    raise
            try:
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
            self.say('PERF TARGET | {} | {} | {} completed tests'.format(name, self.perf_target_mode(name), count))
        except (OSError, ValueError) as error:
            self.say('PERF TARGET unavailable | {} | {}'.format(name, error))

    def report_perf_progress(self):
        """Print every active perf variant and announce one in round-robin order."""
        if time.monotonic() < self.next_perf_report or not self.perf_report_lock.acquire(False):
            return
        try:
            now = time.monotonic()
            if now < self.next_perf_report:
                return
            with self.mutex:
                logs = sorted(self.perf_logs.values(), key=lambda path: path.parent.name)
            with self.mutex:
                pending = set(self.perf_pending)
                completed = dict(self.perf_completed)
            single = len(pending) == 1 and len(logs) == 1 and logs[0].parent.name in pending
            reports = []
            for path in logs:
                try:
                    count, last = perf_log_progress(path)
                except OSError as error:
                    self.say('PERF STATUS unavailable | {} | {}'.format(path, error))
                    continue
                name = path.parent.name
                reference = None
                target = self.perf_target(name)
                if single and '_debug' in name:
                    preferred = name.replace('_debug', '_release', 1)
                    candidates = ([preferred] if preferred in completed else []) + [
                        item for item in sorted(completed) if '_release' in item and item != preferred]
                    for candidate in candidates:
                        try:
                            ref_count, _ = perf_log_progress(completed[candidate])
                            if ref_count:
                                reference = (candidate, ref_count)
                                break
                        except OSError:
                            continue
                if reference:
                    ref_name, total = reference
                    delayed = count > target
                    detail = ('delayed' if delayed else 'approximately {}%'.format(perf_estimated_percent(count, target)))
                    self.say('PERF STATUS | {} | {} completed | {} | target {} | reference {}: {} | debug fast; release {} | last: {}'.format(
                        name, count, detail, target, ref_name, total, 'fast' if self.args.fast else 'normal', last or 'none yet'))
                    disclaimer = not self.args.fast and not self.perf_last_disclaimer
                    first = not self.perf_last_announced
                    if first:
                        self.say('PERF STATUS | all other perf tasks completed; waiting for ' + name)
                        self.perf_last_announced = True
                    reports.append((name, count, ref_name, total, disclaimer, first, target))
                    if disclaimer:
                        self.perf_last_disclaimer = True
                else:
                    detail = 'delayed' if count > target else 'approximately {}% of test elements'.format(perf_estimated_percent(count, target))
                    self.say('PERF STATUS | {} | {} completed | {} | target {} | last: {}'.format(
                        name, count, detail, target, last or 'none yet'))
                    reports.append((name, count, 'target', target))
            self.next_perf_report = now + 600
            if reports:
                report = reports[self.perf_report_index % len(reports)]
                self.perf_report_index += 1
                self.sound('perf-progress', workers=report)
        finally:
            self.perf_report_lock.release()

    def build_progress(self, log):
        """Report aggregate warnings and completed target counts on a growing interval."""
        with self.build_lock:
            for path in self.build_logs or [log]:
                if path.exists():
                    self.build_warnings.observe(path)
                    if path.name == 'build-samples.log':
                        self.sample_warnings.observe(path)
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
                samples_active = any(path.name == 'build-samples.log' for path in self.build_logs)
                self.say(('BUILD STATUS: {} warning(s) | perf: {} ready, {} pending | '
                          'stress: {} ready, {} pending').format(count, *progress) +
                         (' | samples: {} warning(s)'.format(self.sample_warnings.count) if samples_active else ''))
                if samples_active:
                    progress.append(self.sample_warnings.count)
                if 'inv_list' in self.build_targets:
                    targets = self.build_targets['inv_list']
                    ready = sum(path.is_file() and os.access(path, os.X_OK) for path in targets)
                    self.say('BUILD INV_LIST | {} ready, {} pending'.format(ready, len(targets) - ready))
                    if not samples_active:
                        progress.append(0)
                    progress.extend((ready, len(targets) - ready))
                if samples_active:
                    while len(progress) < 7:
                        progress.append(0)
                    progress.append(self.samples_complete)
                if self.args.sound_theme == 'steampunk':
                    completed = {name for name, i in (('perf', 0), ('stress', 2), ('inv_list', 5))
                                 if i + 1 < len(progress) and progress[i] > 0 and progress[i + 1] == 0}
                    if samples_active and self.samples_complete:
                        completed.add('samples')
                    newly_completed = completed - self.engineer_completed_families
                    self.engineer_completed_families.update(newly_completed)
                    event = ('build-first-clear' if count == 0 else 'build-first-warning') if newly_completed else (
                        'build-clear' if count == 0 else 'build-warning')
                    self.sound(event, workers=(count, *progress))
                    for _ in sorted(newly_completed):
                        self.queue_engineer_build_reply()

    def dispatcher_voice(self):
        """Choose a configured eSpeak voice for this announcement."""
        voices = [voice.strip() for voice in self.args.dispatcher_voice.split(',') if voice.strip()]
        return random.choice(voices or ['en-us+f3'])

    def synthesize_announcement(self, path, text, voice, pitch, workers=1, background='build', conductor=False, engineer=False, distress=False, volume_multiplier=1.0, length_scale=1.0, fireman_voice=False, foreman=False):
        """Mix Piper speech, then eSpeak speech on failure; return None if unavailable.

        voice is a callable so eSpeak voice selection happens only when needed.
        Both engines use agent-owned processes and retain the standard sound fallback.
        """
        if self.args.tts != 'auto':
            return None
        if foreman and (not self.foreman_available or self.foreman_failed):
            foreman = False
            engineer = True
        available = (self.piper_conductor_available and not self.piper_conductor_failed if conductor
                     else self.piper_available and not self.piper_failed)
        model = self.piper_conductor_model if conductor else self.piper_model
        if engineer:
            model = self.joker_model
            available = self.joker_available and not self.joker_failed
        if fireman_voice:
            model = self.fireman_model
            available = self.fireman_available and not self.fireman_failed
        if foreman:
            model = self.foreman_model
            available = self.foreman_available and not self.foreman_failed
        if (not self.tts_engine and not available
                and not self.tts_warned):
            self.tts_warned = True
            self.say('TTS unavailable; using prerecorded voice or sound')
        if not self.tts_engine and not available:
            return None
        if callable(text):
            text = text()
        with tempfile.TemporaryDirectory(prefix='tts-', dir=self.temp_dir) as temp:
            speech = Path(temp) / 'speech.wav'
            if available:
                try:
                    command = [str(self.piper_python), '-m', 'piper', '-m', str(model),
                               '-f', str(speech)]
                    if foreman:
                        command += ['-s', '4']
                    if length_scale != 1.0:
                        command += ['--length-scale', str(length_scale)]
                    self.command(command + ['--', text], timeout=60)
                    duration = soundtrack.mix_dynamic_dispatch(speech, path,
                                                               min(1.0, self.args.sound_volume * (1.25 if distress else volume_multiplier)),
                                                               workers, self.jobs, background=background,
                                                               voice_gain=.85, radio=True)
                    if foreman and self.foreman_engine != 'piper':
                        self.foreman_engine = 'piper'
                        self.say('TTS: Foreman Angus using Piper speaker 4 (' + str(model) + ')')
                    if engineer and self.engineer_engine != 'piper':
                        self.engineer_engine = 'piper'
                        self.say('TTS: Engineer using Piper (' + str(model) + ')')
                    if conductor and self.conductor_engine != 'piper':
                        self.conductor_engine = 'piper'
                        self.say('TTS: Conductor using Piper (' + str(model) + ')')
                    return duration
                except (OSError, ValueError, EOFError, wave.Error, Failure) as error:
                    if foreman:
                        self.foreman_failed = True
                        self.say('TTS: Foreman Angus failed ({}); using existing voice fallback'.format(error))
                        return self.synthesize_announcement(path, text, voice, pitch, workers, background,
                            engineer=True, distress=distress, volume_multiplier=volume_multiplier, length_scale=length_scale)
                    if fireman_voice:
                        self.fireman_failed = True
                    elif engineer:
                        self.joker_failed = True
                    elif conductor:
                        self.piper_conductor_failed = True
                    else:
                        self.piper_failed = True
                    fallback = 'eSpeak NG' if self.tts_engine else 'prerecorded audio or sound'
                    self.say('Piper unavailable or failed ({}); falling back to {}'.format(error, fallback))
                    speech.unlink(missing_ok=True)
            if self.tts_engine:
                try:
                    self.command([self.tts_engine, '-v', voice(), '-s', str(round(145 / length_scale)), '-p', str(pitch),
                                  '-w', str(speech), text], timeout=15)
                    duration = soundtrack.mix_dynamic_dispatch(speech, path,
                                                               min(1.0, self.args.sound_volume * (1.25 if distress else volume_multiplier)),
                                                               workers, self.jobs, background=background,
                                                               voice_gain=.85, radio=True)
                    if engineer and self.engineer_engine != 'espeak':
                        self.engineer_engine = 'espeak'
                        self.say('TTS: Engineer using eSpeak NG fallback')
                    if conductor and self.conductor_engine != 'espeak':
                        self.conductor_engine = 'espeak'
                        self.say('TTS: Conductor using eSpeak NG fallback')
                    return duration
                except (OSError, ValueError, EOFError, wave.Error, Failure) as error:
                    if not self.tts_warned:
                        self.tts_warned = True
                        self.say('TTS unavailable or failed ({}); using sound fallback'.format(error))
        return None

    def render_build_report(self, path, report):
        """Synthesize a status snapshot or use the prerecorded warning fallback."""
        if isinstance(report, tuple):
            count, *progress = report
        else:
            count, progress = report, None
        duration = self.synthesize_announcement(path, lambda: build_report_text(count, progress),
                                                lambda: 'en-us+m3', 40, foreman=True)
        if duration is not None:
            return duration
        asset = Path(soundtrack.__file__).resolve().parent / 'assets' / (
            'build-warning.wav' if count else 'build-clear.wav')
        return soundtrack.mix_dynamic_dispatch(asset, path, self.args.sound_volume,
                                               1, self.jobs, background='build', radio=True)

    def render_build_summary(self, path, result):
        """Announce final success and exact diagnostic counts, with offline fallback."""
        errors, warnings, success = result
        text = ('Depot foreman to Dispatch. Build {}. {} compiler errors. '
                '{} compiler warnings.'.format('complete' if success else 'failed', errors, warnings))
        duration = self.synthesize_announcement(path, text, lambda: 'en-us+m3', 40, foreman=True)
        if duration is not None:
            return duration
        return soundtrack.render_build_summary_fallback(path, errors, warnings, success,
                                                        self.args.sound_volume, temp_dir=self.temp_dir)

    def finish_build(self, success):
        """Print and announce final recognized error/warning line counts."""
        with self.audio_mutex:
            self.audio_build_generation += 1
        counts = (self.build_warnings.errors, self.build_warnings.count, success)
        self.say('BUILD {} | diagnostic errors: {} | warnings: {}'.format(
            'COMPLETE' if success else 'FAILED', counts[0], counts[1]))
        if self.args.sound_theme == 'steampunk':
            self.sound('build-final-ok' if success else 'build-final-failed', workers=counts)
            if success:
                if not self.foreman_build_replies:
                    self.foreman_build_replies = list(FOREMAN_BUILD_REPLIES)
                    random.shuffle(self.foreman_build_replies)
                self.sound('foreman-build-reply', workers=self.foreman_build_replies.pop())
                self.queue_engineer_build_reply()

    def queue_engineer_build_reply(self):
        """Queue a rotating Joe reply after a build report."""
        if not self.engineer_build_replies:
            self.engineer_build_replies = list(ENGINEER_BUILD_REPLIES)
            random.shuffle(self.engineer_build_replies)
        self.sound('engineer-build-reply', workers=self.engineer_build_replies.pop())

    def render_dispatch(self, path, workers):
        """Create a live crew-count announcement, falling back to the bundled WAV.

        path is the output file; workers is the last observed active test-process
        count. A missing engine, bad voice, timeout or invalid WAV never stops tests.
        """
        count = max(1, int(workers))
        crew = 'crew' if count == 1 else 'crews'
        text = 'Dispatch to BitMagic. Freight rolling. {} test {} at work.'.format(count, crew)
        duration = self.synthesize_announcement(path, text, lambda: self.args.tts_voice, 40,
                                                workers=count, background='train')
        if duration is not None:
            return duration
        return soundtrack.render_dispatch(path, self.args.sound_volume)

    def play_sound(self, event, workers):
        """Render and play one event through local audio, then console fallback.

        Called by the audio thread; event and workers are captured on scheduling.
        All playback processes belong to the agent's normal cleanup registry.
        """
        rendered = False
        try:
            suffix = '-' + str(workers) if event == 'work' else ''
            path = self.audio_dir / ('sound-' + self.args.sound_theme + '-' + event + suffix + '.wav')
            if event == 'fireman':
                message = workers if isinstance(workers, str) else 'Fireman Sam to Engineer Joe. We are pulling hard up this grade! I am feeding the fire, but this engine wants every pound of steam we can give her!'
                duration = self.synthesize_announcement(path, message, lambda: 'en-us+m2', 40, background='train', fireman_voice=True, length_scale=1.08)
                if duration is None:
                    return False
                if self.depot is None and (self.args.sound_demo or self.args.sound_preview):
                    depot.mix_preview_background(path, self.args.sound_volume)
            elif event == 'joker':
                duration = self.synthesize_announcement(path, workers, lambda: 'en-us+m3', 40,
                                                        background='train', engineer=True)
                if duration is None:
                    return
            elif event == 'interrupt-engineer':
                text = ('Train engineer to Dispatch. Emergency brakes applied. '
                        'BitMagic freight is stopping. Standing by for restart instructions. Over.')
                duration = self.synthesize_announcement(path, text, lambda: 'en-us+m3', 40,
                                                        background='train', engineer=True)
                if duration is None:
                    return False
            elif event == 'interrupt-dispatch':
                text = ('Dispatch to BitMagic engineer. Brakes acknowledged. '
                        'Hold position and stand by for restart instructions. Over.')
                duration = self.synthesize_announcement(path, text, self.dispatcher_voice, 50,
                                                        background='train')
                if duration is None:
                    return False
            elif event == 'distress-engineer':
                text = ('Train engineer to Dispatch. Mayday, mayday, mayday. '
                        'BitMagic freight train derailed. Requesting an emergency response team. '
                        'Please confirm. Mayday, mayday, over.')
                duration = self.synthesize_announcement(path, text, lambda: 'en-us+m3', 40,
                                                        background='train', engineer=True, distress=True)
                if duration is None:
                    return False
            elif event == 'distress-dispatch':
                text = ('Dispatch to BitMagic engineer. Mayday received and confirmed. '
                        'Emergency response team is on the way. Stand by. Over.')
                duration = self.synthesize_announcement(path, text, self.dispatcher_voice, 50,
                                                        background='train', distress=True)
                if duration is None:
                    return False
            elif event == 'engineer-build-reply':
                duration = self.synthesize_announcement(path, workers, lambda: 'en-us+m3', 40,
                                                        background='build', engineer=True)
                if duration is None:
                    return False
            elif event in ('cache-progress', 'cache-complete', 'foreman-build-reply'):
                self.say('DEPOT FOREMAN: ' + workers)
                duration = self.synthesize_announcement(path, workers, lambda: 'en-us+m3', 40,
                                                        background='build', foreman=True)
                if duration is None:
                    return False
            elif event == 'upgrade-dispatch':
                duration = self.synthesize_announcement(path, workers, self.dispatcher_voice, 50, background='train', volume_multiplier=1.2)
                if duration is None:
                    return False
            elif event == 'upgrade-engineer':
                duration = self.synthesize_announcement(path, workers, lambda: 'en-us+m3', 40,
                                                        background='train', engineer=True, volume_multiplier=1.2)
                if duration is None:
                    return False
            elif event == 'upgrade-conductor':
                duration = self.synthesize_announcement(path, workers, lambda: 'en-us+f3', 50,
                                                        background='train', conductor=True, volume_multiplier=1.2)
                if duration is None:
                    return False
            elif event == 'upgrade-depart':
                duration = soundtrack.render(path, self.args.sound_theme, 'depart',
                                             min(1.0, self.args.sound_volume * 1.2), 1, self.jobs)
            elif event == 'joker-reply':
                duration = self.synthesize_announcement(path, workers, self.dispatcher_voice, 50,
                                                        background='train')
                if duration is None:
                    return False
            elif event in ('inv-announcement', 'inv-progress'):
                self.say('CONDUCTOR: ' + workers)
                duration = self.synthesize_announcement(path, workers, lambda: 'en-us+f3', 50, background='train', conductor=True,
                                                        length_scale=1.18 if event == 'inv-progress' else 1.0)
                if duration is None:
                    return False
            elif event == 'stress-transition':
                message = workers
                self.say('CONDUCTOR: ' + message)
                duration = self.synthesize_announcement(path, message, lambda: 'en-us+f3', 50,
                                                        background='train', conductor=True, volume_multiplier=1.2)
                if duration is None:
                    return False
            elif event in ('station-perf', 'station-stress'):
                duration = self.render_station_report(path, event == 'station-stress')
                if duration is None:
                    return False
            elif event == 'perf-progress':
                duration = self.render_perf_progress(path, workers)
            elif event in ('repo-new', 'repo-build', 'repo-unchanged', 'repo-restart'):
                duration = self.render_repo_report(path, event != 'repo-unchanged',
                                                   event == 'repo-build', event == 'repo-restart')
            elif event in ('compilation-depot-start', 'compilation-depot-end', 'compilation-depot'):
                message = (workers if isinstance(workers, str) else
                           'Dispatch to BitMagic. Welcome to Compilation Depot. The build crews are at work.')
                self.say('DISPATCH: ' + message)
                if event in ('compilation-depot-start', 'compilation-depot'):
                    voice = self.temp_dir / ('depot-station-voice-{}.wav'.format(threading.get_ident()))
                    try:
                        spoken = self.synthesize_announcement(voice, message, self.dispatcher_voice, 50, background='build')
                        duration = soundtrack.render_station_visit(path, voice if spoken is not None else None,
                                                                   self.args.sound_volume, max_seconds=20)
                    finally:
                        voice.unlink(missing_ok=True)
                else:
                    duration = self.synthesize_announcement(path, message, self.dispatcher_voice, 50, background='build')
                    if duration is None:
                        return False
            elif event in ('station-visit', 'station-intermediate', 'toronto'):
                self.say('CONDUCTOR: ' + workers)
                voice = self.temp_dir / ('station-voice-{}.wav'.format(threading.get_ident()))
                try:
                    spoken = self.synthesize_announcement(voice, workers, lambda: 'en-us+f3', 50,
                                                          background='train', conductor=True)
                    duration = soundtrack.render_station_visit(path, voice if spoken is not None else None,
                                                               self.args.sound_volume)
                finally:
                    voice.unlink(missing_ok=True)
            elif event == 'train-stop':
                duration = soundtrack.render_train_stop(path, self.args.sound_volume)
            elif event == 'station-manchuria':
                duration = soundtrack.render_manchuria_station(path, self.args.sound_volume)
            elif event == 'train-normal-demo':
                duration = depot.render_preview(path, self.args.sound_volume, 'train-normal', seconds=30)
            elif event in ('depot', 'train-ambient', 'train-normal-ambient'):
                duration = depot.render_preview(path, self.args.sound_volume,
                                                {'depot': 'depot', 'train-ambient': 'train', 'train-normal-ambient': 'train-normal'}[event])
            elif event in ('depart-perf', 'depart-stress'):
                duration = soundtrack.render(path, self.args.sound_theme, 'depart',
                                             self.args.sound_volume, workers, self.jobs,
                                             announcement=event.split('-')[1])
            elif event in ('build-final-ok', 'build-final-failed'):
                duration = self.render_build_summary(path, workers)
            elif event in ('build-clear', 'build-warning', 'build-first-clear', 'build-first-warning'):
                duration = self.render_build_report(path, workers)
            elif event == 'dispatch':
                duration = self.render_dispatch(path, workers)
            else:
                duration = soundtrack.render(path, self.args.sound_theme, event,
                                             min(1.0, self.args.sound_volume * 1.25) if event == 'fail' else self.args.sound_volume,
                                             workers, self.jobs)
            rendered = True
            # Explicit demos/previews still export WAVs on a headless machine.
            if not self.audio_available:
                return
            for player in ('afplay', 'pw-play', 'paplay', 'aplay'):
                if shutil.which(player):
                    try:
                        # Allow device startup, buffering and scheduler delays on WSL.
                        self.command([player, str(path)], timeout=max(30, duration * 2 + 10))
                        return True
                    except (OSError, Failure) as error:
                        self.say('AUDIO player {} failed for {}: {}'.format(player, event, error))
            if not self.console_notice:
                self.console_notice = True
                self.say('No working audio player; trying the local console bell')
            with open('/dev/console', 'w') as console:
                console.write('\a' * (3 if event == 'fail' else 1))
        except (OSError, ValueError, EOFError, wave.Error, subprocess.SubprocessError) as error:
            if not self.audio_warned:
                self.audio_warned = True
                if rendered:
                    self.say('Local audio unavailable; will retry on the next announcement. Generated WAV files: ' + str(self.audio_dir))
                else:
                    self.say('Sound generation failed: {}; continuing tests'.format(error))

    def render_repo_report(self, path, changed, redirect=False, restart=False):
        """Announce remote status with optional TTS and a distinct sound fallback."""
        message = ('Dispatch to BitMagic. Updates received. Restarting TAgent. Freight is redirected to the build depot for a fresh build and tests.'
                   if restart else 'Dispatch to BitMagic. Updates received. Freight train diverted to the build depot.'
                   if redirect else 'Dispatch to BitMagic. Remote updates detected.' if changed else
                   'Dispatch to BitMagic. No new updates. The line is clear.')
        duration = self.synthesize_announcement(path, message, self.dispatcher_voice, 50)
        if duration is not None:
            return duration
        return soundtrack.render(path, self.args.sound_theme, 'build' if changed else 'idle',
                                 self.args.sound_volume, 1, self.jobs)

    def render_station_report(self, path, stress):
        """Announce the next destination and fast stress service when requested."""
        if stress:
            messages = (
                'Conductor to BitMagic. Our next station is Stress Test Station.',
                'Conductor reporting. Our next destination is Stress Test Station.',
                'Conductor to the crew. We are bound for Stress Test Station.',
                'Conductor reporting. Stress Test Station is ahead.',
            )
        else:
            messages = (
                'Conductor to BitMagic. Our next station is Performance Test Station.',
                'Conductor reporting. Our next destination is Performance Test Station.',
                'Conductor to the crew. We are bound for Performance Test Station.',
                'Conductor reporting. Performance Test Station is ahead.',
            )
        message = random.choice(messages)
        if stress:
            if self.args.fast:
                message += ' Stress tests are running as an express service, with reduced checks.'
            elif self.args.pedantic:
                message += (' We are running the extended inspection schedule. '
                            'Please expect significant delays.')
            else:
                message += (' We are running the regular testing schedule. '
                            'Debug stress variants use reduced checks.')
        elif self.args.fast:
            message += ' Performance tests are running as an express service, with reduced repetitions and fixtures.'
        elif self.args.pedantic:
            message += ' We are running the extended inspection schedule. Please expect significant delays.'
        else:
            message += ' We are running the regular testing schedule. Debug variants use reduced repetitions.'
        self.say('CONDUCTOR: ' + message)
        return self.synthesize_announcement(path, message, lambda: 'en-us+f3', 50,
                                            background='train', conductor=True)

    def render_perf_progress(self, path, report):
        """Speak one variant and its observed count, or play a neutral fallback."""
        name, count = report[:2]
        def female_voice():
            """Choose a configured female eSpeak voice for the conductor."""
            voices = [voice.strip() for voice in self.args.dispatcher_voice.split(',')
                      if '+f' in voice.strip().lower()]
            return random.choice(voices or ['en-us+f3'])
        spoken_name = name.replace('perf_release', 'perf', 1).replace('_', ' ')
        if len(report) == 4 and report[2] == 'target':
            total = report[3]
            if count > total:
                message = 'Conductor reporting. The {} train is delayed. {} tests completed.'.format(spoken_name, count)
            else:
                message = 'Conductor reporting. {} completed {} tests. Approximately {} percent of test elements complete.'.format(
                    spoken_name, count, perf_estimated_percent(count, total))
        elif len(report) > 2:
            ref_name, total, disclaimer = report[2:5]
            target = report[6] if len(report) > 6 else total
            if count > target:
                message = ('Conductor reporting. The {} train is delayed. {} tests completed; {} completed {}.'
                           .format(spoken_name, count, ref_name.replace('_', ' '), total))
            else:
                message = ('Conductor reporting. {} completed {} tests. Approximately {} percent, compared with {} tests on {}.'
                           .format(spoken_name, count, perf_estimated_percent(count, target), total, ref_name.replace('_', ' ')))
            if len(report) > 5 and report[5]:
                message = 'Conductor reporting. All other performance trains have arrived. ' + message.replace('Conductor reporting. ', '', 1)
            if disclaimer:
                message += ' Debug is running fast; the release reference ran on the normal schedule.'
        elif count > PERF_ESTIMATED_TESTS:
            message = ('Conductor reporting. The {} train is delayed. {} tests completed.'
                       .format(spoken_name, count))
        else:
            message = ('Conductor reporting. {} completed {} tests. Approximately {} percent of test elements complete.'
                       .format(spoken_name, count, perf_estimated_percent(count)))
        duration = self.synthesize_announcement(path, message, female_voice, 50, conductor=True)
        if duration is not None:
            return duration
        return soundtrack.render(path, self.args.sound_theme, 'start', self.args.sound_volume, 1, self.jobs)

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
                local = self.git('rev-parse', 'HEAD') if self.commit == 'unknown' else self.commit
                changed = tip != local
                self.say('FETCH OK | ' + ('new commit ' if changed else 'no new commits | ') + tip[:12])
                if not changed:
                    self.sound('repo-unchanged')
                if tip != failed_commit:
                    return tip
                self.say('Waiting for a fix to ' + failed_commit[:12])
            except CommandFailure as error:
                self.say('Fetch failed: ' + str(error))
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                context = ('; test failure details: ' + str(self.home / 'last-test-failure.txt')
                           if failed_commit else '')
                raise Failure('Recovery window exhausted; operator intervention required' + context)
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

    def fast_test_variants(self, family, names):
        """Select scalar release and the supported platform SIMD release."""
        scalar = family + '_release'
        simd = scalar + ('_neon' if self.arm else '_avx2')
        selected = [name for name in (scalar, simd) if name in names and self.supported(name)]
        if scalar not in selected:
            raise Failure('Fast test baseline missing: ' + scalar)
        if simd not in selected:
            self.say('SKIP fast SIMD | {} unavailable or unsupported'.format(simd))
        return selected

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

    def inv_builds(self):
        """Return expected debug/release inverted-list executables for this platform."""
        names = ['bminv_debug', 'bminv_release']
        suffixes = ['neon'] if self.arm else ['sse2', 'sse42', 'avx2']
        if not self.arm and self.args.avx512:
            suffixes.append('avx512')
        for suffix in suffixes:
            names.extend('bminv_' + mode + '_' + suffix for mode in ('debug', 'release'))
        return names

    def report_inv_progress(self, output, name, total):
        """Report completed cases and current phase without blocking validation."""
        now = time.monotonic()
        if now < self.inv_progress_time:
            return
        self.inv_progress_time = now + 600
        spoken_progress = None
        try:
            with (output / 'summary.tsv').open(newline='') as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            complete = sum(row.get('status') == 'passed' for row in rows)
            current = next((row for row in rows if row.get('status') == 'running'), None)
            detail = ('{} level {} {}'.format(current['dataset'], current['level'], current['phase'])
                      if current else 'preparing or finalizing')
            if current and current['phase'] == 'roundtrip-levels' and current.get('log'):
                try:
                    root = output.resolve()
                    phase_log = (root / current['log']).resolve()
                    if phase_log.is_relative_to(root):
                        with phase_log.open('rb') as stream:
                            stream.seek(0, os.SEEK_END)
                            stream.seek(max(0, stream.tell() - 16384))
                            recent = stream.read()
                        matches = re.findall(
                            rb'Progress roundtrip completed=(\d+) total=(\d+) elapsed_s=([0-9.]+)'
                            rb'|Progress verify vectors=(\d+) elapsed_s=([0-9.]+)', recent)
                        if matches:
                            done, vectors, seconds, legacy_done, legacy_seconds = matches[-1]
                            if legacy_done:
                                done, vectors, seconds = legacy_done, current.get('vectors', '0'), legacy_seconds
                            done, vectors = int(done), int(vectors)
                            if vectors > 0 and done <= vectors:
                                percentage = done * 100 // vectors
                                detail += ', {} of {} vectors, {} percent, {} seconds'.format(
                                    done, vectors, percentage, seconds.decode())
                                dataset = current['dataset'].split('.')[-1]
                                spoken_progress = '{} round trip. '.format(dataset.capitalize())
                                if percentage:
                                    spoken_progress += '{} percent complete.'.format(percentage)
                                else:
                                    spoken_progress += 'The train is climbing slowly, but we are making progress.'
                except (OSError, ValueError):
                    pass
        except (OSError, csv.Error, KeyError):
            complete, detail = 0, 'checking reference integrity or preparing'
        message = ('INV PROGRESS | {} | {} of {} cases completed | {}'.format(name, complete, total, detail))
        self.say(message)
        if spoken_progress:
            message = 'Conductor to Dispatch. ' + spoken_progress
        else:
            message = ('Conductor to Dispatch. Inverted list validation. '
                       '{} of {} cases completed. {}.'.format(complete, total, detail.replace('-', ' ')))
        self.sound('inv-progress', workers=message)

    def inv_read_budget(self):
        """Resolve an explicit budget or sample Linux available RAM before each variant."""
        if self.args.inv_read_memory_mib is not None:
            return self.args.inv_read_memory_mib
        if platform.system() != 'Linux':
            return 8192
        try:
            sample = memory_monitor.MemoryMonitor(self).sample()
            available = sample.get('available')
            if available is None:
                raise ValueError('MemAvailable missing')
            limit = sample.get('cgroup_limit')
            current = sample.get('cgroup_current')
            if limit is not None and current is not None:
                available = min(available, max(0, limit - current))
            budget = max(1, min(sys.maxsize // (1024 * 1024), int(available * .7 / (1024 * 1024))))
            self.say('INV READ BUDGET | available {} MiB | 70 percent | selected {} MiB'.format(
                available // (1024 * 1024), budget))
            return budget
        except (OSError, ValueError) as error:
            self.say('INV READ BUDGET | available RAM unavailable: {} | using 512 MiB'.format(error))
            return 512

    def station_visit(self, stage, final=False):
        """Announce a completed stage; limit long station visits to once per 30 minutes."""
        messages = {
            'inv_list': ('Arriving at Validation Junction. Reference checks are complete.',
                         'Pulling into Validation Junction. Reference validation is complete.'),
            'perf': ('Arriving at Performance Terminal. Speed trials are complete.',
                     'Pulling into Performance Terminal. The speed tests are complete.'),
            'stress': ('Arriving at Stress Test Station. The stress run is complete.',
                       'Pulling into Stress Test Station. Stress testing is complete.'),
        }
        message = 'Conductor to Dispatch. ' + (
            'Arriving at Completion Yard. All scheduled tests have passed.' if final else random.choice(messages[stage]))
        with self.station_visit_lock:
            now = time.monotonic()
            if (not self.args.no_sound and self.args.sound_volume > 0 and self.audio_available
                    and now - self.last_station_visit >= 1800):
                self.last_station_visit = now
                self.next_intermediate_station = now + random.uniform(3000, 5400)
                self.sound('station-visit', workers=message)
            else:
                self.sound('inv-announcement', workers=message)

    def report_intermediate_station(self):
        """Schedule an occasional station ambience while test processes are active."""
        if (self.test_started_at is None or self.stopping.is_set() or self.args.no_sound
                or self.args.sound_volume <= 0 or not self.audio_available
                or self.args.sound_theme != 'steampunk'):
            return
        with self.station_visit_lock:
            now = time.monotonic()
            if now < self.next_intermediate_station or now - self.last_station_visit < 1800:
                return
            if self.audio_thread and self.audio_thread.is_alive():
                return
            destination = random.choice(('Long Haul Junction', 'Water Stop', 'Highland Siding',
                                         'Freight Junction', 'Coal Creek Station'))
            self.last_station_visit = now
            self.next_intermediate_station = now + random.uniform(3000, 5400)
            self.sound('station-intermediate', workers=(
                'Conductor to Dispatch. Arriving at {}. Testing continues on schedule.'.format(destination)))

    def test_inv_stage(self, run, manifest):
        """Validate reference collections with supported inv_list binaries in sequence."""
        self.memory_phase = 'inv_list'
        self.joker_clock.start()
        if self.depot:
            self.depot.set_mode('train-normal' if not self.args.fast and not self.args.pedantic else 'train')
        self.sound('inv-announcement', workers='Conductor to Dispatch. Our next station is Inverted List Validation Station. Reference tests are starting.')
        self.sound('depart')
        names = [n for n in manifest['builds']['inv_list'] if self.supported(n)]
        releases = [n for n in names if n == 'bminv_release' or n.startswith('bminv_release_')]
        if not releases:
            raise Failure('No supported inv_list release build is available')
        simd = next((n for n in reversed(releases) if n != 'bminv_release'), None)
        if self.args.fast:
            names = [simd or releases[0]]
        else:
            names = [n for n in ('bminv_release', simd) if n and n in releases]
            if not simd:
                self.say('INV SIMD | no supported SIMD release build available')
        self.say('INV VARIANTS | ' + ', '.join(names))
        for index, name in enumerate(names):
            if index and not self.args.fast and not self.args.pedantic:
                self.sound('train-stop')
            base = run / 'inv_list' / name
            base.mkdir(parents=True)
            output = base / 'results'
            read_budget = self.inv_read_budget()
            argv = ['bash', str(self.repo / 'utils/inv_list/run_inv_test_pl.sh'),
                    '--bin', str(self.repo / 'utils/inv_list' / name),
                    '--reference', str(self.inv_reference), '--out', str(output),
                    '--threads', str(self.jobs), '--silent', '--reserve-gib', str(self.args.min_free_gb),
                    '--read-memory-mib', str(read_budget), '--start-level', '7' if self.args.fast else '1',
                    '--golden-level', '6', '--hash', '--hash-cache-dir', str(self.home / 'inv-hash-cache')]
            if self.args.inv_read_mmap:
                argv.append('--read-mmap')
            if self.args.fast:
                argv.append('--fast')
            elif self.args.pedantic:
                argv.append('--pedantic')
            else:
                argv.extend(['--dataset', 'sorted'])
            self.say('INV START | {} | {} workers | read budget {} MiB | {} | logs {}'.format(
                name, self.jobs, read_budget,
                'mmap' if self.args.inv_read_mmap else 'buffered', base))
            (base / 'command.json').write_text(json.dumps(argv, indent=2) + '\n')
            self.inv_progress_time = 0
            with self.audio_mutex:
                self.audio_phase_generation += 1
            self.inv_active = (output, name, 8 if self.args.fast else 20 if self.args.pedantic else 10)
            try:
                status = self.command(argv, base / 'inv-runner.log', base, test_allocator=True)
            except CommandFailure as error:
                try:
                    details = (base / 'inv-runner.log').read_text(errors='replace')[-4000:]
                except OSError:
                    details = 'Runner log unavailable'
                raise CommandFailure('{} | {}\n{}'.format(name, error, details)) from error
            finally:
                self.inv_active = None
            if status == 'skipped-update' or self.update_waiting():
                return 'updated'
            self.say('INV PASS | {} | {}'.format(name, output / 'summary.tsv'))
        return True

    def test_perf_stage(self, run, manifest):
        """Execute performance variants for this cycle."""
        self.memory_phase = 'perf'
        self.start_test_sound('perf')
        names = ([name for name in manifest['builds']['perf'] if '_release' in name]
                 if self.args.fast else manifest['builds']['perf'])
        with concurrent.futures.ThreadPoolExecutor(max_workers=self.jobs) as pool:
            futures = []
            future_names = {}
            with self.mutex:
                self.perf_pending = {n for n in names if not (self.args.no_debug and '_debug' in n) and self.supported(n)}
                self.perf_completed = {}
                self.perf_last_disclaimer = False
                self.perf_last_announced = False
            for name in names:
                if self.args.no_debug and '_debug' in name:
                    self.say('SKIP debug perf ' + name)
                    continue
                if not self.supported(name):
                    self.say('SKIP perf ' + name)
                    continue
                work = run / name
                work.mkdir()
                perf_command = [str(self.repo / 'tests/perf' / name)]
                if self.args.fast or '_debug' in name:
                    perf_command.append('--fast')
                future = pool.submit(self.command, perf_command, work / 'perf.log', work, test_allocator=True)
                futures.append(future)
                future_names[future] = (name, work / 'perf.log')
            for future in concurrent.futures.as_completed(futures):
                try:
                    result = future.result()
                    with self.mutex:
                        name, path = future_names[future]
                        self.perf_pending.discard(name)
                        if result != 'skipped-update':
                            self.perf_completed[name] = path
                            self.save_perf_target(name, path)
                        if len(self.perf_pending) == 1:
                            self.next_perf_report = 0
                except BaseException as error:
                    self.stopping.set()
                    for pending in futures:
                        pending.cancel()
                    self.stop_children()
                    raise getattr(self, 'first_workload_failure', error)
        if self.update_waiting() or self.perf_skipped:
            self.args.force = True
            return 'updated'
        return True

    def test_stress_stage(self, run, manifest):
        """Execute supported stress variants for this cycle."""
        self.memory_phase = 'stress'
        self.start_test_sound('stress')
        stress_names = []
        names = (self.fast_test_variants('stress', manifest['builds']['stress'])
                 if self.args.fast else manifest['builds']['stress'])
        def stress_mode(name):
            """Choose the actual runner mode for this variant and cycle."""
            if self.args.fast or '_debug' in name or (not self.args.pedantic and name == 'stress_release_gcc'):
                return 'fast'
            return 'pedantic' if self.args.pedantic else 'normal'

        def stress_order(name):
            """In normal cycles put GCC last, after all other release and debug runs."""
            if not self.args.fast and not self.args.pedantic and name == 'stress_release_gcc':
                return 2
            return 1 if '_debug' in name else 0

        for name in sorted(names, key=stress_order):
            if self.args.no_debug and '_debug' in name:
                self.say('SKIP debug stress ' + name)
            elif not self.supported(name):
                self.say('SKIP unsupported ' + name)
            else:
                stress_names.append(name)
        for index, name in enumerate(stress_names):
            if index == 0:
                mode = stress_mode(name)
                self.say('STRESS START | {} | {} mode | {} workers'.format(name, mode, self.jobs))
                self.sound('stress-transition', workers='Conductor to Dispatch. {} stress run is starting. {} mode. {} test crews at work.'.format(spoken_stress_name(name), mode, self.jobs))
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
                    '--order-file', str(self.home / 'stress-order.tsv'),
                    '--stop-file', str(self.stress_stop_file)]
            mode = stress_mode(name)
            if mode == 'fast':
                argv.append('--fast')
            elif mode == 'pedantic':
                argv.append('--pedantic')
            self.say('STRESS ' + name + ' | ' + str(work))
            status = self.command(argv, work / 'runner.log', work, fail_fast=True,
                                  allowed_exit_codes=(75,), test_allocator=True)
            self.stress_stop_file = None
            if status == 'skipped-update':
                self.args.force = True
                return 'updated'
            if status == '75':
                # A branch may have moved back while the workers drained.
                # Rebuild even then, since this stress suite is incomplete.
                self.args.force = True
                self.update_waiting()
                return 'updated'
            if self.update_waiting():
                return 'updated'
            next_name = stress_names[index + 1] if index + 1 < len(stress_names) else None
            completed = '{} stress run completed successfully.'.format(spoken_stress_name(name))
            if next_name:
                if not self.args.fast and not self.args.pedantic:
                    self.sound('train-stop')
                self.sound('station-manchuria')
                core = completed + ' The {} stress run is next.'.format(spoken_stress_name(next_name))
                mode = stress_mode(next_name)
                self.say('STRESS TRANSITION | {} completed | starting {} | {} mode | {} workers'.format(name, next_name, mode, self.jobs))
                message = 'Conductor to Dispatch. {0} {0} {1} The next run is {2} mode, with {3} crews.'.format(core, random.choice(STRESS_TRANSITION_PHRASES), mode, self.jobs)
            else:
                self.say('STRESS COMPLETE | {} | all stress variants completed'.format(name))
                message = 'Conductor to Dispatch. {0} {0} All stress runs are complete.'.format(completed)
            self.sound('stress-transition', workers=message)
        return True

    def cycle(self):
        """Build exact current commit; return False on test failure, True on success."""
        if self.depot:
            self.depot.set_mode('depot')
        run = self.home / (datetime.datetime.now().strftime('%Y%m%d-%H%M%S-%f') + '-' + self.commit[:12])
        run.mkdir()
        self.active_run = run
        marker = self.home / 'active-cycle.json'
        marker.write_text(json.dumps({'commit': self.commit, 'run': str(run), 'stage': 'build', 'status': 'active'}) + '\n')
        self.build_warnings = WarningCounter()
        self.sample_warnings = WarningCounter()
        self.samples_complete = False
        self.engineer_completed_families = set()
        if self.args.restart_on_update:
            self.next_remote_poll = time.monotonic() + self.args.poll_minutes * 60
        self.build_report_interval = self.args.build_report_seconds
        self.next_build_report = time.monotonic() + self.args.build_report_seconds
        manifest = {'tagent_run': 1, 'commit': self.commit, 'platform': platform.platform(),
                    'jemalloc': self.jemalloc, 'jobs': self.jobs, 'pedantic': self.args.pedantic, 'fast': self.args.fast,
                    'no_debug': self.args.no_debug, 'builds': {},
                    'inv_read_memory_mib': self.args.inv_read_memory_mib,
                    'inv_read_mmap': self.args.inv_read_mmap}
        manifest['builds'] = {family: self.builds(family) for family in ('perf', 'stress')}
        (run / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        self.build_logs = [run / ('build-' + family + '.log') for family in ('perf', 'stress', 'samples', 'inv_list')]
        build_tasks = []
        for family in ('perf', 'stress'):
            directory = self.repo / 'tests' / family
            names = self.builds(family)
            manifest['builds'][family] = names
            for name in names:
                (directory / name).unlink(missing_ok=True)
            argv = ['bash', str(directory / 'build_all.sh')]
            if family == 'stress':
                argv += ['--jobs', str(self.jobs), '--build-dir', str(run / 'stress-builds')]
            if self.args.avx512:
                argv.append('--avx512')
            if family == 'perf' and self.args.pedantic:
                argv.append('--pedantic')
            build_tasks.append((family, directory, names, argv, run / ('build-' + family + '.log')))
        inv_dir = self.repo / 'utils' / 'inv_list'
        inv_names = self.inv_builds()
        manifest['builds']['inv_list'] = inv_names
        (run / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        for name in inv_names:
            (inv_dir / name).unlink(missing_ok=True)
        inv_command = ['bash', str(inv_dir / 'build_all.sh')]
        if self.args.avx512:
            inv_command.append('--avx512')
        build_tasks.append(('inv_list', inv_dir, inv_names, inv_command, run / 'build-inv_list.log'))
        sample_jobs = max(1, self.jobs - 3)
        sample_dir = self.repo / 'samples'
        build_tasks.append(('samples', sample_dir, [],
                            ['make', '-C', str(sample_dir), '-B', '-j', str(sample_jobs), 'all'],
                            run / 'build-samples.log'))
        self.build_targets = {family: [directory / name for name in names]
                              for family, directory, names, _, _ in build_tasks}
        self.memory_phase = 'build'
        self.test_started_at = None
        self.sound('building')
        self.sound('compilation-depot-start', workers='Dispatch to BitMagic. Welcome to Compilation Depot. The build crews are at work.')

        def run_build(task):
            """Run one build family and retain its separate diagnostic log."""
            family, _, _, argv, build_log = task
            self.say('BUILD ' + family + ' | ' + str(run))
            try:
                result = self.command(argv, build_log)
                if family == 'samples' and result != 'skipped-update':
                    with self.build_lock:
                        self.samples_complete = True
                return result
            finally:
                with self.build_lock:
                    if build_log.exists():
                        self.build_warnings.observe(build_log, final=True)
                        if family == 'samples':
                            self.sample_warnings.observe(build_log, final=True)

        failures = []
        interrupted_build = False
        with concurrent.futures.ThreadPoolExecutor(max_workers=min(4, self.jobs)) as pool:
            futures = [pool.submit(run_build, task) for task in build_tasks]
            for future in concurrent.futures.as_completed(futures):
                try:
                    if future.result() == 'skipped-update':
                        interrupted_build = True
                except Failure as error:
                    failures.append(error)
        for path in self.build_logs:
            if path.exists():
                self.build_warnings.observe(path, final=True)
        if interrupted_build and not failures:
            self.say('Build interrupted for remote update; fetching before restart')
            return 'updated'
        if failures:
            if self.memory_monitor:
                self.memory_monitor.snapshot('build-failure')
            self.finish_build(False)
            raise getattr(self, 'first_workload_failure', failures[0])
        for _, directory, names, _, _ in build_tasks:
            for name in names:
                if not os.access(directory / name, os.X_OK):
                    self.finish_build(False)
                    raise Failure('Expected build missing: ' + str(directory / name))
        (run / 'manifest.json').write_text(json.dumps(manifest, indent=2))
        self.say('SAMPLES BUILD COMPLETE | warnings: {}'.format(self.sample_warnings.count))
        self.finish_build(True)
        self.test_started_at = time.monotonic()
        self.fireman_clock.start()
        with self.station_visit_lock:
            self.next_intermediate_station = max(self.test_started_at, self.last_station_visit) + random.uniform(3000, 5400)
        self.sound('build')
        self.perf_skipped = False
        # Always check at the build/test boundary; do not defer an overdue check
        # or discard an update or retry state detected during compilation.
        if self.poll_failure_since is None and not self.pending_tip:
            self.next_remote_poll = 0
        if self.update_waiting():
            self.args.force = True
            return 'updated'
        self.sound('compilation-depot-end', workers='Dispatch to BitMagic. Compilation Depot reports all builds complete. The freight is cleared to begin testing.')
        # Test execution order is configured here; each stage owns its scheduling.
        stages = (('inv_list', self.test_inv_stage), ('perf', self.test_perf_stage),
                  ('stress', self.test_stress_stage))
        requested_fast = self.args.fast
        try:
            self.say('TEST MODE | ' + ('Fast!' if self.args.fast else 'SloMo' if self.args.pedantic else 'Norm'))
            previous_stage = None
            enabled_stages = [stage for stage, _ in stages if stage not in self.args.suppress]
            for stage, execute in stages:
                if stage in self.args.suppress:
                    self.say('TEST STAGE SUPPRESSED | ' + stage)
                    continue
                if previous_stage and not self.args.fast and not self.args.pedantic:
                    self.sound('train-stop')
                previous_stage = stage
                marker.write_text(json.dumps({'commit': self.commit, 'run': str(run), 'stage': stage, 'status': 'active'}) + '\n')
                if execute(run, manifest) == 'updated':
                    self.args.force = True
                    return 'updated'
                self.station_visit(stage, final=not requested_fast and stage == enabled_stages[-1])
            if requested_fast:
                # Require a fresh remote check before investing in the longer run.
                self.next_remote_poll = 0
                if self.update_waiting():
                    self.args.force = True
                    return 'updated'
                self.args.fast = False
                if self.depot:
                    self.depot.set_mode('train-normal')
                self.say('TEST MODE | Norm | express cycle complete; reusing builds for normal testing')
                self.sound('upgrade-dispatch', workers='Dispatch to BitMagic crew. Express freight delivered safely. The line is clear for the regular service.')
                self.sound('upgrade-conductor', workers='Conductor to Dispatch. Good news, the express delivery is complete. We are ready for the full regular schedule.')
                self.sound('upgrade-engineer', workers='Engineer Joe to Dispatch. The freight is delivered and the engine is eager for another run. Is the crew bonus on the manifest, or is that still in transit?')
                self.sound('upgrade-depart')
                normal_run = run / 'normal'
                normal_run.mkdir()
                previous_stage = None
                for stage, execute in stages:
                    if stage in self.args.suppress:
                        continue
                    if previous_stage:
                        self.sound('train-stop')
                    previous_stage = stage
                    marker.write_text(json.dumps({'commit': self.commit, 'run': str(normal_run), 'stage': stage, 'mode': 'normal', 'status': 'active'}) + '\n')
                    if execute(normal_run, manifest) == 'updated':
                        self.args.force = True
                        return 'updated'
                    self.station_visit(stage, final=stage == enabled_stages[-1])
                self.next_remote_poll = 0
                if self.update_waiting():
                    self.args.force = True
                    return 'updated'
            return True
        except CommandFailure as error:
            if self.memory_monitor:
                self.memory_monitor.snapshot('test-failure')
            print('TEST FAILURE | commit ' + self.commit + ' | ' + str(error), file=sys.stderr, flush=True)
            try:
                (self.home / 'last-test-failure.txt').write_text(
                    'Commit: ' + self.commit + '\n' + str(error) + '\n')
            except OSError as report_error:
                print('Cannot save failure report: ' + str(report_error), file=sys.stderr, flush=True)
            raise
        finally:
            self.args.fast = requested_fast
            if self.depot:
                self.depot.set_mode('depot')
            self.next_perf_report = math.inf
            with self.audio_mutex:
                self.audio_phase_generation += 1
            with self.station_report_lock:
                self.station_phase = None
                self.next_station_report = math.inf
            self.active_run = None
            self.test_started_at = None
            with self.station_visit_lock:
                self.next_intermediate_station = math.inf
            self.stress_stop_file = None
            self.next_remote_poll = math.inf
            self.next_train_at = math.inf
            self.stopping.clear()

    def run(self):
        """Run update/build/test/wait loop; return exit status, or exec on self-update."""
        self.say('Log directory: ' + str(self.home))
        self.prepare_generated_audio()
        self.say('Process context: ' + self.process_context())
        self.report_audio()
        self.start_piper_download()
        self.start_conductor_download()
        if not self.args.no_sound and self.args.tts == 'auto' and not self.args.check and not self.foreman_available:
            def prepare_foreman():
                """Wait for Piper setup before preparing the optional Angus voice."""
                if self.piper_download_thread:
                    while self.piper_download_thread.is_alive() and not self.stopping.is_set():
                        self.piper_download_thread.join(timeout=.25)
                if not self.stopping.is_set():
                    self.prepare_foreman_voice()
            self.foreman_download_thread = threading.Thread(target=prepare_foreman, daemon=True)
            self.foreman_download_thread.start()
        if not self.args.no_sound and self.args.tts == 'auto' and not self.args.check:
            def prepare_engineer():
                """Reuse the shared engine after its startup provisioning completes."""
                if self.piper_download_thread:
                    self.piper_download_thread.join()
                if not self.stopping.is_set():
                    voice = self.joker_model.stem
                    if not self.joker_available and (self.joker_model.suffix != '.onnx' or not re.fullmatch(
                            r'[A-Za-z]{2,3}_[A-Za-z]{2}-[A-Za-z0-9_]+-(?:x_low|low|medium|high)', voice)):
                        self.say('Engineer Piper voice path is not a downloadable voice name; using fallback')
                        return
                    self.download_piper_voice(voice, self.joker_model, engineer=True)
            engineer_thread = threading.Thread(target=prepare_engineer, daemon=True)
            engineer_thread.start()
            if self.args.sound_preview in ('joker', 'fail', 'interrupt'):
                while engineer_thread.is_alive() and not self.stopping.is_set():
                    engineer_thread.join(timeout=.25)
        if not self.args.no_sound and self.args.tts == 'auto' and not self.args.check:
            def prepare_fireman():
                if self.piper_download_thread:
                    self.piper_download_thread.join()
                if not self.stopping.is_set() and not self.fireman_available:
                    if self.fireman_model.suffix != '.onnx' or not re.fullmatch(
                            r'[A-Za-z]{2,3}_[A-Za-z]{2}-[A-Za-z0-9_]+-(?:x_low|low|medium|high)', self.fireman_model.stem):
                        self.say('Fireman Piper voice path is not a downloadable voice name; using fallback')
                        return
                    self.download_piper_voice(self.fireman_model.stem, self.fireman_model, fireman_voice=True)
            fireman_thread = threading.Thread(target=prepare_fireman, daemon=True)
            fireman_thread.start()
            if self.args.sound_preview == 'fireman' or self.args.sound_demo:
                while fireman_thread.is_alive() and not self.stopping.is_set():
                    fireman_thread.join(timeout=.25)
        if (self.args.sound_demo or self.args.sound_preview) and self.piper_download_thread:
            while self.piper_download_thread.is_alive() and not self.stopping.is_set():
                self.piper_download_thread.join(timeout=.25)
        if self.foreman_download_thread and (self.args.sound_demo or self.args.sound_preview in
                ('cache-progress', 'foreman-build-reply', 'build-clear', 'build-warning', 'build-final-ok', 'build-final-failed')):
            deadline = time.monotonic() + 120
            while self.foreman_download_thread.is_alive() and not self.stopping.is_set() and time.monotonic() < deadline:
                self.foreman_download_thread.join(timeout=.25)
        if self.args.sound_demo or self.args.sound_preview:
            if self.conductor_download_thread:
                while self.conductor_download_thread.is_alive() and not self.stopping.is_set():
                    self.conductor_download_thread.join(timeout=.25)
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
                sequence.extend([('repo-new', 1), ('repo-build', 1), ('repo-unchanged', 1),
                                 ('perf-progress', ('perf_release_64', 12)),
                                 ('build-clear', (0, 2, 4, 3, 7)),
                                 ('build-warning', (3, 5, 1, 7, 3)),
                                 ('build-final-ok', (0, 3, True)),
                                 ('build-final-failed', (2, 3, False))])
                if self.args.full_sound_preview:
                    sequence.extend([('depot', 1), ('train-ambient', 1), ('train-normal-ambient', 1)])
                else:
                    sequence.append(('train-normal-demo', 1))
                sequence.append(('compilation-depot', 'Dispatch to BitMagic. Welcome to Compilation Depot. The build crews are at work.'))
                sequence.append(('cache-progress', 'Depot foreman to Dispatch. Reference supplies are being loaded at the depot. 42 of 167 gibibytes copied. Over.'))
                sequence.append(('foreman-build-reply', FOREMAN_BUILD_REPLIES[0]))
                sequence.append(('train-stop', 1))
                sequence.append(('station-manchuria', 1))
                sequence.append(('fireman', 'Fireman Sam to Engineer Joe. We are pulling hard up this grade! The engine needs every pound of steam we can give her!'))
                sequence.append(('toronto', 'Conductor to Dispatch. Arriving at Long Haul Junction. Testing continues on schedule.'))
            clips = []
            for event, workers in sequence:
                self.say('DEMO ' + event + (' | workers ' + str(workers) if event == 'work' else ''))
                self.sound(event, force=True, workers=workers)
                if self.audio_thread:
                    self.wait_for_audio(timeout=360 if self.args.full_sound_preview else 300)
                suffix = '-' + str(workers) if event == 'work' else ''
                clip = self.audio_dir / ('sound-' + self.args.sound_theme + '-' + event + suffix + '.wav')
                if not self.args.no_sound and self.args.sound_volume > 0 and clip.exists():
                    with wave.open(str(clip), 'rb') as stream:
                        rate = stream.getframerate()
                        frames = stream.readframes(stream.getnframes())
                    if rate == soundtrack.RATE // 2:
                        # The background assets are 22.05 kHz mono; duplicate
                        # samples for the 44.1 kHz combined demonstration.
                        source = array.array('h', frames)
                        converted = array.array('h')
                        for value in source:
                            converted.extend((value, value))
                        frames = converted.tobytes()
                    elif rate != soundtrack.RATE:
                        raise Failure('Unsupported demo WAV sample rate: {}'.format(rate))
                    clips.append(frames)
            if clips:
                with wave.open(str(self.audio_dir / ('demo-' + self.args.sound_theme + '.wav')), 'wb') as stream:
                    stream.setparams((1, 2, soundtrack.RATE, 0, 'NONE', 'not compressed'))
                    stream.writeframes((b'\0' * (soundtrack.RATE // 3 * 2)).join(clips))
            self.say('Sound demo complete | WAV files: ' + str(self.audio_dir))
            return 0
        if self.args.sound_preview:
            preview_count = ('Conductor to Dispatch. Arriving at Performance Terminal. Speed trials are complete.' if self.args.sound_preview in ('station-visit', 'toronto') else
                             'Depot foreman to Dispatch. Reference supplies are being loaded at the depot. 42 of 167 gibibytes copied. Over.' if self.args.sound_preview == 'cache-progress' else
                             random.choice(FOREMAN_BUILD_REPLIES) if self.args.sound_preview == 'foreman-build-reply' else
                             'Conductor to Dispatch. Sorted round trip. The train is climbing slowly, but we are making progress.' if self.args.sound_preview == 'inv-progress' else
                             random.choice(DISPATCH_JOKER_REPLIES) if self.args.sound_preview == 'joker-reply' else
                             'Train engineer Joe to BitMagic. ' + random.choice(self.joker_clock.lines)
                             if self.args.sound_preview == 'joker' else ('perf_release_64', 12) if self.args.sound_preview == 'perf-progress' else
                             (0, 2, 4, 3, 7) if self.args.sound_preview == 'build-clear' else
                             (3, 5, 1, 7, 3) if self.args.sound_preview == 'build-warning' else
                             (0, 3, True) if self.args.sound_preview == 'build-final-ok' else
                             (2, 3, False) if self.args.sound_preview == 'build-final-failed' else self.jobs)
            self.sound(self.args.sound_preview, force=True, workers=preview_count)
            if self.args.sound_preview == 'train-stop':
                message = ('Conductor to Dispatch. SSE two stress run completed successfully. '
                           'The AVX two stress run is next. '
                           'SSE two stress run completed successfully. '
                           'The AVX two stress run is next. '
                           'Taking on water before the next leg. Over.')
                self.sound('stress-transition', force=True, workers=message)
            if self.args.sound_preview == 'joker':
                self.wait_for_audio()
                if self.joker_reply_thread:
                    self.joker_reply_thread.join(timeout=8)
                self.wait_for_audio()
            return 0
        marker = self.home / 'active-cycle.json'
        if marker.exists():
            self.say('Previous cycle did not complete successfully | evidence ' + str(marker) + ' | ' + marker.read_text())
        self.validate()
        if self.args.check:
            return 0
        self.cleanup_logs()
        self.sound('start')
        failed = None
        state_path = self.home / 'success.json'
        startup_pending = True
        while True:
            updated = False
            tip = self.fetch(failed)
            self.commit = self.git('rev-parse', 'HEAD')
            if tip != self.commit:
                if self.git('status', '--porcelain', '--untracked-files=no'):
                    raise Failure('Tracked modifications prevent update')
                agent_changed = bool(self.git('diff', '--name-only', self.commit, tip, '--', 'tests/tagent'))
                self.git('merge', '--ff-only', tip)
                self.commit = tip
                updated = True
                if agent_changed or self.args.restart_on_update:
                    self.say('Updates applied; diverting freight to the build depot')
                    self.sound('repo-restart')
                    self.finish_audio()
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
            if self.inv_reference is not None:
                key['inv_reference'] = str(self.inv_reference)
                key['inv_read_memory_mib'] = self.args.inv_read_memory_mib
                key['inv_read_mmap'] = self.args.inv_read_mmap
            if self.args.suppress:
                key['suppress'] = sorted(self.args.suppress)
            if self.args.no_debug:
                key['no_debug'] = True
            if self.args.fast:
                key['fast'] = True
            previous = json.loads(state_path.read_text()) if state_path.exists() else None
            start_cycle = startup_pending and (self.args.startup == 'build' or self.args.resume_cycle)
            startup_pending = False
            if start_cycle or self.args.force or previous != key or failed:
                if updated:
                    self.say('Updates applied; diverting freight to the build depot')
                    self.sound('repo-build')
                result = self.cycle()
                if result == 'updated':
                    self.args.force = True
                    failed = None
                    continue
                if not result:
                    if self.args.once:
                        return 3
                    return 3
                temp = state_path.with_suffix('.tmp')
                temp.write_text(json.dumps(key))
                temp.replace(state_path)
                (self.home / 'active-cycle.json').unlink(missing_ok=True)
                self.sound('success')
                self.say('SUCCESS ' + self.commit)
                self.args.force = False
                failed = None
            if self.args.once:
                return 0
            self.memory_phase = 'waiting'
            self.sound('idle')
            self.pause(self.args.poll_minutes * 60, 'WAIT')


def main():
    """Start agent, handle signals, and guarantee child cleanup and lock release."""
    try:
        agent = Agent(arguments())
    except (Failure, OSError, ValueError) as error:
        print('TAgent startup ERROR: ' + str(error), file=sys.stderr, flush=True)
        return 1
    def interrupted(signum, frame):
        """Capture the signal and let normal exception cleanup stop children."""
        agent.stopping.set()
        raise ReceivedSignal(signum, agent.commit, agent.active_run, len(agent.children),
                             frame.f_code.co_name if frame else 'unknown')
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    signal.signal(signal.SIGHUP, interrupted)
    try:
        return agent.run()
    except ReceivedSignal as error:
        if agent.memory_monitor:
            agent.memory_monitor.snapshot('signal')
        name = signal.Signals(error.signum).name
        agent.say('SIGNAL {} ({}) received; stopping all child jobs | {} | commit {} | cycle {} | children {} | interrupted {} | sender unavailable'.format(
            name, error.signum, agent.process_context(), error.commit, error.cycle or 'none',
            error.children, error.location))
        agent.announce_interruption()
        return 128 + error.signum
    except KeyboardInterrupt:
        agent.say('Interrupted; stopping all child jobs')
        return 130
    except (Failure, OSError, ValueError) as error:
        error = getattr(agent, 'first_workload_failure', error)
        print('TAgent ERROR ' + str(error) + ' | logs ' + str(agent.home), file=sys.stderr, flush=True)
        try:
            report = agent.home / 'last-test-failure.txt.tmp'
            report.write_text('Commit: ' + agent.commit + '\n' + str(error) + '\n')
            report.replace(agent.home / 'last-test-failure.txt')
        except OSError as report_error:
            print('Cannot save failure report: ' + str(report_error), file=sys.stderr, flush=True)
        agent.stop_children()
        agent.stopping.clear()
        try:
            agent.say('ERROR ' + str(error) + ' | logs ' + str(agent.home))
        except OSError:
            pass
        agent.sound('fail')
        return 3 if 'Recovery window exhausted' in str(error) else 1
    finally:
        if agent.memory_monitor:
            agent.memory_monitor.close()
        if not agent.stopping.is_set():
            agent.finish_audio()
        agent.stopping.set()
        if agent.depot:
            agent.depot.close()
        agent.stop_children()
        if agent.lock_file:
            agent.lock_file.close()


if __name__ == '__main__':
    sys.exit(main())
