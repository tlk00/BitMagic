#!/usr/bin/env python3
"""Build stress variants concurrently using the existing platform Make flags.

Each variant owns its objects, dependencies, linked binary and compiler log.
Final binaries are published atomically only after successful linking. The
supervisor forwards complete diagnostic lines once and stops owned process
groups on the first failure or a termination signal.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import platform
import signal
import shutil
import tempfile
import subprocess
import sys
import time


def variants(arm, avx512):
    """Return the established output names and unchanged Make variable arguments."""
    result = [('stress_debug', ['DEBUG=YES']), ('stress_release', [])]
    kernels = [('neon', '-DBMNEONOPT')] if arm else [
        ('sse2', '-DBMSSE2OPT'), ('sse42', '-DBMSSE42OPT'), ('avx2', '-DBMAVX2OPT')]
    if not arm:
        result.append(('stress_release_gcc', ['BMOPTFLAGS=-DBM_USE_GCC_BUILD']))
        if avx512:
            kernels.append(('avx512', '-DBMAVX512OPT'))
    for kernel, flags in kernels:
        result.extend([
            ('stress_debug_' + kernel, ['DEBUG=YES', 'BMOPTFLAGS=' + flags]),
            ('stress_release_' + kernel, ['BMOPTFLAGS=' + flags])])
    return result


def relay(job, final=False):
    """Forward newly written complete log lines, retaining partial diagnostics."""
    chunk = job['reader'].read()
    lines = (job['partial'] + chunk).split(b'\n')
    job['partial'] = lines.pop()
    if final and job['partial']:
        lines.append(job['partial'])
        job['partial'] = b''
    for line in lines:
        sys.stdout.buffer.write(line + b'\n')
    sys.stdout.flush()


def stop_jobs(active):
    """Terminate and reap every owned Make process group, including compilers."""
    for sig in (signal.SIGTERM, signal.SIGKILL):
        for job in active:
            try:
                os.killpg(job['process'].pid, sig)
            except ProcessLookupError:
                pass
        if sig == signal.SIGTERM and active:
            time.sleep(.5)
    for job in active:
        job['process'].wait()


def interrupted(signum, frame):
    """Unwind through cleanup and preserve the conventional signal exit status."""
    raise SystemExit(128 + signum)


def publish(binary, destination):
    """Copy across filesystems, then atomically replace the final executable."""
    descriptor, temporary = tempfile.mkstemp(prefix=destination.name + '-', suffix='.tmp',
                                             dir=destination.parent)
    os.close(descriptor)
    temporary = Path(temporary)
    try:
        shutil.copy2(binary, temporary)
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def build(args, source, matrix):
    """Schedule at most jobs variants; return nonzero on the first failed build."""
    root = args.build_dir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    # The final executable names are shared even when build directories differ.
    with (source / 'build-stress.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another stress build is already running')
        (root / 'failure.json').unlink(missing_ok=True)
        for name, _ in matrix:
            (source / name).unlink(missing_ok=True)
        active = []
        next_variant = 0
        try:
            while active or next_variant < len(matrix):
                # Reap failures before starting more work.
                for job in list(active):
                    relay(job)
                    status = job['process'].poll()
                    if status is None:
                        continue
                    relay(job, final=True)
                    name = job['name']
                    if status:
                        reason = signal.Signals(-status).name if status < 0 else 'exit ' + str(status)
                        print('BUILD FAILED | {} | {} | log {}'.format(name, reason, job['log']), flush=True)
                        (root / 'failure.json').write_text(json.dumps({
                            'variant': name, 'exit_code': status, 'log': str(job['log']),
                            'command': job['command']}) + '\n')
                        return 1
                    publish(job['binary'], source / name)
                    print('BUILD READY | {} | log {}'.format(name, job['log']), flush=True)
                    job['reader'].close()
                    active.remove(job)
                while len(active) < args.jobs and next_variant < len(matrix):
                    name, flags = matrix[next_variant]
                    directory = root / name
                    directory.mkdir(parents=True, exist_ok=True)
                    binary = directory / name
                    log = directory / 'build.log'
                    command = ['make', '-j1', 'BUILD_ROOT=' + str(directory),
                               'DESTINATION=' + str(binary), *flags, 'rebuild']
                    print('BUILD START | {} | log {}'.format(name, log), flush=True)
                    print('COMMAND ' + repr(command), flush=True)
                    # Defer termination until the new process is registered for cleanup.
                    old_mask = signal.pthread_sigmask(signal.SIG_BLOCK,
                                                      {signal.SIGINT, signal.SIGTERM, signal.SIGHUP})
                    try:
                        with log.open('wb') as output:
                            reader = log.open('rb')
                            try:
                                process = subprocess.Popen(command, cwd=source, stdout=output,
                                                           stderr=subprocess.STDOUT, start_new_session=True,
                                                           preexec_fn=lambda: signal.pthread_sigmask(signal.SIG_SETMASK, old_mask))
                            except BaseException:
                                reader.close()
                                raise
                        active.append({'name': name, 'process': process, 'binary': binary,
                                       'log': log, 'command': command, 'reader': reader, 'partial': b''})
                    finally:
                        signal.pthread_sigmask(signal.SIG_SETMASK, old_mask)
                    next_variant += 1
                time.sleep(.1)
            return 0
        finally:
            stop_jobs(active)
            for job in active:
                try:
                    relay(job, final=True)
                finally:
                    job['reader'].close()


def main():
    """Parse build concurrency and platform options and supervise all variants."""
    source = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--jobs', type=int, default=1, help='concurrent variant builds (default: 1)')
    parser.add_argument('--avx512', action='store_true', help='also build x86 AVX-512 variants')
    parser.add_argument('--build-dir', type=Path, default=source / 'build',
                        help='isolated build outputs and logs (default: tests/stress/build)')
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    system, machine = platform.system(), platform.machine()
    arm = system == 'Darwin' and machine in ('arm64', 'aarch64')
    if not arm and (system != 'Linux' or machine not in ('x86_64', 'amd64')):
        parser.error('Unsupported platform: {} {}'.format(system, machine))
    if arm and args.avx512:
        parser.error('--avx512 requires x86 Linux')
    for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
        signal.signal(sig, interrupted)
    try:
        return build(args, source, variants(arm, args.avx512))
    except (OSError, RuntimeError) as error:
        print('Stress build ERROR: ' + str(error), file=sys.stderr, flush=True)
        return 1


if __name__ == '__main__':
    sys.exit(main())
