#!/usr/bin/env python3
"""Generate immutable, fully validated Gov2 references at compression levels 1–10."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time

HERE = Path(__file__).resolve().parent
APP = HERE.parent
REPO = APP.parents[1]
FIELDS = ['dataset', 'compression_level', 'status', 'stream_file', 'vector_count',
          'integer_count', 'input_bytes', 'input_sha256', 'compressed_bytes', 'index_bytes',
          'bits_per_id', 'compressed_to_input_ratio', 'compression_wall_seconds',
          'threaded_verify_wall_seconds', 'serial_verify_wall_seconds', 'decode_wall_seconds',
          'worker_count', 'batch_vectors', 'logical_cpu_count', 'simd_backend', 'platform',
          'architecture', 'compiler', 'build_flags', 'git_revision', 'working_tree_dirty',
          'source_sha256', 'binary_sha256', 'stream_sha256', 'completed_utc']
DIRECTORY_README = '''# Gov2 reference compression — November 2026

| Content | Meaning |
|---|---|
| `gov2.{sorted,unsorted}.level-NN.bv` | Validated full-collection stream, levels 01–10. |
| Matching `.idx` | BLOB boundaries for threaded reads; excluded from stream bits/ID. |
| `references.tsv` | One row per validated entry: sizes, weighted bits/ID, timings and build provenance. |
| `reports/<entry>/` | Compression/opcode log, threaded and serial raw-ID verification, decode log, opcode TSV, manifest and checksums. |
| `.tools/` | Frozen generator executable and build metadata. |
| `.partial/` | Incomplete work; never a validated reference. |
| `matrix.json` | Input hashes, generator identity and run policy. |

Final names are published only after both readers compare every ID against the
raw collection, independent decoding succeeds, and opcode/index accounting
reconciles. SHA-256 covers streams, indexes and retained reports. Reruns verify
completed entries before skipping them. Existing references are never overwritten.
Bits/ID is `8 * stream_bytes / total_IDs`; timings include profiling and I/O,
exclude hashing, and are historical observations rather than isolated codec benchmarks.

Executable sources and usage: `utils/inv_list/reference/` in the BitMagic repo.
Build separately, then run `run_gov2_reference.sh`; `-j` selects worker count.
'''


def digest(path):
    """Hash a file with bounded memory."""
    result = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def write_json(path, value):
    """Atomically replace a metadata file within its existing directory."""
    temporary = path.with_name(path.name + '.tmp')
    with temporary.open('w') as stream:
        json.dump(value, stream, indent=2)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def number(log, label):
    """Read one unambiguous integer summary value from a phase log."""
    values = re.findall(r'^' + re.escape(label) + r'(\d+)\s*$', log, re.M)
    if len(values) != 1:
        raise RuntimeError(f'Missing or duplicate summary: {label}')
    return int(values[0])


def identity(path):
    """Capture a file identity and hash, rejecting changes during hashing."""
    before = path.stat()
    value = dict(path=str(path), bytes=before.st_size, sha256=digest(path))
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise RuntimeError(f'Input changed while hashing: {path}')
    value['mtime_ns'] = after.st_mtime_ns
    return value


def git_info():
    """Capture revision/dirty state and actual C++/header/build contents, including untracked files."""
    def git(*args):
        return subprocess.check_output(['git', *args], cwd=REPO, text=True).strip()
    source = hashlib.sha256()
    paths = set((REPO / 'src').glob('*.h'))
    paths.update(APP.glob('*.h'))
    paths.update(APP.glob('*.cpp'))
    paths.update([APP / 'Makefile', REPO / 'makefile.in'])
    for path in sorted(paths):
        source.update(str(path.relative_to(REPO)).encode() + b'\0')
        source.update(bytes.fromhex(digest(path)))
    return dict(git_revision=git('rev-parse', 'HEAD'),
                working_tree_dirty=bool(git('status', '--porcelain')),
                source_sha256=source.hexdigest())


def check_volume(args):
    """Reject missing/replaced production mounts; local output is explicitly for fixture testing."""
    if args.allow_local_output:
        return
    if not args.volume.is_mount() or args.volume.stat().st_dev != args.volume_device:
        raise RuntimeError(f'Required reference volume unavailable: {args.volume}')
    if not args.out.is_relative_to(args.volume):
        raise RuntimeError('Output directory must be on the specified volume')


def space(args, required=0):
    """Check current free bytes against a safety reserve plus the next-case estimate."""
    check_volume(args)
    free = shutil.disk_usage(args.out).free
    if free < args.reserve_bytes + required:
        raise RuntimeError(f'Insufficient free space: {free:,} bytes; need {required:,} + reserve {args.reserve_bytes:,}')
    return free


def run_phase(args, command, directory, phase, renderer):
    """Run a phase with live progress, retained logs, elapsed time, and a free-space watchdog."""
    print(f'[{datetime.now().astimezone():%H:%M:%S}] {directory.parent.name}.{phase}', flush=True)
    write_json(directory / (phase + '.command.json'), command)
    errors = []
    start = time.monotonic()
    with (directory / (phase + '.log')).open('wb') as log:
        console = subprocess.Popen([sys.executable, str(renderer), '--hide-opcode-profile'], stdin=subprocess.PIPE)
        child = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        def relay():
            try:
                for line in child.stdout:
                    log.write(line)
                    log.flush()
                    console.stdin.write(line)
                    console.stdin.flush()
            except Exception as error:
                errors.append(error)
            finally:
                child.stdout.close()
                console.stdin.close()
        thread = threading.Thread(target=relay, daemon=True)
        thread.start()
        try:
            while child.poll() is None:
                space(args)
                if errors or console.poll() not in (None, 0):
                    raise RuntimeError('Progress/log relay failed')
                try:
                    child.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    pass
            thread.join()
            console.wait()
            if child.returncode or console.returncode or errors:
                raise RuntimeError(f'{phase} failed: child={child.returncode}, console={console.returncode}, relay={errors}')
        except BaseException:
            child.terminate()
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()
            thread.join(timeout=10)
            console.terminate()
            console.wait()
            raise
    elapsed = time.monotonic() - start
    (directory / (phase + '.time')).write_text(f'{elapsed:.6f} 0 0\n')
    return elapsed


def validate_logs(profile, input_info):
    """Require whole-input accounting and equality via both readers, plus decode completion."""
    logs = {phase: (profile / (phase + '.log')).read_text(errors='replace')
            for phase in ('compress', 'verify', 'verify-serial', 'decode')}
    count = number(logs['compress'], 'Total vectors=')
    ids = number(logs['compress'], 'Total ints=')
    if not count or not ids:
        raise RuntimeError('Empty reference collection')
    for phase in ('compress', 'verify', 'verify-serial'):
        text = logs[phase]
        if number(text, 'Total vectors=') != count or number(text, 'Total ints=') != ids:
            raise RuntimeError(f'{phase}: vector/ID count mismatch')
        if number(text, 'Input consumed bytes=') != input_info['bytes']:
            raise RuntimeError(f'{phase}: not the full input')
        if phase != 'compress' and 'Verification complete.' not in text:
            raise RuntimeError(f'{phase}: missing completion')
    if 'Decode complete.' not in logs['decode'] or number(logs['decode'], 'Total vectors=') != count:
        raise RuntimeError('Incomplete independent decode')
    return count, ids


def completed(out, key, input_info, binary_hash):
    """Verify immutable artifacts and recorded validation, or reject an inconsistent entry."""
    path = out / 'reports' / key / 'metadata.json'
    if not path.exists():
        if (out / (key + '.bv')).exists() or (out / (key + '.idx')).exists() or (out / 'reports' / key).exists():
            raise RuntimeError(f'Uncommitted published artifacts for {key}; inspect before recovery')
        return None
    manifest = json.loads(path.read_text())
    if manifest['status'] != 'validated' or manifest['input_sha256'] != input_info['sha256'] or manifest['binary_sha256'] != binary_hash:
        raise RuntimeError(f'Incompatible completed entry: {key}')
    for relative, checksum in manifest['artifacts'].items():
        artifact = out / relative
        if not artifact.resolve().is_relative_to(out) or not artifact.is_file() or digest(artifact) != checksum:
            raise RuntimeError(f'Reference integrity failure: {artifact}')
    count, ids = validate_logs(out / 'reports' / key, input_info)
    if (count, ids) != (manifest['vector_count'], manifest['integer_count']):
        raise RuntimeError(f'Manifest scope mismatch: {key}')
    if (out / manifest['stream_file']).stat().st_size != manifest['compressed_bytes']:
        raise RuntimeError(f'Manifest size mismatch: {key}')
    expected_stream = path.parent.name + '.bv'
    if manifest['stream_file'] != expected_stream or manifest['status'] != 'validated':
        raise RuntimeError(f'Manifest identity mismatch: {key}')
    return manifest


def write_tsv(out, entries):
    """Regenerate the matrix table atomically from independently validated manifests."""
    temporary = out / 'references.tsv.tmp'
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(sorted(entries, key=lambda row: (row['dataset'], row['compression_level'])))
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, out / 'references.tsv')


def parser():
    """Define production defaults and bounded fixture-test controls."""
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument('--out', type=Path, default=Path('/Volumes/USB_1TB/gov2-reference-Nov-2026'))
    result.add_argument('--volume', type=Path, default=Path('/Volumes/USB_1TB'))
    result.add_argument('--input-dir', type=Path, default=Path('/Volumes/USB_1TB/int_set_bench/Gov2Flat'))
    result.add_argument('--bin', type=Path, help='Optimized binary from build_reference.py; existing matrix reuses its frozen binary')
    result.add_argument('-j', type=int, default=max(1, (os.cpu_count() or 1) - 2))
    result.add_argument('--levels', default='1-10', help='1-10, or a comma-separated subset for staged runs')
    result.add_argument('--datasets', choices=['both', 'sorted', 'unsorted'], default='both')
    result.add_argument('--progress', type=int, default=1000)
    result.add_argument('--reserve-gib', type=float, default=8)
    result.add_argument('--allow-local-output', action='store_true', help='Fixture testing only: bypass production mount enforcement')
    return result


def main(argv=None):
    """Lock the matrix, freeze tools/input identities, then generate missing validated entries."""
    args = parser().parse_args(argv)
    args.out = args.out.resolve()
    args.volume = args.volume.resolve()
    args.input_dir = args.input_dir.resolve()
    if not 1 <= args.j <= 256 or args.progress <= 0 or args.reserve_gib < 0:
        raise RuntimeError('Require 1..256 workers, positive progress and a nonnegative reserve')
    levels = list(range(1, 11)) if args.levels == '1-10' else sorted(set(map(int, args.levels.split(','))))
    if not levels or any(level < 1 or level > 10 for level in levels):
        raise RuntimeError('Levels must be within 1..10')
    sets = ['sorted', 'unsorted'] if args.datasets == 'both' else [args.datasets]
    if not args.allow_local_output and (not args.volume.is_mount() or not args.out.is_relative_to(args.volume)):
        raise RuntimeError('Reference USB volume must be mounted and contain --out')
    args.volume_device = args.volume.stat().st_dev if not args.allow_local_output else None
    args.reserve_bytes = int(args.reserve_gib * 1024**3)
    args.out.mkdir(parents=True, exist_ok=True)
    lock = args.out / '.reference-lock'
    try:
        lock.mkdir()  # Atomic lock: never steal an existing/stale lock automatically.
    except FileExistsError as error:
        raise RuntimeError(f'Matrix already locked: {lock}; inspect owner.json before removing a stale lock') from error
    try:
        write_json(lock / 'owner.json', dict(pid=os.getpid(), host=platform.node()))
        execute(args, sets, levels)
    finally:
        (lock / 'owner.json').unlink(missing_ok=True)
        lock.rmdir()


def execute(args, sets, levels):
    """Execute under an exclusive directory lock; publish a manifest only after all checks."""
    space(args)
    matrix_path = args.out / 'matrix.json'
    inputs = {}
    for dataset in sets:
        path = args.input_dir / ('gov2.' + dataset)
        print(f'Hashing input {path}', flush=True)
        inputs[dataset] = identity(path)
    if matrix_path.exists():
        matrix = json.loads(matrix_path.read_text())
        for dataset, info in inputs.items():
            if dataset in matrix['inputs'] and matrix['inputs'][dataset]['sha256'] != info['sha256']:
                raise RuntimeError(f'Input differs from reference matrix: {dataset}')
        binary = args.out / matrix['binary_file']
        if digest(args.out / '.tools/build.log') != matrix['build']['build_log_sha256']:
            raise RuntimeError('Frozen build log checksum mismatch')
        if digest(binary) != matrix['build']['binary_sha256']:
            raise RuntimeError('Frozen generator checksum mismatch')
        if args.bin and digest(args.bin.resolve()) != matrix['build']['binary_sha256']:
            raise RuntimeError('Requested binary differs from frozen matrix generator; use a new output directory')
        if matrix['build']['architecture'] != platform.machine().lower():
            raise RuntimeError('Frozen generator architecture differs from this host')
        matrix['inputs'].update(inputs)
    else:
        if not args.bin:
            backend = 'neon' if platform.machine().lower() in ('arm64', 'aarch64') else 'avx2'
            args.bin = APP / ('bminv_reference_' + backend)
        original = args.bin.resolve()
        build = json.loads(original.with_suffix('.build.json').read_text())
        if build['architecture'] != platform.machine().lower():
            raise RuntimeError('Build architecture differs from this host')
        if build['backend'] not in ('NEON', 'AVX2') or digest(original) != build['binary_sha256']:
            raise RuntimeError('Optimized build metadata does not match binary; run build_reference.py')
        tools = args.out / '.tools'
        tools.mkdir(exist_ok=True)
        binary = tools / ('bminv-' + build['binary_sha256'])
        if binary.exists():
            raise RuntimeError('Uncommitted tool snapshot exists; inspect output directory')
        shutil.copy2(original, binary)
        build_log = original.with_suffix('.build.log')
        if digest(build_log) != build['build_log_sha256']:
            raise RuntimeError('Build log differs from recorded build metadata')
        shutil.copy2(build_log, tools / 'build.log')
        matrix = dict(schema=1, binary_file=str(binary.relative_to(args.out)), build=build,
                      provenance=build['provenance'], inputs=inputs)
    write_json(matrix_path, matrix)
    readme = args.out / 'readme.md'
    if not readme.exists():
        readme.write_text(DIRECTORY_README)
    entries = []
    # Validate all committed entries before doing work, including unrequested levels.
    for path in sorted((args.out / 'reports').glob('*/metadata.json')):
        stored = json.loads(path.read_text())
        ds = stored['dataset'].split('.')[-1]
        if ds not in inputs:
            inputs[ds] = identity(args.input_dir / ('gov2.' + ds))
        entries.append(completed(args.out, path.parent.name, inputs[ds], matrix['build']['binary_sha256']))
    write_tsv(args.out, entries)
    with tempfile.TemporaryDirectory(prefix='bm-reference-tools-') as temporary:
        frozen = Path(temporary)
        # Freeze the script dependencies as well as the executable for this invocation.
        for directory in ('current', 'opcode_profile'):
            (frozen / directory).mkdir()
        for relative in ('current/progress_console.py', 'current/validate_stream.py', 'opcode_profile/report_opcode_profile.py'):
            shutil.copy2(APP / 'experiments' / relative, frozen / relative)
        for dataset in sets:
            for level in levels:
                key = f'gov2.{dataset}.level-{level:02d}'
                saved = completed(args.out, key, inputs[dataset], matrix['build']['binary_sha256'])
                if saved is not None:
                    print(f'Integrity checked; skipping {key}', flush=True)
                    continue
                # Input bytes are a conservative planning estimate, not a codec bound.
                estimated = inputs[dataset]['bytes']
                free = space(args, estimated)
                remaining = sum(not any(row['stream_file'] == f'gov2.{ds}.level-{lv:02d}.bv' for row in entries)
                                for ds in sets for lv in levels)
                planned = sum(inputs[ds]['bytes'] for ds in sets for lv in levels
                              if not any(row['stream_file'] == f'gov2.{ds}.level-{lv:02d}.bv' for row in entries))
                print(f'Remaining-matrix raw-size planning estimate: {planned / 1024**3:.2f} GiB (not a bound)', flush=True)
                print(f'{key}: free={free / 1024**3:.2f} GiB; {remaining} entries remaining; next-case estimate={estimated / 1024**3:.2f} GiB', flush=True)
                generate(args, binary, frozen, matrix, inputs[dataset], dataset, level, key, entries)
    print(f'Reference table: {args.out / "references.tsv"}', flush=True)


def generate(args, binary, frozen, matrix, info, dataset, level, key, entries):
    """Generate one full-collection case in staging and commit validated artifacts."""
    work = args.out / '.partial' / key
    if work.exists():
        raise RuntimeError(f'Incomplete case exists: {work}; inspect/move it before retrying')
    profile = work / 'profile'
    profile.mkdir(parents=True)
    source = Path(info['path'])
    try:
        read_args = ['-t', '-s', '-bvin', str(profile / 'candidate.bv'), '-bv-index-in', str(profile / 'candidate.idx'),
                     '-threads', str(args.j), '-batch-vectors', '512', '-progress-vectors', str(args.progress)]
        commands = {
            'compress': ['-t', '-s', '-level', str(level), '-threads', str(args.j), '-batch-vectors', '512',
                         '-u32in', str(source), '-bvout', str(profile / 'candidate.bv'),
                         '-bv-index-out', str(profile / 'candidate.idx'), '-profile-opcodes', '-progress-vectors', str(args.progress)],
            'verify': ['-verify', '-u32in', str(source), *read_args],
            'verify-serial': ['-t', '-s', '-verify', '-u32in', str(source), '-bvin', str(profile / 'candidate.bv'),
                              '-progress-vectors', str(args.progress)],
            'decode': ['-decode', *read_args],
        }
        timings = {phase: run_phase(args, [str(binary), *command], profile, phase,
                                   frozen / 'current/progress_console.py') for phase, command in commands.items()}
        count, ids = validate_logs(profile, info)
        for row in entries:
            if row['dataset'] == 'Gov2.' + dataset and (row['vector_count'], row['integer_count']) != (count, ids):
                raise RuntimeError('Collection scope differs across levels')
        if identity(source)['sha256'] != info['sha256']:
            raise RuntimeError('Input content changed during validation')
        with (profile / 'report.txt').open('w') as stream:
            subprocess.run([sys.executable, str(frozen / 'opcode_profile/report_opcode_profile.py'), str(work), '0', 'profile'],
                           stdout=stream, stderr=subprocess.STDOUT, check=True)
        measurement = json.loads((profile / 'measurement.json').read_text())
        with (profile / 'opcodes.tsv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, ['code', 'name', 'count', 'opcode_bytes', 'payload_bytes', 'bytes', 'ids',
                                             'bits_per_int', 'payload_bits_per_int'], delimiter='\t')
            writer.writeheader()
            writer.writerows(measurement['profile']['rows'])
        # Record checksums before promotion. Metadata is the final commit marker.
        stream_path = args.out / (key + '.bv')
        index_path = args.out / (key + '.idx')
        report_path = args.out / 'reports' / key
        if stream_path.exists() or index_path.exists() or report_path.exists():
            raise RuntimeError('Refusing to overwrite existing reference artifacts')
        stream_hash = digest(profile / 'candidate.bv')
        size = (profile / 'candidate.bv').stat().st_size
        row = dict(dataset='Gov2.' + dataset, compression_level=level, status='validated', stream_file=stream_path.name,
                   vector_count=count, integer_count=ids, input_bytes=info['bytes'], input_sha256=info['sha256'],
                   compressed_bytes=size, index_bytes=(profile / 'candidate.idx').stat().st_size,
                   bits_per_id=8 * size / ids, compressed_to_input_ratio=size / info['bytes'],
                   compression_wall_seconds=timings['compress'], threaded_verify_wall_seconds=timings['verify'],
                   serial_verify_wall_seconds=timings['verify-serial'], decode_wall_seconds=timings['decode'],
                   worker_count=args.j, batch_vectors=512, logical_cpu_count=os.cpu_count(),
                   simd_backend=matrix['build']['backend'], platform=platform.platform(), architecture=platform.machine(),
                   compiler=matrix['build']['compiler'].splitlines()[0], build_flags=matrix['build']['build_flags'],
                   binary_sha256=matrix['build']['binary_sha256'], stream_sha256=stream_hash,
                   completed_utc=datetime.now(timezone.utc).isoformat(), **matrix['provenance'])
        artifacts = {stream_path.name: stream_hash, index_path.name: digest(profile / 'candidate.idx')}
        artifacts.update({str((report_path / path.name).relative_to(args.out)): digest(path)
                          for path in profile.iterdir() if path.name not in ('candidate.bv', 'candidate.idx')})
        row['artifacts'] = artifacts
        os.rename(profile / 'candidate.bv', stream_path)
        os.rename(profile / 'candidate.idx', index_path)
        report_path.parent.mkdir(exist_ok=True)
        os.rename(profile, report_path)
        write_json(report_path / 'metadata.json', row)
        work.rmdir()
        entries.append(row)
        write_tsv(args.out, entries)
        print(f'VALIDATED {key}: {size:,} bytes; {8 * size / ids:.6f} bits/ID', flush=True)
    except BaseException as error:
        if work.exists():
            (work / 'FAILURE.txt').write_text(str(error) + '\n')
        raise


def interrupted(signum, frame):
    """Convert termination into normal child cleanup and lock release."""
    raise KeyboardInterrupt(f'Received signal {signum}')


if __name__ == '__main__':
    signal.signal(signal.SIGTERM, interrupted)
    try:
        main()
    except (Exception, KeyboardInterrupt) as error:
        print(f'FAILED: {error}', file=sys.stderr)
        sys.exit(1)
