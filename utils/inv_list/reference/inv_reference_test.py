#!/usr/bin/env python3
"""Validate a target bminv against both Gov2 datasets and saved levels 1–10."""
import argparse
import csv
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import signal
import subprocess
import sys
import threading
import time

from gov2_reference import completed, digest, identity, number, write_json

HERE = Path(__file__).resolve().parent
FIELDS = ['dataset', 'level', 'status', 'phase', 'exit_code', 'reference_bytes',
          'candidate_bytes', 'candidate_minus_reference_bytes', 'candidate_bits_per_id',
          'serialized_bytes_identical', 'vectors', 'ids', 'workers', 'wall_seconds', 'log']


class TestFailure(Exception):
    """Failure carrying the shell exit status of the tested operation."""

    def __init__(self, message, code=1):
        """Retain a contextual message and an exit status in the shell range."""
        super().__init__(message)
        self.code = code


def protect_paths(reference, output):
    """Reject overlapping directory trees, resolving symlink aliases first."""
    if reference.is_relative_to(output) or output.is_relative_to(reference):
        raise TestFailure('Reference and output directories must be disjoint', 2)


def check_space(out, reserve, estimate=0):
    """Enforce free-space reserve before rebuilding and during each operation."""
    if shutil.disk_usage(out).free < reserve + estimate:
        raise TestFailure('Insufficient free output space (including safety reserve)', 2)


def child_limits():
    """Prevent a crashing target from creating large core files in its working directory."""
    resource.setrlimit(resource.RLIMIT_CORE, (0, 0))


def run_phase(binary, flags, out, case, phase, reserve, silent=False):
    """Run one target process; preserve logs and propagate positive exits or signals."""
    directory = out / case
    directory.mkdir(exist_ok=True)
    logfile = directory / (phase + '.log')
    command = [str(binary), *map(str, flags)]
    write_json(directory / (phase + '.command.json'), command)
    if not silent:
        print(f'[{datetime.now().astimezone():%H:%M:%S}] {case}.{phase}', flush=True)
    start = time.monotonic()
    errors = []
    with logfile.open('wb') as log:
        renderer = subprocess.Popen([sys.executable, str(out / '.tools/progress_console.py'),
                                     '--hide-opcode-profile', *(['--silent'] if silent else [])], stdin=subprocess.PIPE, cwd=out)
        # preexec_fn is safe here: the relay thread starts after Popen.
        child = subprocess.Popen(command, cwd=out, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                 preexec_fn=child_limits)
        def relay():
            """Retain every target line and forward a compact console view."""
            try:
                for line in child.stdout:
                    log.write(line)
                    log.flush()
                    renderer.stdin.write(line)
                    renderer.stdin.flush()
            except Exception as error:
                errors.append(str(error))
            finally:
                child.stdout.close()
                try:
                    renderer.stdin.close()
                except BrokenPipeError:
                    pass
        thread = threading.Thread(target=relay, daemon=True)
        thread.start()
        try:
            while child.poll() is None:
                check_space(out, reserve)
                if errors or renderer.poll() not in (None, 0):
                    raise TestFailure(f'Console/log relay failed: {errors}')
                try:
                    child.wait(timeout=1)
                except subprocess.TimeoutExpired:
                    pass
            thread.join(timeout=10)
            if thread.is_alive():
                raise TestFailure('Log relay did not finish after target exit')
            renderer.wait(timeout=10)
        except BaseException:
            child.terminate()
            try:
                child.wait(timeout=10)
            except subprocess.TimeoutExpired:
                child.kill()
                child.wait()
            thread.join(timeout=10)
            renderer.terminate()
            try:
                renderer.wait(timeout=10)
            except subprocess.TimeoutExpired:
                renderer.kill()
                renderer.wait()
            raise
    elapsed = time.monotonic() - start
    result = dict(exit_code=child.returncode, wall_seconds=elapsed, log=str(logfile.relative_to(out)))
    write_json(directory / (phase + '.result.json'), result)
    if child.returncode:
        if child.returncode < 0:
            sig = -child.returncode
            reason = f'terminated by {signal.Signals(sig).name}'
            code = 128 + sig
        else:
            reason = f'returned exit {child.returncode}'
            code = child.returncode
        raise TestFailure(f'{case} phase={phase}: bminv {reason}; log={logfile}', code)
    if renderer.returncode or errors:
        raise TestFailure(f'{case} phase={phase}: console/log relay failed; log={logfile}')
    return logfile.read_text(errors='replace'), result


def validate_read(text, row, phase, raw_bytes):
    """Require expected completion, full vector count, and expected-input scope where applicable."""
    marker = 'Decode complete.' if 'decode' in phase else 'Verification complete.'
    if marker not in text or number(text, 'Total vectors=') != row['vector_count']:
        raise TestFailure(f'{phase}: incomplete or incorrect decoded vector count')
    if 'decode' not in phase:
        if number(text, 'Total ints=') != row['integer_count']:
            raise TestFailure(f'{phase}: expected ID count mismatch')
        if number(text, 'Input consumed bytes=') != raw_bytes:
            raise TestFailure(f'{phase}: incomplete expected-input comparison')


def save_summary(out, rows):
    """Replace only this test run's summary TSV, never the reference manifest."""
    temporary = out / 'summary.tsv.tmp'
    with temporary.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, FIELDS, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, out / 'summary.tsv')


def validate_index(path, size, count):
    """Check an entire generated BLOB index without assuming serialized byte identity."""
    with path.open() as stream:
        if stream.readline().split() != ['BM_INV_BLOB_INDEX_V1', str(size), str(count)]:
            raise TestFailure('Rebuilt index header mismatch')
        previous = -1
        entries = 0
        for line in stream:
            offset = int(line)
            if offset <= previous or (entries == 0 and offset != 0):
                raise TestFailure('Rebuilt index offsets are invalid')
            previous = offset
            entries += 1
    if entries != count + 1 or previous != size:
        raise TestFailure('Rebuilt index scope mismatch')


def execute(args, references, binary):
    """Validate each golden against raw IDs before comparing selected reference/candidate streams."""
    rows = []
    reserve = int(args.reserve_gib * 1024**3)
    progress = ['-progress-vectors', str(args.progress)]
    baselines = {row['dataset']: row for row in references
                 if row['compression_level'] == args.golden_level}
    validated = set()
    for row in references:
        if row['compression_level'] < args.start_level:
            continue
        dataset = row['dataset'].split('.')[-1]
        case = Path(row['stream_file']).stem
        source = args.reference / 'Gov2Flat' / ('gov2.' + dataset)
        original = args.reference / row['stream_file']
        original_index = original.with_suffix('.idx')
        baseline = baselines[row['dataset']]
        golden = args.reference / baseline['stream_file']
        golden_index = golden.with_suffix('.idx')
        entry = dict(dataset=row['dataset'], level=row['compression_level'], status='running',
                     workers=args.j, reference_bytes=row['compressed_bytes'],
                     vectors=row['vector_count'], ids=row['integer_count'])
        rows.append(entry)
        def phase(name, flags):
            """Record the current phase before launch so failures have a precise summary row."""
            entry['phase'] = name
            entry['log'] = str(Path(case) / (name + '.log'))
            save_summary(args.out, rows)
            text, result = run_phase(binary, flags, args.out, case, name, reserve, args.silent)
            entry.update(result)
            return text
        def reads(stream, index, prefix):
            """Compare with validated golden; full mode also checks raw IDs through serial traversal."""
            flags = ['-t', '-s', '-bvin', stream, '-bv-index-in', index, '-threads', str(args.j),
                     '-batch-vectors', '512', '-read-memory-mib', str(args.read_memory_mib),
                     *(['-read-mmap'] if args.read_mmap else []), *progress]
            # The golden stream itself already passed raw verification in this run.
            if stream != golden:
                text = phase(prefix + '-verify-threaded', ['-verify', '-golden-bvin', golden,
                             '-golden-index-in', golden_index, *flags])
                validate_read(text, row, 'verify', baseline['compressed_bytes'])
            if not args.fast:
                text = phase(prefix + '-verify-serial', ['-t', '-s', '-verify', '-u32in', source,
                                                       '-bvin', stream, *progress])
                validate_read(text, row, 'verify-serial', row['input_bytes'])
                text = phase(prefix + '-decode', ['-decode', *flags])
                validate_read(text, row, 'decode', row['input_bytes'])
        try:
            if row['dataset'] not in validated:
                golden_flags = ['-t', '-s', '-verify', '-u32in', source, '-bvin', golden,
                                '-bv-index-in', golden_index, '-threads', str(args.j),
                                '-batch-vectors', '512', '-read-memory-mib', str(args.read_memory_mib),
                                *(['-read-mmap'] if args.read_mmap else []), *progress]
                text = phase('golden-verify-threaded', golden_flags)
                validate_read(text, baseline, 'verify', baseline['input_bytes'])
                validated.add(row['dataset'])
            reads(original, original_index, 'reference')
            if not args.fast:
                check_space(args.out, reserve, row['input_bytes'])
                candidate = args.out / case / (case + '.bv')
                index = candidate.with_suffix('.idx')
                if candidate.exists() or index.exists():
                    raise TestFailure('Refusing to overwrite a candidate stream')
                text = phase('rebuild', ['-t', '-s', '-level', str(row['compression_level']),
                                        '-threads', str(args.j), '-batch-vectors', '512', '-u32in', source,
                                        '-bvout', candidate, '-bv-index-out', index, *progress])
                size = candidate.stat().st_size
                if number(text, 'Total vectors=') != row['vector_count'] or number(text, 'Total ints=') != row['integer_count']:
                    raise TestFailure('Rebuild scope mismatch')
                if number(text, 'Input consumed bytes=') != row['input_bytes'] or number(text, 'BV size = ') != size:
                    raise TestFailure('Rebuild byte accounting mismatch')
                validate_index(index, size, row['vector_count'])
                reads(candidate, index, 'candidate')
                entry.update(candidate_bytes=size, candidate_minus_reference_bytes=size-row['compressed_bytes'],
                             candidate_bits_per_id=8*size/row['integer_count'],
                             serialized_bytes_identical=(digest(candidate) == row['stream_sha256']))
            entry.update(status='passed', exit_code=0)
            save_summary(args.out, rows)
            if not args.silent:
                print(f'PASS {case}: logical equivalence against validated level {args.golden_level}', flush=True)
        except BaseException as error:
            entry.update(status='failed', exit_code=getattr(error, 'code', 130 if isinstance(error, KeyboardInterrupt) else 1))
            save_summary(args.out, rows)
            raise
    return rows


def parse():
    """Parse platform binary, immutable reference tree, execution budget, and test mode."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bin', required=True, help='Target executable path or executable name on PATH')
    parser.add_argument('--reference', required=True, type=Path)
    parser.add_argument('--out', required=True, type=Path, help='Fresh test directory outside reference tree')
    parser.add_argument('-j', '--threads', dest='j', type=int, default=max(1, (os.cpu_count() or 1)-2))
    parser.add_argument('--silent', action='store_true', help='Suppress console progress; retain full logs, final results and errors')
    parser.add_argument('--fast', action='store_true', help='Validate golden against raw IDs, then compare selected levels; no checksum audit or rebuild')
    parser.add_argument('--golden-level', type=int, choices=range(1, 11), default=6,
                        help='Golden compression level validated against raw IDs first (default 6)')
    parser.add_argument('--start-level', type=int, choices=range(1, 11), default=1,
                        help='First compression level to test, inclusive through level 10 (default 1)')
    parser.add_argument('--read-memory-mib', type=int, default=512, help='Combined staged input budget for indexed reads in positive MiB (default 512)')
    parser.add_argument('--read-mmap', action='store_true', help='Use Unix read-only BLOB mapping windows')
    parser.add_argument('--progress', type=int, default=1000)
    parser.add_argument('--reserve-gib', type=float, default=8)
    return parser.parse_args()


def main():
    """Audit reference inputs, snapshot the target, and run without writing in the reference tree."""
    args = parse()
    args.reference = args.reference.resolve(strict=True)
    args.out = args.out.resolve()
    protect_paths(args.reference, args.out)
    if args.out.exists():
        raise TestFailure('Require a fresh --out directory', 2)
    if not 1 <= args.read_memory_mib <= sys.maxsize // (1024 * 1024) or (args.read_mmap and os.name != 'posix'):
        raise TestFailure('Read budget must be positive and fit host address limits; mapping requires Unix', 2)
    if not 1 <= args.j <= 256 or args.progress <= 0 or args.reserve_gib < 0:
        raise TestFailure('Invalid threads, progress or reserve', 2)
    found = shutil.which(args.bin)
    if not found:
        raise TestFailure(f'Target executable not found: {args.bin}', 2)
    target = Path(found).resolve()
    if (args.reference / '.reference-lock').exists():
        raise TestFailure('Reference generator lock exists; finish or inspect that run first', 2)
    matrix = None if args.fast else json.loads((args.reference / 'matrix.json').read_text())
    table = list(csv.DictReader((args.reference / 'references.tsv').open(), delimiter='\t'))
    expected = {(f'Gov2.{dataset}', level) for dataset in ('sorted', 'unsorted') for level in range(1, 11)}
    if len(table) != 20 or {(row['dataset'], int(row['compression_level'])) for row in table} != expected:
        raise TestFailure('Reference matrix must contain exactly both datasets at levels 1–10', 2)
    table = [row for row in table if int(row['compression_level']) >= args.start_level
             or int(row['compression_level']) == args.golden_level]
    # Explicit local logs are created only after path checks. No reference writes occur.
    args.out.mkdir(parents=True)
    snapshots = {}
    references = []
    try:
        check_space(args.out, int(args.reserve_gib * 1024**3))
        if args.fast:
            if not args.silent:
                print('Fast mode: skipping reference checksum and historical validation audits', flush=True)
        else:
            for dataset in ('sorted', 'unsorted'):
                source = args.reference / 'Gov2Flat' / ('gov2.' + dataset)
                if not args.silent:
                    print(f'Checking copied raw input: {source}', flush=True)
                info = identity(source)
                if info['sha256'] != matrix['inputs'][dataset]['sha256']:
                    raise TestFailure(f'Copied raw input checksum mismatch: {dataset}', 2)
                snapshots[source] = info['sha256']
        for table_row in table:
            dataset = table_row['dataset'].split('.')[-1]
            expected_file = f'gov2.{dataset}.level-{int(table_row["compression_level"]):02d}.bv'
            if table_row['stream_file'] != expected_file:
                raise TestFailure('Unexpected reference filename in TSV', 2)
            key = Path(expected_file).stem
            if args.fast:
                row = dict(table_row)
                for field in ('compression_level', 'vector_count', 'integer_count', 'input_bytes', 'compressed_bytes'):
                    row[field] = int(row[field])
                required = [args.reference / expected_file,
                            (args.reference / expected_file).with_suffix('.idx'),
                            args.reference / 'Gov2Flat' / ('gov2.' + dataset)]
                if not all(path.is_file() for path in required):
                    raise TestFailure(f'Missing stream, index or raw input: {key}', 2)
                references.append(row)
            else:
                if not args.silent:
                    print(f'Checking reference integrity: {key}', flush=True)
                row = completed(args.reference, key, matrix['inputs'][dataset], matrix['build']['binary_sha256'])
                if row is None:
                    raise TestFailure(f'Missing committed reference: {key}', 2)
                for name, value in table_row.items():
                    if name not in row or str(row[name]) != value:
                        raise TestFailure(f'TSV/manifest mismatch: {key}, {name}', 2)
                references.append(row)
                for name, checksum in row['artifacts'].items():
                    snapshots[args.reference / name] = checksum
        if not args.fast:
            # Include manifests and auxiliary reference files in the immutability check.
            for path in args.reference.rglob('*'):
                if path.is_file() and path not in snapshots:
                    snapshots[path] = digest(path)
        # Snapshot target so a concurrent rebuild cannot change later phases.
        tools = args.out / '.tools'
        tools.mkdir()
        binary = tools / 'bminv'
        shutil.copy2(target, binary)
        shutil.copy2(HERE.parent / 'experiments/current/progress_console.py', tools / 'progress_console.py')
        write_json(args.out / 'test.json', dict(binary=str(target), binary_sha256=digest(binary),
                   reference=str(args.reference), mode='fast' if args.fast else 'normal', workers=args.j,
                   reference_integrity_audit=not args.fast, silent_progress=args.silent,
                   read_memory_mib=args.read_memory_mib, read_mmap=args.read_mmap,
                   start_level=args.start_level, golden_level=args.golden_level, end_level=10,
                   case_count=2 * (11 - args.start_level),
                   platform=platform.platform(), started_utc=datetime.now(timezone.utc).isoformat()))
        execute(args, references, binary)
        if not args.fast:
            if not args.silent:
                print('Rechecking reference integrity after target execution', flush=True)
            for path, checksum in snapshots.items():
                if digest(path) != checksum:
                    raise TestFailure(f'Reference changed during testing: {path}')
        case_count = 2 * (11 - args.start_level)
        message = (f'All {case_count} collection validations passed against a raw-validated golden; reference integrity audit omitted.\n' if args.fast else
                   f'All {case_count} reference cases passed; reference checksums unchanged.\n')
        (args.out / 'PASS.txt').write_text(message)
        print(f'ALL PASS: {case_count} cases; reports: {args.out / "summary.tsv"}', flush=True)
    except BaseException as error:
        (args.out / 'FAILURE.txt').write_text(str(error) + '\n')
        raise


def interrupted(signum, frame):
    """Route external termination through normal child cleanup."""
    raise TestFailure(f'Runner received {signal.Signals(signum).name}', 128 + signum)


if __name__ == '__main__':
    signal.signal(signal.SIGTERM, interrupted)
    try:
        main()
    except KeyboardInterrupt as error:
        print(f'FAILED: {error}', file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        print(f'FAILED: {error}', file=sys.stderr)
        sys.exit(getattr(error, 'code', 1))
