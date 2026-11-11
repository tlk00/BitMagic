#!/usr/bin/env python3
"""Validate a target bminv against both Gov2 datasets and saved levels 1–10."""
import argparse
import csv
import hashlib
import uuid
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
          'serialized_bytes_identical', 'roundtrip_log', 'vectors', 'ids', 'workers', 'wall_seconds', 'log']


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


def ensure_hash_cache(args, baseline, source, golden, golden_index, phase):
    """Validate raw equality before first cache build; reuse only matching provenance.

    A per-cache mkdir lock serializes publication. Temporary outputs have unique
    names and are cleaned on failure. The ready cache and JSON are immutable.
    Reuse is keyed by frozen TSV hashes and scope, making cache files portable.
    File metadata detects mutation during the first validation/build only.
    """
    def stamp(path):
        """Return cheap source identity without reading its data bytes."""
        info = path.stat()
        return dict(path=str(path), bytes=info.st_size, mtime_ns=info.st_mtime_ns)
    expected = dict(format='BMBLOCKSHA256V1', dataset=baseline['dataset'],
                    level=args.golden_level, vectors=baseline['vector_count'],
                    ids=baseline['integer_count'], golden_sha256=baseline['stream_sha256'],
                    raw_sha256=baseline['input_sha256'], raw_bytes=baseline['input_bytes'])
    key = hashlib.sha256(json.dumps(expected, sort_keys=True).encode()).hexdigest()[:20]
    args.hash_cache_dir.mkdir(parents=True, exist_ok=True)
    cache = args.hash_cache_dir / (Path(baseline['stream_file']).stem + '.' + key + '.sha256')
    manifest = cache.with_suffix(cache.suffix + '.json')
    lock = cache.with_suffix(cache.suffix + '.lock')
    try:
        lock.mkdir()
    except FileExistsError:
        raise TestFailure(f'Hash cache is being prepared, or a previous run left a lock: {lock}', 2)
    staging = cache.with_name(cache.name + '.building-' + uuid.uuid4().hex)
    try:
        if cache.exists() or manifest.exists():
            if not cache.is_file() or not manifest.is_file():
                raise TestFailure(f'Incomplete hash cache publication: {cache}', 2)
            saved = json.loads(manifest.read_text())
            if saved.get('provenance') != expected or saved.get('cache_sha256') != digest(cache):
                raise TestFailure(f'Hash cache provenance/integrity mismatch: {cache}', 2)
            if not args.silent:
                print(f'Reusing SHA-256 cache: {cache}', flush=True)
            return cache
        origin = dict(raw=stamp(source), golden=stamp(golden), index=stamp(golden_index))
        common = ['-t', '-s', '-bvin', golden, '-bv-index-in', golden_index,
                  '-threads', str(args.j), '-batch-vectors', '512',
                  '-read-memory-mib', str(args.read_memory_mib),
                  *(['-read-mmap'] if args.read_mmap else []), '-progress-vectors', str(args.progress)]
        text = phase('golden-verify-threaded', ['-verify', '-u32in', source, *common])
        validate_read(text, baseline, 'verify', baseline['input_bytes'])
        text = phase('hash-cache-build', ['-hash-cache-out', staging, *common])
        if ('Verification complete.' not in text or
            number(text, 'Total vectors=') != baseline['vector_count'] or
            number(text, 'Total ints=') != baseline['integer_count']):
            raise TestFailure('Incomplete hash cache generation')
        with staging.open('rb') as stream:
            header = stream.read(64)
        if (header[:8] != b'BMSHA001' or len(header) != 64 or
            int.from_bytes(header[8:16], 'little') != baseline['vector_count'] or
            header[32:64].hex() != baseline['stream_sha256']):
            raise TestFailure('Hash cache source checksum differs from golden reference')
        if origin != dict(raw=stamp(source), golden=stamp(golden), index=stamp(golden_index)):
            raise TestFailure('Hash cache source changed during validation/build')
        saved = dict(provenance=expected, origin=origin, cache_sha256=digest(staging))
        staging.rename(cache)
        try:
            write_json(manifest, saved)
        except BaseException:
            cache.unlink()
            temporary_manifest = manifest.with_name(manifest.name + '.tmp')
            if temporary_manifest.exists():
                temporary_manifest.unlink()
            raise
        return cache
    finally:
        for temporary in (staging, Path(str(staging) + '.partial')):
            if temporary.exists():
                temporary.unlink()
        staging_lock = Path(str(staging) + '.lock')
        if staging_lock.is_dir():
            staging_lock.rmdir()
        lock.rmdir()


def execute(args, references, binary):
    """Validate each golden against raw IDs before comparing selected reference/candidate streams."""
    rows = []
    reserve = int(args.reserve_gib * 1024**3)
    progress = ['-progress-vectors', str(args.progress)]
    baselines = {row['dataset']: row for row in references
                 if row['compression_level'] == args.golden_level}
    validated = set()
    hash_caches = {}
    roundtrips = {}
    for row in references:
        if row['compression_level'] < args.start_level:
            continue
        dataset = row['dataset'].split('.')[-1]
        case = Path(row['stream_file']).stem
        source = args.reference / 'Gov2Flat' / ('gov2.' + dataset)
        original = args.reference / row['stream_file']
        original_index = original.with_suffix('.idx')
        baseline = baselines[row['dataset']]
        if any(row[field] != baseline[field] for field in
               ('vector_count', 'integer_count', 'input_bytes', 'input_sha256')):
            raise TestFailure('Selected levels have inconsistent dataset scope', 2)
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
            if args.hash:
                text = phase(prefix + '-verify-hash', ['-verify', '-hash-cache-in',
                             hash_caches[row['dataset']], *flags])
                if ('Verification complete.' not in text or
                    number(text, 'Total vectors=') != row['vector_count'] or
                    number(text, 'Total ints=') != row['integer_count']):
                    raise TestFailure('Incomplete cached SHA-256 verification')
            elif stream != golden:
                text = phase(prefix + '-verify-threaded', ['-verify', '-golden-bvin', golden,
                             '-golden-index-in', golden_index, *flags])
                validate_read(text, row, 'verify', baseline['compressed_bytes'])
        try:
            if args.hash and row['dataset'] not in hash_caches:
                hash_caches[row['dataset']] = ensure_hash_cache(
                    args, baseline, source, golden, golden_index, phase)
            if not args.hash and row['dataset'] not in validated:
                golden_flags = ['-t', '-s', '-verify', '-u32in', source, '-bvin', golden,
                                '-bv-index-in', golden_index, '-threads', str(args.j),
                                '-batch-vectors', '512', '-read-memory-mib', str(args.read_memory_mib),
                                *(['-read-mmap'] if args.read_mmap else []), *progress]
                text = phase('golden-verify-threaded', golden_flags)
                validate_read(text, baseline, 'verify', baseline['input_bytes'])
                validated.add(row['dataset'])
            if not args.fast and row['dataset'] not in roundtrips:
                flags = ['-t', '-s', '-verify', '-roundtrip', '-bvin', golden,
                         '-bv-index-in', golden_index, '-threads', str(args.j),
                         '-batch-vectors', '512', '-read-memory-mib', str(args.read_memory_mib),
                         '-roundtrip-start-level', str(args.start_level), '-roundtrip-end-level', '10',
                         *(['-read-mmap'] if args.read_mmap else []), *progress]
                if args.hash:
                    flags += ['-hash-cache-in', hash_caches[row['dataset']]]
                text = phase('roundtrip-levels', flags)
                if ('Roundtrip complete.' not in text or
                    number(text, 'Total vectors=') != baseline['vector_count'] or
                    number(text, 'Total ints=') != baseline['integer_count'] or
                    number(text, 'Total round trips=') != baseline['vector_count'] * (11 - args.start_level)):
                    raise TestFailure('Incomplete in-memory round-trip verification')
                levels = {}
                for line in text.splitlines():
                    if line.startswith('Roundtrip level='):
                        values = dict(field.split('=', 1) for field in line.split()[1:])
                        level = int(values['level'])
                        if level in levels or int(values['vectors']) != baseline['vector_count']:
                            raise TestFailure('Duplicate/incomplete round-trip level result')
                        levels[level] = values
                if set(levels) != set(range(args.start_level, 11)):
                    raise TestFailure('Missing in-memory round-trip levels')
                roundtrips[row['dataset']] = (levels, entry['log'])
            if not args.fast:
                levels, roundtrip_log = roundtrips[row['dataset']]
                rebuilt_bytes = int(levels[row['compression_level']]['bytes'])
                entry.update(candidate_bytes=rebuilt_bytes,
                             candidate_minus_reference_bytes=rebuilt_bytes-row['compressed_bytes'],
                             candidate_bits_per_id=8*rebuilt_bytes/row['integer_count'],
                             roundtrip_log=roundtrip_log)
                if not args.pedantic:
                    entry.update(phase='roundtrip-levels', log=roundtrip_log)
            if args.fast or args.pedantic:
                reads(original, original_index, 'reference')
                if args.pedantic:
                    text = phase('reference-verify-serial', ['-t', '-s', '-verify', '-u32in', source,
                                 '-bvin', original, *progress])
                    validate_read(text, row, 'verify-serial', row['input_bytes'])
                    text = phase('reference-decode', ['-decode', '-bvin', original,
                                 '-bv-index-in', original_index, '-threads', str(args.j),
                                 '-read-memory-mib', str(args.read_memory_mib),
                                 *(['-read-mmap'] if args.read_mmap else []), *progress])
                    validate_read(text, row, 'decode', row['input_bytes'])
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
    parser.add_argument('--pedantic', action='store_true', help='Add saved-level decoding/verification and full reference integrity audits to normal round trips')
    parser.add_argument('--fast', action='store_true', help='Validate golden against raw IDs, then compare selected levels; no checksum audit or rebuild')
    parser.add_argument('--hash', action='store_true', help='Opt in to cached canonical-block SHA-256 validation')
    parser.add_argument('--hash-cache-dir', type=Path, help='Persistent cache directory (default <reference>/validation-cache)')
    parser.add_argument('--dataset', choices=('both', 'sorted', 'unsorted'), default='both',
                        help='Select Gov2 input collection (default: both)')
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
    if args.fast and args.pedantic:
        raise TestFailure('--fast and --pedantic are separate validation modes', 2)
    args.reference = args.reference.resolve(strict=True)
    args.out = args.out.resolve()
    protect_paths(args.reference, args.out)
    if args.hash_cache_dir and not args.hash:
        raise TestFailure('--hash-cache-dir requires --hash', 2)
    args.hash_cache_dir = (args.hash_cache_dir or args.reference / 'validation-cache').resolve()
    if args.hash and (args.hash_cache_dir == args.reference or
                      args.hash_cache_dir.is_relative_to(args.out) or
                      args.out.is_relative_to(args.hash_cache_dir)):
        raise TestFailure('Hash cache must be a separate persistent subdirectory disjoint from --out', 2)
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
    matrix = None if not args.pedantic else json.loads((args.reference / 'matrix.json').read_text())
    table = list(csv.DictReader((args.reference / 'references.tsv').open(), delimiter='\t'))
    expected = {(f'Gov2.{dataset}', level) for dataset in ('sorted', 'unsorted') for level in range(1, 11)}
    if len(table) != 20 or {(row['dataset'], int(row['compression_level'])) for row in table} != expected:
        raise TestFailure('Reference matrix must contain exactly both datasets at levels 1–10', 2)
    datasets = ('sorted', 'unsorted') if args.dataset == 'both' else (args.dataset,)
    table = [row for row in table if row['dataset'].split('.')[-1] in datasets and
             (int(row['compression_level']) >= args.start_level
              or int(row['compression_level']) == args.golden_level)]
    # Explicit local logs are created only after path checks. No reference writes occur.
    args.out.mkdir(parents=True)
    snapshots = {}
    references = []
    try:
        check_space(args.out, int(args.reserve_gib * 1024**3))
        if not args.pedantic:
            if not args.silent:
                print('Using reference metadata without full checksum/history audits', flush=True)
        else:
            for dataset in datasets:
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
            if not args.pedantic:
                row = dict(table_row)
                for field in ('compression_level', 'vector_count', 'integer_count', 'input_bytes', 'compressed_bytes'):
                    row[field] = int(row[field])
                required = []
                if args.fast or row['compression_level'] == args.golden_level:
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
        if args.pedantic:
            # Include manifests and auxiliary reference files in the immutability check.
            for path in args.reference.rglob('*'):
                if (path.is_file() and path not in snapshots and
                    not (args.hash and path.is_relative_to(args.hash_cache_dir))):
                    snapshots[path] = digest(path)
        # Snapshot target so a concurrent rebuild cannot change later phases.
        tools = args.out / '.tools'
        tools.mkdir()
        binary = tools / 'bminv'
        shutil.copy2(target, binary)
        shutil.copy2(HERE.parent / 'experiments/current/progress_console.py', tools / 'progress_console.py')
        write_json(args.out / 'test.json', dict(binary=str(target), binary_sha256=digest(binary),
                   reference=str(args.reference), mode='fast' if args.fast else ('pedantic' if args.pedantic else 'normal'), workers=args.j,
                   reference_integrity_audit=args.pedantic, silent_progress=args.silent,
                   read_memory_mib=args.read_memory_mib, read_mmap=args.read_mmap,
                   start_level=args.start_level, golden_level=args.golden_level, hash_validation=args.hash,
                   hash_cache_dir=str(args.hash_cache_dir) if args.hash else None, end_level=10,
                   case_count=len(datasets) * (11 - args.start_level), datasets=list(datasets),
                   platform=platform.platform(), started_utc=datetime.now(timezone.utc).isoformat()))
        if args.hash or not args.fast:
            support = subprocess.run([str(binary), '--help'], capture_output=True, text=True, timeout=30)
            (args.out / 'hash-support.log').write_text(support.stdout + support.stderr)
            if support.returncode or (not args.fast and '-roundtrip ' not in support.stdout):
                raise TestFailure('Target binary lacks in-memory round-trip support; rebuild inv_list', 2)
            if args.hash and ('-hash-cache-in' not in support.stdout or
                              'SHA-256 cache support disabled' in support.stdout):
                raise TestFailure('Target binary lacks SHA-256 cache support; rebuild with crypto enabled', 2)
        execute(args, references, binary)
        if args.pedantic:
            if not args.silent:
                print('Rechecking reference integrity after target execution', flush=True)
            for path, checksum in snapshots.items():
                if digest(path) != checksum:
                    raise TestFailure(f'Reference changed during testing: {path}')
        case_count = len(datasets) * (11 - args.start_level)
        message = (f'All {case_count} collection validations passed against a raw-validated golden; reference integrity audit omitted.\n' if args.fast else
                   f'All {case_count} reference cases passed; reference checksums unchanged.\n' if args.pedantic else
                   f'All {case_count} in-memory round-trip cases passed.\n')
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
