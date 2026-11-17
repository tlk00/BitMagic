#!/usr/bin/env python3
"""Build and validate immutable xsample05 compression references."""
import argparse
import concurrent.futures
import hashlib
import fcntl
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time


def digest(path, cancelled=lambda: False):
    """Hash a file in bounded memory."""
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            if cancelled():
                raise RuntimeError("Validation cancelled")
            h.update(chunk)
    return h.hexdigest()


def cases(levels):
    """Return stable case names and their input/mode settings."""
    for level in levels:
        for remap in (False, True):
            for xor in (False, True):
                yield (f'L{level:02d}-remap{int(remap)}-xor{int(xor)}', level, remap, xor)


def run(command, log):
    """Save output and fail on a nonzero child status."""
    with log.open('w') as output:
        result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f'Exit {result.returncode}: {log}')
    return log.read_text()


def check_case(case, root, inputs, cancelled=lambda: False):
    """Verify the published manifest and files before resume or validation."""
    name, level, remap, xor = case
    folder = root / 'cases' / name
    metadata = json.loads((folder / 'metadata.json').read_text())
    if metadata['case'] != name or metadata['level'] != level or metadata['remap'] != remap or metadata['xor'] != xor:
        raise RuntimeError(f'Case metadata mismatch: {name}')
    for filename in ('reference.sv', 'opcodes.txt', 'create.log', 'command.json', 'validation.log'):
        if digest(folder / filename, cancelled) != metadata['sha256'][filename]:
            raise RuntimeError(f'Case checksum mismatch: {name}/{filename}')
    if (folder / 'reference.sv').stat().st_size != metadata['bytes']:
        raise RuntimeError(f'Case size mismatch: {name}')
    if metadata['input_sha256'] != inputs['remap' if remap else 'plain']['sha256']:
        raise RuntimeError(f'Case input mismatch: {name}')
    return metadata


def create_case(case, args, inputs):
    """Create one case privately, validate it and publish by directory rename."""
    name, level, remap, xor = case
    target = args.reference_dir / 'cases' / name
    if target.exists():
        check_case(case, args.reference_dir, inputs)
        return name, 'existing'
    source = Path(inputs['remap' if remap else 'plain']['path'])
    work = Path(tempfile.mkdtemp(prefix=name + '-', dir=args.reference_dir / '.partial'))
    try:
        command = [str(args.bin), '-svin', str(source), '-svout', str(work / 'reference.sv'),
                   '-level', str(level), '-profile-opcodes']
        if xor:
            command.append('-xor')
        (work / 'command.json').write_text(json.dumps(command, indent=2) + '\n')
        start = time.monotonic()
        log = run(command, work / 'create.log')
        elapsed = time.monotonic() - start
        match = re.search(r'Serialization measurement bytes=(\d+) serialize_ms=([\d.]+) deserialize_ms=([\d.]+)', log)
        if not match or 'reconciled=1' not in log:
            raise RuntimeError(f'Missing measurement or opcode reconciliation: {name}')
        size = int(match.group(1))
        if not size or (work / 'reference.sv').stat().st_size != size:
            raise RuntimeError(f'Serialized file size mismatch: {name}')
        (work / 'opcodes.txt').write_text(log[log.index('Opcode profile blobs='):])
        run([str(args.bin), '-svin', str(work / 'reference.sv'), '-verify-sv', str(source)],
            work / 'validation.log')
        metadata = dict(case=name, level=level, remap=remap, xor=xor,
                        input_sha256=inputs['remap' if remap else 'plain']['sha256'],
                        bytes=size, serialize_ms=float(match.group(2)),
                        deserialize_ms=float(match.group(3)), wall_seconds=elapsed,
                        sha256={f: digest(work / f) for f in ('reference.sv', 'opcodes.txt', 'create.log', 'command.json', 'validation.log')})
        (work / 'metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
        # Existing references are never replaced, including when two writers race.
        if target.exists():
            check_case(case, args.reference_dir, inputs)
            shutil.rmtree(work)
            return name, 'existing'
        work.rename(target)
    except Exception:
        print(f'Incomplete case retained at {work}', file=sys.stderr, flush=True)
        raise
    return name, 'created'


def bundle_inputs(root, inputs):
    """Copy and verify prepared inputs under their portable names without rewriting metadata."""
    folder = root / 'inputs'
    folder.mkdir(exist_ok=True)
    for key, item in inputs.items():
        target = folder / (key + '.sv')
        if target.exists():
            if digest(target) != item['sha256']:
                raise RuntimeError(f'Bundled input checksum mismatch: {target}')
            continue
        source = Path(item['path'])
        temporary = target.with_suffix('.sv.tmp')
        shutil.copyfile(source, temporary)
        if digest(temporary) != item['sha256']:
            raise RuntimeError(f'Input checksum mismatch: {source}')
        temporary.replace(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('create', 'validate', 'bundle'))
    parser.add_argument('--reference-dir', type=Path, required=True)
    parser.add_argument('--bin', type=Path, help='Reader/creator executable; required except for bundle')
    parser.add_argument('--out', type=Path, default=Path('xsample05-validation-' + time.strftime('%Y%m%d-%H%M%S')), help='New output directory for validation logs, outside reference')
    parser.add_argument('--plain', type=Path, help='Prepared unmapped input; optional validation/bundle override')
    parser.add_argument('--remap', type=Path, help='Prepared remapped input; optional validation/bundle override')
    parser.add_argument('--levels', default='0-10', help='0-10 or comma-separated levels')
    parser.add_argument('-j', type=int, default=1, help='Maximum simultaneous xsample05 processes')
    args = parser.parse_args()
    if args.j < 1:
        parser.error('-j must be positive')
    levels = list(range(11)) if args.levels == '0-10' else sorted(set(map(int, args.levels.split(','))))
    if not levels or any(level < 0 or level > 10 for level in levels):
        parser.error('--levels must be within 0..10')
    if args.action != 'bundle':
        if args.bin is None:
            parser.error('--bin is required')
        args.bin = args.bin.resolve(strict=True)
        if not args.bin.is_file() or not os.access(args.bin, os.X_OK):
            parser.error('--bin must be executable')
    args.reference_dir = args.reference_dir.resolve()
    lock = None
    if args.action == 'create':
        args.reference_dir.mkdir(parents=True, exist_ok=True)
        lock = (args.reference_dir / '.create.lock').open('a+')
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Another reference creator is active')
        if not args.plain or not args.remap:
            parser.error('create requires --plain and --remap')
        args.reference_dir.mkdir(parents=True, exist_ok=True)
        inputs = {key: dict(path=str(path.resolve(strict=True)), sha256=digest(path.resolve(strict=True)))
                  for key, path in (('plain', args.plain), ('remap', args.remap))}
        matrix_path = args.reference_dir / 'matrix.json'
        matrix = dict(schema=1, inputs=inputs, binary_sha256=digest(args.bin))
        if matrix_path.exists():
            if json.loads(matrix_path.read_text()) != matrix:
                raise RuntimeError('Existing matrix input or executable differs')
        else:
            preflight = args.reference_dir / 'input-comparison.log'
            run([str(args.bin), '-svin', inputs['remap']['path'],
                 '-verify-sv', inputs['plain']['path']], preflight)
            with matrix_path.open('x') as output:
                json.dump(matrix, output, indent=2)
                output.write('\n')
        (args.reference_dir / 'cases').mkdir(exist_ok=True)
        (args.reference_dir / '.partial').mkdir(exist_ok=True)
        bundle_inputs(args.reference_dir, inputs)
        task = lambda case: create_case(case, args, inputs)
    else:
        matrix = json.loads((args.reference_dir / 'matrix.json').read_text())
        if matrix.get('schema') != 1 or set(matrix.get('inputs', {})) != {'plain', 'remap'}:
            raise RuntimeError('Invalid reference matrix schema or inputs: ' + str(args.reference_dir / 'matrix.json'))
        inputs = matrix['inputs']
        for key, replacement in (('plain', args.plain), ('remap', args.remap)):
            if replacement:
                inputs[key] = dict(path=str(replacement.resolve(strict=True)), sha256=inputs[key]['sha256'])
        for key, item in inputs.items():
            bundled = args.reference_dir / 'inputs' / (key + '.sv')
            explicit = args.plain if key == 'plain' else args.remap
            if explicit is None and bundled.is_file():
                item['path'] = str(bundled)
            elif not Path(item['path']).is_absolute():
                item['path'] = str(args.reference_dir / item['path'])
            if digest(Path(item['path'])) != item['sha256']:
                raise RuntimeError(f"Prepared input checksum mismatch: {key} | {item['path']}")
        if args.action == 'bundle':
            bundle_inputs(args.reference_dir, inputs)
            print('Portable inputs ready: ' + str(args.reference_dir / 'inputs'))
            return
        from reference_validation import validate
        validate(args, inputs, list(cases(levels)), check_case)
        return
    failed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.j) as pool:
        future_map = {pool.submit(task, case): case[0] for case in cases(levels)}
        for future in concurrent.futures.as_completed(future_map):
            name = future_map[future]
            try:
                print(*future.result(), flush=True)
            except Exception as error:
                failed.append(name)
                print(f'FAILED {name}: {error}', file=sys.stderr, flush=True)
    if failed:
        raise RuntimeError(f'{len(failed)} cases failed: {", ".join(sorted(failed))}')
    print(f'{args.action} complete: {4 * len(levels)} cases', flush=True)
    if lock:
        lock.close()


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, KeyError, TypeError) as error:
        print(f'Error: {error}', file=sys.stderr)
        sys.exit(1)
