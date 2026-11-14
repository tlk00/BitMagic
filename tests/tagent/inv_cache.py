"""Optional verified local copies of immutable inv_list reference trees."""
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
import time
import uuid


CHUNK_BYTES = 4 * 1024 * 1024


def inventory(root):
    """Return relative file sizes and mtimes; reject symlinks and locked trees."""
    if (root / '.reference-lock').exists():
        raise ValueError('reference generation lock exists')
    files = {}
    def fail_walk(error):
        """Propagate directory traversal failures rather than accepting a partial tree."""
        raise error
    for directory, dirs, names in os.walk(root, followlinks=False, onerror=fail_walk):
        base = Path(directory)
        for name in dirs:
            if (base / name).is_symlink():
                raise ValueError('reference contains symlink: ' + str(base / name))
        for name in names:
            path = base / name
            if path.is_symlink():
                raise ValueError('reference contains symlink: ' + str(path))
            if not path.is_file():
                raise ValueError('reference contains unsupported file: ' + str(path))
            info = path.stat()
            files[path.relative_to(root).as_posix()] = [info.st_size, info.st_mtime_ns]
    return files


def matches(root, files):
    """Check the cached file list and every length without rereading file contents."""
    actual = inventory(root)
    return {name: data[0] for name, data in actual.items()} == {
        name: data[0] for name, data in files.items()}


def remove_old_copies(parent, current, report):
    """Remove recognized published copies after current was verified and published.

    Leave unknown directories, symlinks, staging and quarantined files untouched.
    Log deletion failures without invalidating the newly published reference.
    """
    if parent.is_symlink():
        report('INV CACHE | cleanup skipped: cache parent is a symlink')
        return
    try:
        candidates = list(parent.iterdir())
    except OSError as error:
        report('INV CACHE | cleanup unavailable: ' + str(error))
        return
    for candidate in candidates:
        if candidate == current or not re.fullmatch(r'[0-9a-f]{24}', candidate.name):
            continue
        try:
            if candidate.is_symlink() or not candidate.is_dir():
                continue
            if {child.name for child in candidate.iterdir()} != {'reference', 'inventory.json'}:
                continue
            manifest = candidate / 'inventory.json'
            reference = candidate / 'reference'
            if manifest.is_symlink() or reference.is_symlink() or not reference.is_dir():
                continue
            saved = json.loads(manifest.read_text())
            if not isinstance(saved, dict) or saved.get('format') != 1 or not isinstance(saved.get('files'), dict):
                continue
            key = hashlib.sha256(json.dumps(saved, sort_keys=True).encode()).hexdigest()[:24]
            if key != candidate.name:
                continue
            size = sum(item[0] for item in saved['files'].values())
            report('INV CACHE | removing old copy {} | {} bytes recorded'.format(candidate, size))
            started = time.monotonic()
            shutil.rmtree(candidate)
            report('INV CACHE | old copy removed {} | {:.1f}s'.format(candidate, time.monotonic() - started))
        except (OSError, ValueError, TypeError, IndexError, json.JSONDecodeError) as error:
            report('INV CACHE | cannot remove old copy {}: {}'.format(candidate, error))


def copy_reference(source, cache_root, reserve_bytes, report, progress=None):
    """Return a validated cached reference, or raise on failure for caller fallback.

    source must be immutable during copying. cache_root must already exist and be
    writable. Each file is SHA-256 verified by rereading its completed copy; only
    an entirely verified tree with unchanged source metadata is published.
    Existing copies use their saved source inventory and per-file sizes.
    report receives operator-readable status and per-file timing messages.
    progress, when supplied, receives phase, aggregate bytes copied, total bytes
    and the current filename. It must not delay the copy.
    """
    source = source.resolve()
    cache_root = cache_root.resolve()
    if not cache_root.is_dir():
        raise ValueError('cache directory does not exist: ' + str(cache_root))
    if cache_root == source or cache_root.is_relative_to(source) or source.is_relative_to(cache_root):
        raise ValueError('cache and reference directories must not contain one another')
    with tempfile.TemporaryDirectory(prefix='.tagent-probe-', dir=cache_root) as probe:
        (Path(probe) / 'write-test').write_bytes(b'TAgent')
    started = time.monotonic()
    files = inventory(source)
    if not files:
        raise ValueError('reference directory is empty')
    source_info = {'format': 1, 'source': str(source), 'files': files}
    key = hashlib.sha256(json.dumps(source_info, sort_keys=True).encode()).hexdigest()[:24]
    parent = cache_root / 'inv-reference-v1'
    parent.mkdir(exist_ok=True)
    destination = parent / key
    reference = destination / 'reference'
    total = sum(item[0] for item in files.values())
    report('INV CACHE | source {} | {} files | {} bytes'.format(source, len(files), total))
    if destination.exists():
        saved = None
        try:
            saved = json.loads((destination / 'inventory.json').read_text())
            if saved == source_info and matches(reference, files):
                report('INV CACHE | reuse {} | {} bytes | {:.1f}s'.format(reference, total, time.monotonic() - started))
                return reference
        except (OSError, ValueError, json.JSONDecodeError):
            pass
        if saved != source_info:
            raise ValueError('unrecognized existing cache directory: ' + str(destination))
        quarantine = parent / ('.invalid-' + key + '-' + uuid.uuid4().hex)
        os.rename(destination, quarantine)
        report('INV CACHE | altered cache preserved at ' + str(quarantine))
    free = shutil.disk_usage(parent).free
    if free < total + reserve_bytes:
        raise ValueError('insufficient cache space: {} bytes free; {} bytes required including reserve'.format(free, total + reserve_bytes))
    report('INV CACHE | copying to {} | free {} bytes | reserve {} bytes'.format(cache_root, free, reserve_bytes))
    completed_bytes = 0
    if progress:
        progress('copying', 0, total, '')
    with tempfile.TemporaryDirectory(prefix='.copy-', dir=parent) as temp:
        staging = Path(temp)
        tree = staging / 'reference'
        tree.mkdir()
        for name, (size, mtime) in sorted(files.items()):
            before = time.monotonic()
            report('INV CACHE | copying {} | {} bytes'.format(name, size))
            progress_at = before + 30
            copied_bytes = 0
            source_file = source / name
            target = tree / name
            target.parent.mkdir(parents=True, exist_ok=True)
            checksum = hashlib.sha256()
            with source_file.open('rb') as incoming, target.open('xb') as output:
                while True:
                    data = incoming.read(CHUNK_BYTES)
                    if not data:
                        break
                    checksum.update(data)
                    output.write(data)
                    copied_bytes += len(data)
                    if time.monotonic() >= progress_at:
                        report('INV CACHE | copying {} | {} of {} bytes | {:.1f}s'.format(
                            name, copied_bytes, size, time.monotonic() - before))
                        progress_at = time.monotonic() + 30
                        if progress:
                            progress('copying', completed_bytes + copied_bytes, total, name)
            copied = time.monotonic()
            completed_bytes += copied_bytes
            report('INV CACHE | verifying {} | {} bytes'.format(name, size))
            if progress:
                progress('verifying', completed_bytes, total, name)
            actual = hashlib.sha256()
            checked_bytes = 0
            progress_at = copied + 30
            with target.open('rb') as stream:
                for data in iter(lambda: stream.read(CHUNK_BYTES), b''):
                    actual.update(data)
                    checked_bytes += len(data)
                    if time.monotonic() >= progress_at:
                        report('INV CACHE | verifying {} | {} of {} bytes | {:.1f}s'.format(
                            name, checked_bytes, size, time.monotonic() - copied))
                        progress_at = time.monotonic() + 30
                        if progress:
                            progress('verifying', completed_bytes, total, name)
            current = source_file.stat()
            if (current.st_size != size or current.st_mtime_ns != mtime or
                    target.stat().st_size != size or actual.digest() != checksum.digest()):
                raise ValueError('copy verification or source stability failed: ' + name)
            report('INV CACHE | file {} | {} bytes | copy {:.1f}s | verify {:.1f}s'.format(
                name, size, copied - before, time.monotonic() - copied))
        if inventory(source) != files:
            raise ValueError('reference changed during copying')
        (staging / 'inventory.json').write_text(json.dumps(source_info, sort_keys=True))
        os.rename(staging, destination)
    report('INV CACHE | ready {} | {} bytes | {:.1f}s'.format(reference, total, time.monotonic() - started))
    if progress:
        progress('ready', total, total, '')
    remove_old_copies(parent, destination, report)
    return reference
