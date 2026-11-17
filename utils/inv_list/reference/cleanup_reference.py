#!/usr/bin/env python3
"""Remove generated AppleDouble metadata and empty reference staging directories.

Never follows directory symlinks, deletes ordinary dotfiles, or removes frozen
tools and failed/incomplete work. Payload files and their checksums are unchanged.
"""
import argparse
import errno
import os
from pathlib import Path
import struct


def is_appledouble(path):
    """Recognize a regular AppleDouble v2 file by its header and bounded entry table.

    path is a candidate ._ file. Symlinks, malformed files and other formats are
    rejected. Reads only metadata headers; returns True for a valid container.
    """
    if path.is_symlink() or not path.is_file():
        return False
    size = path.stat().st_size
    with path.open('rb') as stream:
        header = stream.read(26)
        if len(header) != 26 or header[:8] != b'\x00\x05\x16\x07\x00\x02\x00\x00':
            return False
        count = struct.unpack('>H', header[24:26])[0]
        table_end = 26 + count * 12
        if not count or table_end > size:
            return False
        for _ in range(count):
            entry = stream.read(12)
            if len(entry) != 12:
                return False
            _, offset, length = struct.unpack('>III', entry)
            if offset < table_end or offset > size or length > size - offset:
                return False
    return True


def cleanup_output(root):
    """Clean generated metadata below root, returning file count, bytes and dirs.

    The caller must own the output directory or hold its reference lock. Only
    validated ._ metadata files and the empty top-level .partial directory are
    removed. Failed work, tools, user dotfiles and symlink targets are preserved.
    Filesystem errors propagate so callers can report incomplete cleanup.
    """
    root = Path(root)
    removed = 0
    total_bytes = 0
    for directory, dirs, files in os.walk(root, followlinks=False):
        parent = Path(directory)
        dirs[:] = [name for name in dirs if not (parent / name).is_symlink()]
        for name in files:
            if not name.startswith('._'):
                continue
            path = parent / name
            if is_appledouble(path):
                total_bytes += path.stat().st_size
                path.unlink()
                removed += 1
    empty_dirs = 0
    partial = root / '.partial'
    if partial.is_dir() and not partial.is_symlink():
        try:
            partial.rmdir()
            empty_dirs = 1
        except OSError as error:
            if error.errno not in (errno.ENOTEMPTY, errno.EEXIST):
                raise
    return removed, total_bytes, empty_dirs


def main():
    """Clean one existing matrix under its exclusive lock; print removed totals."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    args = parser.parse_args()
    root = args.reference.resolve(strict=True)
    if not (root / 'matrix.json').is_file():
        parser.error('Expected an existing reference directory with matrix.json')
    lock = root / '.reference-lock'
    lock.mkdir()  # Never steal an active or stale reference lock.
    try:
        files, size, directories = cleanup_output(root)
    finally:
        cleanup_output(lock)
        lock.rmdir()
        # Removing a directory can leave its own AppleDouble companion behind.
        cleanup_output(root)
    print(f'Removed {files} metadata files ({size} bytes) and {directories} empty staging directories')


if __name__ == '__main__':
    main()
