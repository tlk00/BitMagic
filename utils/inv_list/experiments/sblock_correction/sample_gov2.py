#!/usr/bin/env python3
"""Paired systematic midpoint sample of whole raw Gov2 vectors."""
import argparse
import hashlib
import json
import shutil
import struct
from pathlib import Path


def records(stream, size):
    offset = 0
    index = 0
    while offset < size:
        stream.seek(offset)
        header = stream.read(4)
        if len(header) != 4:
            raise ValueError(f'Truncated length at {offset}')
        count, = struct.unpack('<I', header)
        end = offset + 4 + 4 * count
        if not count or end > size:
            raise ValueError(f'Invalid vector at {offset}')
        yield index, offset, end, count
        index += 1
        offset = end


def sample(input_dir, output_dir, count):
    if count < 1:
        raise ValueError('Sample count must be positive')
    sources = {}
    for dataset in ('sorted', 'unsorted'):
        source = input_dir / ('gov2.' + dataset)
        stamp = source.stat()
        total = ids = 0
        with source.open('rb') as stream:
            for _, _, _, n in records(stream, stamp.st_size):
                total += 1
                ids += n
        sources[dataset] = (source, stamp, total, ids)
        print(f'Scanned {dataset}: {total:,} vectors, {ids:,} IDs', flush=True)
    totals = {(s[2], s[3]) for s in sources.values()}
    if len(totals) != 1:
        raise ValueError('Sorted and unsorted source counts differ')
    total, full_ids = totals.pop()
    if count > total:
        raise ValueError(f'Requested {count} vectors, collection has {total}')
    # One central record in each equal-width index stratum; no forced endpoints.
    indices = [(2*i+1)*total//(2*count) for i in range(count)]
    selected = set(indices)
    assert len(selected) == count
    output_dir.mkdir(parents=True, exist_ok=False)
    manifest = dict(method='systematic midpoint by vector index',
                    indices_zero_based=indices, sampled_vectors=count,
                    total_vectors=total, total_ids=full_ids, datasets={})
    lengths = None
    for dataset, (source, stamp, _, _) in sources.items():
        output = output_dir / ('gov2.' + dataset)
        entries = []
        digest = hashlib.sha256()
        with source.open('rb') as stream, output.open('xb') as dest:
            for index, offset, end, n in records(stream, stamp.st_size):
                if index not in selected:
                    continue
                remaining = end-offset
                if shutil.disk_usage(output_dir).free < remaining + 64*1024*1024:
                    raise OSError('Insufficient space for sampled raw vectors')
                stream.seek(offset)
                while remaining:
                    chunk = stream.read(min(remaining, 1024*1024))
                    if not chunk:
                        raise ValueError('Source truncated while copying')
                    dest.write(chunk)
                    digest.update(chunk)
                    remaining -= len(chunk)
                entries.append(dict(index=index, source_offset=offset, ids=n))
        after = source.stat()
        if (stamp.st_size, stamp.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError('Source changed during sampling')
        current_lengths = [r['ids'] for r in entries]
        if len(entries) != count or (lengths is not None and lengths != current_lengths):
            raise ValueError('Selected sorted/unsorted vector cardinalities differ')
        lengths = current_lengths
        manifest['datasets'][dataset] = dict(source=str(source.resolve()),
            source_bytes=stamp.st_size, source_mtime_ns=stamp.st_mtime_ns,
            bytes=output.stat().st_size, ids=sum(current_lengths),
            sha256=digest.hexdigest(), records=entries)
        print(f'Sampled {dataset}: {count:,} vectors, {sum(current_lengths):,} IDs; {output.stat().st_size:,} bytes', flush=True)
    (output_dir/'sample_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--count', type=int, default=10000)
    args = parser.parse_args()
    sample(args.input_dir, args.out, args.count)
