#!/usr/bin/env python3
"""Run exact/SHA pair checks, build raw/BV caches, then scan levels without golden reads."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
import platform
from pathlib import Path
import subprocess
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--bin', required=True, type=Path)
p.add_argument('--reference', required=True, type=Path)
p.add_argument('--out', required=True, type=Path)
p.add_argument('--dataset', choices=['sorted', 'unsorted'], default='sorted')
p.add_argument('--threads', type=int, default=8)
p.add_argument('--memory-mib', type=int, default=1024)
p.add_argument('--batch-vectors', type=int, default=512)
p.add_argument('--max-vectors', type=int, default=0, help='Explicit prefix experiment; default whole collection')
p.add_argument('--mmap', action='store_true')
a = p.parse_args()
a.out = a.out.resolve()
a.reference = a.reference.resolve(strict=True)
if a.out == a.reference or a.reference in a.out.parents:
    p.error('--out must be outside reference tree')
a.out.mkdir(parents=True, exist_ok=False)
binary = a.bin.resolve(strict=True)
common = ['--threads', str(a.threads), '--memory-mib', str(a.memory_mib), '--batch-vectors', str(a.batch_vectors)]
if a.max_vectors:
    common += ['--max-vectors', str(a.max_vectors)]
if a.mmap:
    common.append('--mmap')
results = []
(a.out / 'run.json').write_text(json.dumps(dict(binary=str(binary), binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
    reference=str(a.reference), dataset=a.dataset, threads=a.threads, memory_mib=a.memory_mib,
    max_vectors=a.max_vectors, mapped=a.mmap, batch_vectors=a.batch_vectors,
    platform=platform.platform(), logical_cpus=os.cpu_count(), started=datetime.now(timezone.utc).isoformat()), indent=2))


def collection(level, role):
    """Construct read-only collection and independent-index arguments."""
    path = a.reference / f'gov2.{a.dataset}.level-{level:02d}.bv'
    return ['--' + role, str(path), '--' + role + '-index', str(path.with_suffix('.idx'))]


def phase(name, mode, flags):
    """Run one measured phase and save its exact command, output, and result."""
    command = [str(binary), '--mode', mode, *flags, *common]
    (a.out / (name + '.command.json')).write_text(json.dumps(command, indent=2))
    print(name, flush=True)
    start = time.monotonic()
    with (a.out / (name + '.log')).open('w') as log:
        result = subprocess.run(command, stdout=log, stderr=log)
    text = (a.out / (name + '.log')).read_text()
    row = dict(phase=name, wall_seconds=time.monotonic()-start, exit_code=result.returncode)
    passes = [line for line in text.splitlines() if line.startswith('PASS ')]
    if passes:
        row.update(dict(field.split('=', 1) for field in passes[-1].split()[1:]))
    results.append(row)
    (a.out / 'results.json').write_text(json.dumps(results, indent=2))
    if result.returncode or not passes:
        (a.out / 'FAILURE.txt').write_text(name + '\n' + text[-4000:])
        raise SystemExit(f'FAILED {name}, inspect {a.out / (name + ".log")}')
    print(passes[-1], flush=True)


phase('exact-6-vs-10', 'exact', collection(10, 'candidate') + collection(6, 'golden'))
phase('sha-6-vs-10', 'pair-sha', collection(10, 'candidate') + collection(6, 'golden'))
bv_cache = a.out / 'level6.sha'
raw_cache = a.out / 'raw.sha'
phase('build-bv-cache', 'cache-bv', collection(6, 'candidate') + ['--cache', str(bv_cache)])
phase('build-raw-cache', 'cache-raw', ['--raw', str(a.reference / 'Gov2Flat' / ('gov2.' + a.dataset)), '--cache', str(raw_cache)])
# Compare record streams with constant memory; source-provenance headers differ.
with bv_cache.open('rb') as left, raw_cache.open('rb') as right:
    left.seek(8)
    right.seek(8)
    count = int.from_bytes(left.read(8), 'little')
    if right.read(8) != count.to_bytes(8, 'little'):
        (a.out / 'FAILURE.txt').write_text('Raw/BV cache vector counts differ\n')
        raise SystemExit('Raw/BV cache vector counts differ')
    left.seek(64)
    right.seek(64)
    remaining = count * 40
    while remaining:
        take = min(remaining, 65536)
        if left.read(take) != right.read(take):
            (a.out / 'FAILURE.txt').write_text('Raw and BV digest records differ\n')
            raise SystemExit('Raw and BV digest records differ')
        remaining -= take
for level in [7, 8, 9, 10]:
    phase(f'cached-level-{level}', 'cached', collection(level, 'candidate') + ['--cache', str(raw_cache)])
(a.out / 'PASS.txt').write_text('Exact/SHA comparisons, matching raw/BV caches and cached candidate levels passed.\n')
