#!/usr/bin/env python3
"""Summarize measured GAP serialization routes from one verified collection."""
import argparse
import pathlib
import re

p = argparse.ArgumentParser()
p.add_argument('run_dir', type=pathlib.Path)
a = p.parse_args()
log = (a.run_dir / 'compress.log').read_text(errors='replace')
verify = (a.run_dir / 'verify.log').read_text(errors='replace')
decode = (a.run_dir / 'decode.log').read_text(errors='replace')

def fields(line):
    return {k: int(v) for k, v in re.findall(r'([a-z_]+)=(\d+)', line)}

rows = [x for x in log.splitlines() if x.startswith('GAP route ')]
if len(rows) != 11:
    raise SystemExit(f'Expected 11 GAP route rows, found {len(rows)}')
totals = fields(rows[0])
levels = [fields(x) for x in rows[1:5]]
buckets = [fields(x) for x in rows[5:]]
if [r.get('level') for r in levels] != list(range(4)):
    raise SystemExit('Malformed level rows')
if 'Validation: PASS' not in verify and 'Verification: PASS' not in verify:
    # inv_list versions differ in their final wording; the harness also checks
    # both process exit codes and full vector counts.
    if 'Total vectors=' not in verify:
        raise SystemExit('Verification summary missing')
if 'Total vectors=' not in decode:
    raise SystemExit('Decode summary missing')

def one(pattern):
    m = re.search(pattern, log)
    if not m:
        raise SystemExit(f'Missing compression statistic: {pattern}')
    return int(m.group(1))

vectors = one(r'Total vectors=(\d+)')
stream = one(r'BV size = (\d+)')
for phase, body in (('verify', verify), ('decode', decode)):
    m = re.search(r'Total vectors=(\d+)', body)
    if not m or int(m.group(1)) != vectors:
        raise SystemExit(f'{phase} vector count does not match compression')
direct = sum(r['direct_blocks'] for r in levels)
grouped = sum(r['grouped_blocks'] for r in levels)
direct_bytes = sum(r['direct_bytes'] for r in levels)
if direct + grouped == 0 or sum(r['direct_blocks'] for r in buckets) != levels[0]['direct_blocks']:
    raise SystemExit('GAP route counters failed reconciliation')
if sum(r['grouped_blocks'] for r in buckets) != levels[0]['grouped_blocks']:
    raise SystemExit('Grouped level-0 counters failed reconciliation')
if totals['superblock_bytes'] + direct_bytes + totals['direct_bit_bytes'] > stream:
    raise SystemExit('Record bytes exceed stream size')
residual = stream - totals['superblock_bytes'] - direct_bytes - totals['direct_bit_bytes']

print(f'Correctness: PASS; vectors={vectors:,}; stream_bytes={stream:,}')
print('Route                  Records/blocks      Emitted bytes    % of stream')
print(f"Sparse superblocks     {totals['superblocks']:>14,} {totals['superblock_bytes']:>18,} {100*totals['superblock_bytes']/stream:>13.2f}")
print(f'Individual GAP blocks  {direct:>14,} {direct_bytes:>18,} {100*direct_bytes/stream:>13.2f}')
print(f"Individual BIT blocks  {totals['direct_bit']:>14,} {totals['direct_bit_bytes']:>18,} {100*totals['direct_bit_bytes']/stream:>13.2f}")
print(f'Other framing/spans                 {residual:>18,} {100*residual/stream:>13.2f}')
print(f'Grouped GAP blocks: {grouped:,} ({100*grouped/(direct+grouped):.2f}% of all resident GAP blocks)')
print(f"Grouped BIT/full/empty: {totals['grouped_bit']:,}/{totals['grouped_full']:,}/{totals['grouped_empty']:,}")
print('GAP level     Direct blocks    Direct bytes   Grouped blocks    Grouped words')
for r in levels:
    print(f"{r['level']:>9} {r['direct_blocks']:>17,} {r['direct_bytes']:>15,} {r['grouped_blocks']:>16,} {r['grouped_words']:>16,}")
print('Level-0 GAP length (16-bit words):')
print('Length        Direct blocks    Direct bytes   Grouped blocks')
for label, r in zip(('1-4', '5-8', '9-16', '17-32', '33-64', '65+'), buckets):
    print(f"{label:>9} {r['direct_blocks']:>17,} {r['direct_bytes']:>15,} {r['grouped_blocks']:>16,}")
print('Superblock bytes include every constituent block; they cannot be attributed to one GAP block.')
