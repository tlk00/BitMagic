#!/usr/bin/env python3
"""Validate and summarize the indexed BIC superblock length decode survey."""
import argparse
import csv
from pathlib import Path
import re

parser = argparse.ArgumentParser()
parser.add_argument('run_dir', type=Path)
args = parser.parse_args()
log = (args.run_dir / 'sblock-length.decode.log').read_text(errors='replace')
route = (args.run_dir / 'report.txt').read_text(errors='replace')

def kv(line):
    return {key: value for key, value in re.findall(r'([a-z0-9_]+)=([0-9.]+)', line)}

summary_lines = [line for line in log.splitlines() if line.startswith('SBlock BIC summary ')]
bin_lines = [line for line in log.splitlines() if line.startswith('SBlock BIC bin=')]
if len(summary_lines) != 1 or len(bin_lines) != 64 or 'Decode complete.' not in log:
    raise SystemExit('Incomplete BIC length decode survey')
summary = kv(summary_lines[0])
rows = [kv(line) for line in bin_lines]
count = int(summary['records'])
encoded_bytes = int(summary['bytes'])
if sum(int(r['records']) for r in rows) != count or sum(int(r['bytes']) for r in rows) != encoded_bytes:
    raise SystemExit('Histogram does not sum to the decoder totals')
for i, r in enumerate(rows):
    if int(r['bin']) != i or int(r['lo']) > int(r['hi']):
        raise SystemExit('Invalid histogram bin')
    if i and int(rows[i-1]['hi']) + 1 != int(r['lo']):
        raise SystemExit('Histogram intervals are not contiguous')
if int(rows[0]['lo']) != 1 or int(rows[-1]['hi']) != 65535:
    raise SystemExit('Histogram does not cover valid BIC lengths')
route_line = next((line for line in route.splitlines() if line.startswith('Sparse superblocks ')), None)
if not route_line:
    raise SystemExit('GAP route totals missing')
route_match = re.search(r'Sparse superblocks\s+([\d,]+)\s+([\d,]+)', route_line)
if not route_match:
    raise SystemExit('Malformed GAP route totals')
if count != int(route_match.group(1).replace(',', '')) or encoded_bytes != int(route_match.group(2).replace(',', '')):
    raise SystemExit('Decode histogram does not match previously measured emitted records')

csv_path = args.run_dir / 'sblock-length.csv'
with csv_path.open('w', newline='') as output:
    writer = csv.writer(output)
    writer.writerow(('bin', 'min_length', 'max_length', 'records', 'record_pct',
                     'encoded_bytes', 'byte_pct', 'avg_bytes_per_record'))
    for r in rows:
        n = int(r['records'])
        b = int(r['bytes'])
        writer.writerow((r['bin'], r['lo'], r['hi'], n, f'{100*n/count:.6f}',
                         b, f'{100*b/encoded_bytes:.6f}', f'{b/n:.3f}' if n else ''))

print(f"BIC superblocks: {count:,}; encoded bytes: {encoded_bytes:,}")
print(f"Length (set positions): min {summary['min']}, max {summary['max']}, "
      f"mean {summary['mean']}, exact p75 {summary['p75']}")
print('The 64 intervals cover lengths 1..65535; lengths 1..16 have individual bins,')
print('and higher lengths use four intervals per doubling.')
print(f'Full 64-bin histogram: {csv_path}')
print('Range       Records     % records    Encoded bytes     % bytes')
for lo, hi in ((1, 16), (17, 255), (256, 1023), (1024, 4095),
               (4096, 16383), (16384, 65535)):
    selected = [r for r in rows if int(r['lo']) >= lo and int(r['hi']) <= hi]
    n = sum(int(r['records']) for r in selected)
    b = sum(int(r['bytes']) for r in selected)
    print(f'{lo:>5}-{hi:<5} {n:>12,} {100*n/count:>11.2f} '
          f'{b:>16,} {100*b/encoded_bytes:>10.2f}')
