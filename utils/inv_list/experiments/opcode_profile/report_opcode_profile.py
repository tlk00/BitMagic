#!/usr/bin/env python3
"""Validate measured opcode accounting and show a compact collection report."""
from pathlib import Path
import json
import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "current"))
from validate_stream import read

root = Path(sys.argv[1])
limit = int(sys.argv[2])
run = read('profile', root, limit)
log = run.pop('log')
summary = re.findall(r'^Opcode profile blobs=(\d+) bytes=(\d+) header_bytes=(\d+) '
                     r'id_list_bytes=(\d+) records=(\d+) opcode_bytes=(\d+) '
                     r'payload_bytes=(\d+) reconciled=(\d+)$', log, re.M)
if len(summary) != 1:
    raise SystemExit('Expected one complete opcode profile summary')
blobs, total, header, idlist, records, opcode, payload, reconciled = map(int, summary[0])
if (blobs, total, opcode, reconciled) != (run['vectors'], run['size'], records, 1):
    raise SystemExit('Opcode profile size, count or reconciliation differs')
if header + idlist + opcode + payload != total:
    raise SystemExit('Opcode profile parts do not add to collection size')
if 'Opcode profile scope=top_level embedded_formula_bytes=in_parent_payload' not in log:
    raise SystemExit('Opcode scope is missing')
groups = []
section = re.split(r'Opcode group                 Records     (?:Opcode|Marker) bytes    Payload bytes      Total bytes   Percent\n', log, maxsplit=1)
if len(section) != 2:
    raise SystemExit('Missing opcode group table')
group_text = section[1].split('Opcode ID  Name', 1)[0]
for line in group_text.splitlines():
    match = re.match(r'^(.+?)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)$', line)
    if not match: continue
    name, count, opbytes, paybytes, bytes_, pct = match.groups()
    g = dict(name=name.strip(), count=int(count), opcode_bytes=int(opbytes),
             payload_bytes=int(paybytes), bytes=int(bytes_))
    if g['count'] != g['opcode_bytes'] or g['opcode_bytes']+g['payload_bytes'] != g['bytes']:
        raise SystemExit(f'Invalid opcode group row: {line}')
    groups.append(g)
if not groups or sum(g['bytes'] for g in groups)+header+idlist != total:
    raise SystemExit('Opcode group rows do not reconcile')
rows = []
for line in log.splitlines():
    match = re.match(r'^\s*(\d+)\s+(set_\w+|packed_zero_run|unlabelled_opcode)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)$', line)
    if not match: continue
    code, name, count, opbytes, paybytes, bytes_, pct = match.groups()
    values = tuple(map(int, (code, count, opbytes, paybytes, bytes_)))
    if values[2] != values[1] or values[2] + values[3] != values[4]:
        raise SystemExit(f'Invalid opcode row: {line}')
    rows.append(dict(code=int(code), name=name, count=int(count), opcode_bytes=int(opbytes),
                     payload_bytes=int(paybytes), bytes=int(bytes_)))
if sum(r['bytes'] for r in rows) + header + idlist != total or sum(r['count'] for r in rows) != records:
    raise SystemExit('Opcode rows do not add to collection totals')
# Original-ID cardinalities are collected separately from byte spans.
card = re.findall(r'^Opcode cardinality ids=(\d+) id_list_ids=(\d+) blobs=(\d+) reconciled=(\d+)$', log, re.M)
if len(card) != 1:
    raise SystemExit('Missing exact opcode cardinality summary; rebuild bminv')
ids, id_list_ids, card_blobs, card_ok = map(int, card[0])
if (ids, card_blobs, card_ok) != (run['ids'], blobs, 1):
    raise SystemExit('Opcode cardinalities differ from raw-ID verification')
eff = log.split('Opcode ID       Original IDs  Total bits/int  Payload bits/int\n', 1)
if len(eff) != 2: raise SystemExit('Missing opcode efficiency table')
counts = {}
for line in eff[1].splitlines():
    m = re.fullmatch(r'\s*(\d+)\s+(\d+)\s+(?:[\d.]+|n/a)\s+(?:[\d.]+|n/a)', line)
    if not m: break
    code, n = map(int, m.groups())
    if code in counts: raise SystemExit('Duplicate opcode cardinality')
    counts[code] = n
if set(counts) != {r['code'] for r in rows} or sum(counts.values()) + id_list_ids != ids:
    raise SystemExit('Per-opcode cardinalities do not reconcile')
for r in rows:
    r['ids'] = counts[r['code']]
    r['bits_per_int'] = 8*r['bytes']/r['ids'] if r['ids'] else None
    r['payload_bits_per_int'] = 8*r['payload_bytes']/r['ids'] if r['ids'] else None
run.update(profile=dict(blobs=blobs, bytes=total, header_bytes=header, id_list_bytes=idlist,
                        records=records, opcode_bytes=opcode, payload_bytes=payload, ids=ids, id_list_ids=id_list_ids, groups=groups, rows=rows))
(root/'profile'/'measurement.json').write_text(json.dumps(run, indent=2)+'\n')
print(f"Verified and decoded {run['vectors']:,} vectors, {run['ids']:,} original IDs")
print(f"Serialized bytes: {total:,}; weighted bits/ID: {8*total/run['ids']:.6f}")
print(f"Header bytes: {header:,}; ID-list bytes: {idlist:,}; opcode bytes: {opcode:,}; payload bytes: {payload:,}")
print(f"Encode/verify/decode wall seconds: {run['timings'][0]:.3f} / {run['timings'][1]:.3f} / {run['timings'][2]:.3f}")
print('Group                         Records   Marker bytes  Payload bytes    Total bytes   Share %')
for g in sorted(groups, key=lambda g: -g['bytes']):
    print(f"{g['name']:<24} {g['count']:>12,} {g['opcode_bytes']:>12,} "
          f"{g['payload_bytes']:>14,} {g['bytes']:>14,} {100*g['bytes']/total:>9.3f}")
print('Opcode ID  Name                                  Records   Marker bytes  Payload bytes    Total bytes   Share %    Original IDs  Bits/int  Payload b/int')
for r in rows:
    bpi = f"{r['bits_per_int']:.6f}" if r['ids'] else 'n/a'
    payload_bpi = f"{r['payload_bits_per_int']:.6f}" if r['ids'] else 'n/a'
    print(f"{r['code']:>9}  {r['name']:<34} {r['count']:>12,} {r['opcode_bytes']:>12,} "
          f"{r['payload_bytes']:>14,} {r['bytes']:>14,} {100*r['bytes']/total:>9.3f} "
          f"{r['ids']:>15,} {bpi:>10} {payload_bpi:>14}")
print('Marker bytes = one top-level opcode byte per record; Opcode ID is the numeric format identifier.')
print('Payload bytes include metadata, nested opcodes, formulas and padding; total = marker + payload bytes.')
print('Scope: top-level records; embedded formula bytes are included in their parent payload.')
print('Accounting and raw-ID verification: PASS')

print('Bits/int = 8 * complete record bytes / original IDs represented; weighted, not a mean of record ratios.')
print('Payload b/int excludes only the marker byte. Zero-ID records show n/a; collection bits/ID includes headers.')
