#!/usr/bin/env python3
"""Presentation/byte-identity checks using verified streams and opcode profiles."""
import argparse,json
from pathlib import Path
parser=argparse.ArgumentParser();parser.add_argument('root',type=Path);parser.add_argument('--reference',type=Path)
a=parser.parse_args();results=[]
def require(ok,msg):
    if not ok: raise SystemExit(msg)
for dataset in ('sorted','unsorted'):
    p=a.root/dataset/'profile'
    if not p.exists():continue
    d=json.loads((p/'measurement.json').read_text())
    require('Verification complete.' in (p/'verify.log').read_text(),f'{dataset}: missing raw verification')
    require('Decode complete.' in (p/'decode.log').read_text(),f'{dataset}: missing decode')
    require('Accounting and raw-ID verification: PASS' in (p/'report.txt').read_text(),f'{dataset}: missing accounting')
    if a.reference:
        old=a.reference/dataset/'profile';ref=json.loads((old/'measurement.json').read_text())
        require(all(ref[k]==d[k] for k in ('vectors','ids','size')),f'{dataset}: reference scope/size differs')
        require((old/'stream.sha256').read_text().split()[0]==(p/'stream.sha256').read_text().split()[0],f'{dataset}: reference SHA differs')
    results.append((dataset,d))
require(results,'No verified datasets')
print('Compression preset: see recorded -level and overrides in run.txt')
print('Dataset        Vectors               IDs             Bytes    Bits/ID   Encode s   Verify s   Decode s')
for ds,d in results:
    print(f"{ds:<10} {d['vectors']:>11,} {d['ids']:>17,} {d['size']:>17,} {8*d['size']/d['ids']:>10.6f} "+' '.join(f'{t:>10.3f}' for t in d['timings']))
print('Raw-ID verification, independent decode and opcode/stream accounting: PASS')
if a.reference:print('SHA-256 byte identity with pre-cleanup reference: PASS')
print('Configuration and BWT cutoff: recorded per-phase commands in run.txt.')
if any(any(r['code']==98 for r in d['profile']['rows']) for _,d in results):
    print('Optional high-count/polarity superblock opcode 98 selected.')
(a.root/'verification-measurement.json').write_text(json.dumps({ds:d for ds,d in results},indent=2)+'\n')
