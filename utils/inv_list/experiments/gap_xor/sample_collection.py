#!/usr/bin/env python3
"""Deterministic evenly spaced complete-vector sampling, preserving raw IDs."""
from pathlib import Path
import argparse, json, struct
p=argparse.ArgumentParser()
p.add_argument('input',type=Path);p.add_argument('output',type=Path)
p.add_argument('--count',type=int,default=10000)
a=p.parse_args()
if a.count<1: p.error('count must be positive')
if a.output.exists(): p.error('output already exists')
before=a.input.stat()
size=before.st_size
# Scan record headers only; do not load the collection or all offsets into RAM.
def records(f):
    offset=0; index=0
    while offset<size:
        f.seek(offset); header=f.read(4)
        if len(header)!=4: raise SystemExit('Truncated vector length')
        n=struct.unpack('<I',header)[0]
        end=offset+4+4*n
        if not n or end>size: raise SystemExit(f'Invalid vector at offset {offset}')
        yield index,offset,end,n
        offset=end;index+=1
with a.input.open('rb') as f: total=sum(1 for _ in records(f))
if not total: raise SystemExit('Empty collection')
count=min(a.count,total)
selected={i*(total-1)//(count-1) for i in range(count)} if count>1 else {total//2}
manifest=[];ids=0
with a.input.open('rb') as f,a.output.open('xb') as out:
    for index,offset,end,n in records(f):
        if index not in selected:continue
        f.seek(offset);left=end-offset
        while left:
            chunk=f.read(min(left,1024*1024))
            if not chunk:raise SystemExit('Input changed or truncated during copy')
            out.write(chunk);left-=len(chunk)
        manifest.append(dict(index=index,offset=offset,ids=n));ids+=n
after=a.input.stat()
if (after.st_size,after.st_mtime_ns)!=(before.st_size,before.st_mtime_ns):
    raise SystemExit('Input changed during sampling')
meta=dict(source=str(a.input.resolve()),source_bytes=size,total_vectors=total,
          sampled_vectors=count,original_ids=ids,method='evenly spaced vector indices including endpoints',records=manifest)
a.output.with_suffix('.manifest.json').write_text(json.dumps(meta,indent=2)+'\n')
print(f'Sample: {count:,} of {total:,} vectors; {ids:,} IDs; {a.output.stat().st_size:,} bytes')
