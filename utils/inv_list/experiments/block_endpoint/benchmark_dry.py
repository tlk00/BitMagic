#!/usr/bin/env python3
"""ABBA comparison of separately built, single-thread bench_endpoint binaries."""
import argparse,hashlib,json,statistics,subprocess,sys,time
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument('--before',required=True);p.add_argument('--after',required=True)
p.add_argument('--input',default='/Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.unsorted')
p.add_argument('--out',required=True);p.add_argument('--max-vectors',type=int,default=10000)
a=p.parse_args();root=Path(a.out).resolve()
if a.max_vectors <= 0:sys.exit('--max-vectors must be positive')
if not str(root).startswith('/Volumes/USB_1TB/') or not Path('/Volumes/USB_1TB').is_mount():sys.exit('Require mounted USB_1TB output')
root.mkdir(parents=True,exist_ok=False);rows=[]
for i,variant in enumerate(('before','after','after','before')):
 name=f'{i+1}-{variant}';d=root/name;d.mkdir()
 cmd=[str(Path(getattr(a,variant)).resolve()),a.input,str(d/'candidate.bv'),str(a.max_vectors)]
 (d/'command.json').write_text(json.dumps(cmd));print(name,flush=True);start=time.perf_counter()
 with (d/'progress.log').open('w') as log:
  result=subprocess.run(cmd,stdout=subprocess.PIPE,stderr=log,text=True,check=True)
 row=json.loads(result.stdout);row.update(name=name,variant=variant,wall_s=time.perf_counter()-start)
 if row['vectors']!=a.max_vectors or not row['verified']:sys.exit('Incomplete validation')
 digest=hashlib.sha256()
 with (d/'candidate.bv').open('rb') as stream:
  for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
 row['sha256']=digest.hexdigest();rows.append(row)
 (root/'results.json').write_text(json.dumps(rows,indent=2));print(row,flush=True)
if len({r['sha256'] for r in rows})!=1:sys.exit('FAIL: streams differ')
before=statistics.mean(r['encode_s'] for r in rows if r['variant']=='before');after=statistics.mean(r['encode_s'] for r in rows if r['variant']=='after')
report='Variant | Encode s | Decode s | Wall s | Bytes\n'
for r in rows:report+=f"{r['name']} | {r['encode_s']:.3f} | {r['decode_s']:.3f} | {r['wall_s']:.3f} | {r['bytes']}\n"
report+=f'Encode means: before={before:.3f}s after={after:.3f}s reduction={(1-after/before)*100:.2f}% speedup={before/after:.3f}x\nAll four streams round-trip verified and SHA256 identical.\n'
(root/'report.txt').write_text(report);print(report)
