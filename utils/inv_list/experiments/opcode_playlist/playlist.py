#!/usr/bin/env python3
"""Rank a verified per-vector opcode survey; extract/replay serialized playlists."""
import argparse,csv,functools,hashlib,itertools,json,re,shutil,struct,subprocess
from pathlib import Path
FIELDS=['vector_index','source_offset','vector_ids','vector_bytes','records','bytes','ids']

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def rank_compare(a,b):
    # Descending exact ratio, then bytes; ascending index. No rounded-score ties.
    lhs=a['bytes']*b['ids'];rhs=b['bytes']*a['ids']
    if lhs!=rhs:return -1 if lhs>rhs else 1
    if a['bytes']!=b['bytes']:return -1 if a['bytes']>b['bytes'] else 1
    return (a['vector_index']>b['vector_index'])-(a['vector_index']<b['vector_index'])

def byte_rank_key(row):
    return (-row['bytes'],row['vector_index'])

def select(args):
    root=args.survey
    config=json.loads((root/'configuration.json').read_text())
    all_meta={}
    for ds in ('sorted','unsorted'):
        directory=root/ds
        if not (directory/'COMPLETE.json').exists():raise ValueError(f'{ds}: survey is not complete')
        complete=json.loads((directory/'COMPLETE.json').read_text())
        source=Path(complete['source']);stamp=source.stat()
        if (stamp.st_size,stamp.st_mtime_ns)!=(complete['source_bytes'],complete['source_mtime_ns']):
            raise ValueError('Survey source identity changed')
        candidates=[];vectors=ids=stream_bytes=records=0;target_bytes=target_ids=target_records=0
        opcode_totals={}
        with (directory/'opcodes.csv').open() as f:
            reader=csv.DictReader(f)
            expected_offset=0
            for index,group in itertools.groupby(reader,lambda r:int(r['vector_index'])):
                rows=[{k:int(v) for k,v in r.items()} for r in group]
                first=rows[0]
                if index!=vectors or first['source_offset']!=expected_offset:raise ValueError('Nonsequential survey rows')
                if len({r['opcode'] for r in rows})!=len(rows):raise ValueError('Duplicate per-vector opcode')
                for r in rows:
                    if any(r[k]!=first[k] for k in ('vector_index','source_offset','vector_ids','vector_bytes')):raise ValueError('Survey row scope differs')
                    if not -2<=r['opcode']<256:raise ValueError('Invalid opcode value')
                    if r['bytes']<r['records'] or min(r['bytes'],r['ids'],r['records'])<0:raise ValueError('Invalid opcode costs')
                    if r['opcode']>=0:
                        t=opcode_totals.setdefault(r['opcode'],[0,0,0]);t[0]+=r['records'];t[1]+=r['bytes'];t[2]+=r['ids']
                        records+=r['records']
                if sum(r['bytes'] for r in rows)!=first['vector_bytes'] or sum(r['ids'] for r in rows)!=first['vector_ids']:
                    raise ValueError('Per-vector bytes/IDs fail reconciliation')
                vectors+=1;ids+=first['vector_ids'];stream_bytes+=first['vector_bytes']
                expected_offset+=4+4*first['vector_ids']
                target=next((r for r in rows if r['opcode']==args.opcode),None)
                if target:
                    target_bytes+=target['bytes'];target_ids+=target['ids'];target_records+=target['records']
                    if target['bytes']>=args.min_bytes and target['ids']:
                        candidates.append({k:target[k] for k in FIELDS})
        log=(directory/'survey.log').read_text()
        def number(pattern):
            found=re.findall(pattern,log,re.M)
            if len(found)!=1:raise ValueError(f'Missing/duplicate survey result: {pattern}')
            return int(found[0])
        if (vectors,ids,stream_bytes)!=(number(r'^Total vectors=(\d+)$'),number(r'^Total ints=(\d+)$'),number(r'^BV size = (\d+)$')):
            raise ValueError('Survey aggregate scope differs from compression log')
        if (complete.get('source_bytes'),complete.get('source_mtime_ns'))!=(stamp.st_size,stamp.st_mtime_ns):
            raise ValueError('Completed survey identity differs')
        if vectors!=number(r'^Opcode survey verified_vectors=(\d+) rows_file=.*$') or expected_offset!=stamp.st_size:
            raise ValueError('Survey did not verify and consume full input')
        summary=re.findall(r'^Opcode profile blobs=(\d+) bytes=(\d+) header_bytes=(\d+) id_list_bytes=(\d+) records=(\d+) opcode_bytes=(\d+) payload_bytes=(\d+) reconciled=(\d+)$',log,re.M)
        if len(summary)!=1 or tuple(map(int,summary[0]))[0:2]!=(vectors,stream_bytes) or int(summary[0][-1])!=1 or int(summary[0][4])!=records:
            raise ValueError('Aggregate opcode accounting differs')
        # Independently reconcile CSV opcode totals against the printed aggregate.
        packed=[0,0,0]
        normalized={}
        for code,counts in opcode_totals.items():
            if code>=128:
                packed=[x+y for x,y in zip(packed,counts)]
            else:normalized[code]=counts
        if packed[0]:normalized[128]=packed
        printed={}
        for m in re.finditer(r'^\s*(\d+)\s+(?:set_\w+|packed_zero_run|unlabelled_opcode)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)$',log,re.M):
            code,count,marker,payload,total,_=m.groups()
            if int(marker)+int(payload)!=int(total):raise ValueError('Printed opcode row does not reconcile')
            printed[int(code)]=[int(count),int(total)]
        if printed!={code:counts[:2] for code,counts in normalized.items()}:
            raise ValueError('CSV opcode totals differ from profiler')
        cardinalities=re.findall(r'^Opcode cardinality ids=(\d+) id_list_ids=(\d+) blobs=(\d+) reconciled=(\d+)$',log,re.M)
        if len(cardinalities)!=1 or int(cardinalities[0][0])!=ids or int(cardinalities[0][2])!=vectors or int(cardinalities[0][3])!=1:
            raise ValueError('Survey ID reconciliation is missing')
        summary_data=dict(vectors=vectors,ids=ids,bytes=stream_bytes,opcodes={str(c):dict(records=v[0],bytes=v[1],ids=v[2],bits_per_id=8*v[1]/v[2] if v[2] else None) for c,v in opcode_totals.items()})
        (directory/'survey_summary.json').write_text(json.dumps(summary_data,indent=2)+'\n')
        ranking=getattr(args,'ranking','high-bpi')
        candidates.sort(key=byte_rank_key if ranking=='high-bytes' else functools.cmp_to_key(rank_compare));chosen=candidates[:args.top]
        if not chosen:raise ValueError(f'{ds}: no eligible opcode-{args.opcode} vectors')
        name=f'gov2.{ds}.op{args.opcode}.{ranking}.min{args.min_bytes}B.top{args.top}'
        dest=root/'playlists';dest.mkdir(exist_ok=True)
        bv=dest/(name+'.bv');meta_path=dest/(name+'.json');ranked=dest/(name+'.csv')
        if any(p.exists() for p in (bv,meta_path,ranked)):raise ValueError('Playlist already exists')
        txt=dest/(name+'.indices.txt')
        txt.write_text(''.join(str(r['vector_index'])+'\n' for r in chosen))
        subprocess.run([str(args.tool.resolve()),'pack',str(txt),str(bv)],check=True)
        txt.unlink()
        # Bind selected records to their exact original bytes without hashing
        # both 24-GB sources again on every replay.
        with source.open('rb') as raw:
            for r in chosen:
                raw.seek(r['source_offset']);left=4+4*r['vector_ids'];digest=hashlib.sha256()
                while left:
                    chunk=raw.read(min(left,1024*1024))
                    if not chunk:raise ValueError('Selected source record truncated')
                    digest.update(chunk);left-=len(chunk)
                r['source_record_sha256']=digest.hexdigest()
        after=source.stat()
        if (after.st_size,after.st_mtime_ns)!=(stamp.st_size,stamp.st_mtime_ns):raise ValueError('Source changed during selection')
        with ranked.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['rank']+FIELDS+['source_record_sha256','bits_per_id']);w.writeheader()
            for i,r in enumerate(chosen):w.writerow(dict(rank=i+1,**r,bits_per_id=8*r['bytes']/r['ids']))
        sb=sum(r['bytes'] for r in chosen);si=sum(r['ids'] for r in chosen)
        meta=dict(format='BM_OPCODE_PLAYLIST_V1',dataset=ds,opcode=args.opcode,index_base=0,
            survey_directory=str(root.resolve()),survey_summary=summary_data,
            source=complete,configuration=config,ranking=('descending target bytes; ascending original index' if ranking=='high-bytes' else 'descending exact target bytes/IDs; descending target bytes; ascending original index'),
            ranking_mode=ranking,
            min_target_bytes=args.min_bytes,requested_count=args.top,selected_count=len(chosen),eligible_vectors=len(candidates),
            total_vectors=vectors,total_ids=ids,total_stream_bytes=stream_bytes,target_bytes=target_bytes,target_ids=target_ids,target_records=target_records,
            selected_target_bytes=sb,selected_target_ids=si,selected_vector_ids=sum(r['vector_ids'] for r in chosen),
            selected_target_records=sum(r['records'] for r in chosen),playlist_sha256=sha(bv),ranked_csv=ranked.name,ranked_sha256=sha(ranked))
        meta_path.write_text(json.dumps(meta,indent=2)+'\n');all_meta[ds]=meta
        print(f'Gov2.{ds}: {vectors:,} vectors; {ids:,} IDs; full stream equivalent {stream_bytes:,} bytes; {8*stream_bytes/ids:.6f} bits/ID')
        print('Opcode     Records          Bytes    Share %          IDs    Bits/ID')
        for code,vals in sorted(opcode_totals.items(),key=lambda kv:-kv[1][1]):
            n,bs,ns=vals;bpi=f'{8*bs/ns:.6f}' if ns else 'n/a'
            print(f'{code:>6} {n:>11,} {bs:>15,} {100*bs/stream_bytes:>10.4f} {ns:>14,} {bpi:>11}')
        print(f'{ds}: eligible={len(candidates):,}; selected={len(chosen):,}; target bytes={sb:,}/{target_bytes:,} ({100*sb/target_bytes:.3f}%); target IDs={si:,}; weighted target bits/ID={8*sb/si:.6f}')
        print(f'  {bv}')
    summary_name='selection-summary.json' if ranking=='high-bpi' else 'selection-summary-high-bytes.json'
    (root/summary_name).write_text(json.dumps(all_meta,indent=2)+'\n')

def extract(args):
    meta=json.loads(args.playlist.with_suffix('.json').read_text())
    if meta['format']!='BM_OPCODE_PLAYLIST_V1' or sha(args.playlist)!=meta['playlist_sha256']:
        raise ValueError('Playlist metadata/hash mismatch')
    source=Path(meta['source']['source']) if args.source is None else args.source
    stamp=source.stat()
    if (stamp.st_size,stamp.st_mtime_ns)!=(meta['source']['source_bytes'],meta['source']['source_mtime_ns']):
        raise ValueError('Source size/mtime differs from surveyed collection')
    ranked=args.playlist.parent/meta['ranked_csv']
    if sha(ranked)!=meta['ranked_sha256']:raise ValueError('Ranked measurements changed')
    with ranked.open() as f:rows=sorted((dict(**{k:int(r[k]) for k in FIELDS},source_record_sha256=r['source_record_sha256']) for r in csv.DictReader(f)),key=lambda r:r['vector_index'])
    args.out.parent.mkdir(parents=True,exist_ok=True)
    decoded=args.out.with_name(args.out.name+'.indices.txt')
    subprocess.run([str(args.tool.resolve()),'unpack',str(args.playlist),str(decoded)],check=True)
    indices=[int(s) for s in decoded.read_text().splitlines()];decoded.unlink()
    if indices!=[r['vector_index'] for r in rows] or len(rows)!=meta['selected_count']:raise ValueError('Playlist membership differs from measurements')
    if shutil.disk_usage(args.out.parent).free<sum(4+4*r['vector_ids'] for r in rows)+64*1024*1024:
        raise OSError('Insufficient free space for replay inputs')
    with source.open('rb') as f,args.out.open('xb') as out:
        for r in rows:
            f.seek(r['source_offset']);header=f.read(4)
            if len(header)!=4 or struct.unpack('<I',header)[0]!=r['vector_ids']:raise ValueError('Source record framing differs')
            out.write(header);left=4*r['vector_ids'];digest=hashlib.sha256(header)
            while left:
                chunk=f.read(min(left,1024*1024))
                if not chunk:raise ValueError('Truncated replay input')
                out.write(chunk);digest.update(chunk);left-=len(chunk)
            if digest.hexdigest()!=r['source_record_sha256']:raise ValueError('Source record bytes differ from selected survey input')
    after=source.stat()
    if (stamp.st_size,stamp.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise ValueError('Source changed while extracting')
    mapping=args.out.with_name(args.out.name+'.mapping.csv')
    with mapping.open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['replay_index','original_vector_index','source_offset','ids'])
        for i,r in enumerate(rows):w.writerow([i,r['vector_index'],r['source_offset'],r['vector_ids']])
    print(f'Extracted {meta["dataset"]}: {len(rows):,} vectors; {args.out.stat().st_size:,} bytes')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    s=sub.add_parser('select');s.add_argument('--survey',type=Path,required=True);s.add_argument('--opcode',type=int,default=68);s.add_argument('--top',type=int,default=10000);s.add_argument('--min-bytes',type=int,default=256);s.add_argument('--tool',type=Path,default=Path('./gov2_playlist'))
    s.add_argument('--ranking',choices=('high-bpi','high-bytes'),default='high-bpi')
    e=sub.add_parser('extract');e.add_argument('--playlist',type=Path,required=True);e.add_argument('--source',type=Path);e.add_argument('--out',type=Path,required=True);e.add_argument('--tool',type=Path,default=Path('./gov2_playlist'))
    a=p.parse_args()
    if a.action=='select':
        if not 0<=a.opcode<256 or a.top<1 or a.min_bytes<1:p.error('Invalid opcode/top/min-bytes')
        select(a)
    else:extract(a)
