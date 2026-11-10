#!/usr/bin/env python3
"""Compare verified complete streams and show prefix/full opcode efficiency."""
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
full = {
    'unsorted': Path('/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-unsorted-opcode-bpi-20261102-091640/unsorted/profile/measurement.json'),
    'sorted': Path('/Volumes/USB_1TB/BitMagic-experiments/2026-11-01/runs/inv_list/gov2-highest-clean-20261101-160939/sorted/profile/measurement.json'),
}
def load(path):
    return json.loads(path.read_text())
def metric(row, key, digits=4):
    value = row.get(key)
    return 'n/a' if value is None else f'{value:.{digits}f}'
manifest_path = root.parent/'sample'/'sample_manifest.json'
manifest = load(manifest_path) if manifest_path.exists() else None
playlist_path = root.parent/'inputs'/'playlist_manifest.json'
playlist = load(playlist_path) if playlist_path.exists() else None
label = 'Playlist' if playlist else ('Sample' if manifest else 'Prefix')
if playlist:
    print('Workload: independent opcode-ranked playlists; see inputs/*.mapping.csv for original vector indices.')
if manifest:
    print('Sampling: '+manifest['method']+'; identical indices for sorted and unsorted.')
    print(f"Index range: {manifest['indices_zero_based'][0]:,}..{manifest['indices_zero_based'][-1]:,}; full vector count: {manifest['total_vectors']:,}.")
results = {}
for dataset in ('sorted', 'unsorted'):
    modes = {m: load(root/m/dataset/'profile/measurement.json') for m in ('baseline','selected')}
    a,b = modes['baseline'],modes['selected']
    if (a['vectors'],a['ids']) != (b['vectors'],b['ids']):
        raise SystemExit(f'{dataset}: comparison scope differs')
    for mode,d in modes.items():
        p=root/mode/dataset/'profile'
        if 'Accounting and raw-ID verification: PASS' not in (p/'report.txt').read_text():
            raise SystemExit(f'{dataset}/{mode}: verification report missing')
        if sum(r['ids'] for r in d['profile']['rows'])+d['profile']['id_list_ids'] != d['ids']:
            raise SystemExit('ID accounting differs')
    if playlist:
        pm=playlist['datasets'][dataset]
        if (a['vectors'],a['ids']) != (pm['selected_count'],pm['selected_vector_ids']):
            raise SystemExit('Encoded scope differs from playlist metadata')
        # Baseline must reproduce the exact opcode population that defined ranking.
        pr=next((r for r in a['profile']['rows'] if r['code']==pm['opcode']),{})
        if (pr.get('bytes'),pr.get('ids'),pr.get('count')) != (pm['selected_target_bytes'],pm['selected_target_ids'],pm['selected_target_records']):
            raise SystemExit('Replay baseline target costs differ from frozen survey')
    if manifest and (a['vectors'],a['ids']) != (manifest['sampled_vectors'],manifest['datasets'][dataset]['ids']):
        raise SystemExit('Encoded scope differs from sample manifest')
    saved=a['size']-b['size']
    if saved<0: raise SystemExit(f'{dataset}: selected stream grew by {-saved} bytes')
    print(f'\nGov2.{dataset}: {a["vectors"]:,} vectors, {a["ids"]:,} IDs')
    print('Mode          Bytes          Bits/ID    Encode s   Verify s   Decode s')
    for mode,d in modes.items():
        print(f'{mode:<10} {d["size"]:>14,} {8*d["size"]/d["ids"]:>12.6f} '+
              ' '.join(f'{t:>10.3f}' for t in d['timings']))
    print(f'Final savings: {saved:,} bytes ({100*saved/a["size"]:.6f}%).')
    ref=load(full[dataset]) if full[dataset].exists() else None
    reference_name=str(full[dataset])
    if playlist:
        survey=playlist['datasets'][dataset].get('survey_summary')
        if survey:
            ref=dict(vectors=survey['vectors'],ids=survey['ids'],size=survey['bytes'],
                profile=dict(rows=[dict(code=int(c),count=r['records'],bytes=r['bytes'],ids=r['ids'],bits_per_int=r['bits_per_id']) for c,r in survey['opcodes'].items()]))
            reference_name=playlist['datasets'][dataset]['survey_directory']+'/'+dataset+'/survey_summary.json'
    refrows={r['code']:r for r in ref['profile']['rows']} if ref else {}
    ar={r['code']:r for r in a['profile']['rows']};br={r['code']:r for r in b['profile']['rows']}
    print(f'Opcode   Base bytes   Selected bytes   Base bits/ID   Selected bits/ID   {label} share%   Full share%   Full bits/ID')
    for code in sorted(ar.keys()|br.keys(),key=lambda c:-ar.get(c,{}).get('bytes',0)):
        x=ar.get(code,{}); y=br.get(code,{}); f=refrows.get(code,{})
        full_share=f'{100*f.get("bytes",0)/ref["size"]:.4f}' if ref else 'n/a'
        print(f'{code:>6} {x.get("bytes",0):>12,} {y.get("bytes",0):>16,} '
              f'{metric(x,"bits_per_int"):>14} {metric(y,"bits_per_int"):>18} '
              f'{100*x.get("bytes",0)/a["size"]:>15.4f} {full_share:>13} {metric(f,"bits_per_int"):>14}')
    if ref:
        print(f'{label} contains {100*a["vectors"]/ref["vectors"]:.3f}% of vectors but {100*a["ids"]/ref["ids"]:.3f}% of full IDs.')
        print(f'Collection bits/ID: {label.lower()} {8*a["size"]/a["ids"]:.6f}, full {8*ref["size"]/ref["ids"]:.6f}.')
    if ref:
        print('Opcode-68 coverage:')
        x=ar.get(68,{}) ; f=refrows.get(68,{})
        if f.get('count'):
            print(f"  records {x.get('count',0):,}/{f['count']:,} ({100*x.get('count',0)/f['count']:.4f}%); bytes {100*x.get('bytes',0)/f['bytes']:.4f}% of full opcode-68 bytes")
    if not ref or any('bits_per_int' not in r for r in refrows.values()):
        print('Full per-opcode bits/ID unavailable for this reference; n/a is not evidence of no skew.')
    results[dataset]={'baseline':a,'selected':b,'saved_bytes':saved,'full_reference':reference_name,'full':ref,'sampling':manifest,'playlist':playlist['datasets'][dataset] if playlist else None}
print(f'\n{label} results are workload-specific. Opcode shares and bits/ID can reveal skew, not prove representativeness.')
print('Both arms include all retained selectors and opcode profiling; encode timings include profiling overhead.')
print('Complete stream verification, independent decode, byte and ID accounting: PASS')
(root/'comparison.json').write_text(json.dumps(results,indent=2)+'\n')
