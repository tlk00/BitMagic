#!/usr/bin/env python3
"""Evaluate both opcode playlists using the existing verified correction comparison."""
import argparse,json,subprocess
from pathlib import Path
from common import usb_output
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--survey',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
p.add_argument('--bin',type=Path,default=Path('./bminv'));p.add_argument('--tool',type=Path,default=Path('./gov2_playlist'))
p.add_argument('--threads',type=int,default=10);p.add_argument('--progress',type=int,default=1000)
p.add_argument('--opcode',type=int,default=68);p.add_argument('--top',type=int,default=10000);p.add_argument('--min-bytes',type=int,default=256)
p.add_argument('--input-dir',type=Path,help='Optional same-identity raw input location')
p.add_argument('--ranking',choices=('high-bpi','high-bytes'),default='high-bpi')
a=p.parse_args()
if not 1<=a.threads<=256 or min(a.top,a.min_bytes,a.progress)<1:p.error('Invalid numeric option')
usb_output(a.out)
try:
    inputs=a.out/'inputs';inputs.mkdir();manifest=dict(method=f'opcode {a.ranking} playlists',datasets={})
    for ds in ('sorted','unsorted'):
        playlist=a.survey/'playlists'/f'gov2.{ds}.op{a.opcode}.{a.ranking}.min{a.min_bytes}B.top{a.top}.bv'
        meta=json.loads(playlist.with_suffix('.json').read_text())
        # Artifact membership/hashes stay frozen; reports use the actual survey
        # location after a completed survey has been archived.
        meta['survey_directory_created']=meta.get('survey_directory')
        meta['survey_directory']=str(a.survey.resolve())
        manifest['datasets'][ds]=meta
        cmd=['python3',str(Path(__file__).with_name('playlist.py')),'extract','--playlist',str(playlist),
             '--tool',str(a.tool.resolve()),'--out',str(inputs/('gov2.'+ds))]
        if a.input_dir:cmd+=['--source',str(a.input_dir/('gov2.'+ds))]
        subprocess.run(cmd,check=True)
    (inputs/'playlist_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    # No vector limit: independent lists may contain different numbers of entries.
    harness=Path(__file__).resolve().parent.parent/'sblock_correction/run_gov2_correction.sh'
    cmd=['bash',str(harness),'--bin',str(a.bin.resolve()),'--input-dir',str(inputs),
         '--threads',str(a.threads),'--progress',str(a.progress),'--max-vectors','0','--out',str(a.out/'comparison')]
    subprocess.run(cmd,check=True)
    report=(a.out/'comparison/report.txt').read_text();(a.out/'report.txt').write_text(report)
    print(f'Playlist replay reports: {a.out}')
except BaseException as e:
    (a.out/'FAILURE.txt').write_text(str(e)+'\n');raise
