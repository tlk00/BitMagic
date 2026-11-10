#!/usr/bin/env python3
"""Full retained-codec Gov2 survey, then reusable top-N opcode playlists."""
import argparse,json,subprocess
from pathlib import Path
from common import usb_output,logged
from playlist import sha
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--bin',type=Path,default=Path('./bminv'));p.add_argument('--tool',type=Path,default=Path('./gov2_playlist'))
p.add_argument('--input-dir',type=Path,default=Path('/Volumes/USB_1TB/int_set_bench/Gov2Flat'))
p.add_argument('--out',type=Path,required=True);p.add_argument('--threads',type=int,default=10);p.add_argument('--progress',type=int,default=1000)
p.add_argument('--top',type=int,default=10000);p.add_argument('--opcode',type=int,default=68);p.add_argument('--min-bytes',type=int,default=256)
a=p.parse_args()
if not 1<=a.threads<=256 or min(a.progress,a.top,a.min_bytes)<1 or not 0<=a.opcode<256:p.error('Invalid numeric option')
for binary in (a.bin,a.tool):
    if not binary.is_file():p.error(f'Missing built binary: {binary}')
for ds in ('sorted','unsorted'):
    if not (a.input_dir/('gov2.'+ds)).is_file():p.error(f'Missing Gov2.{ds}')
usb_output(a.out)
try:
    config=dict(format='BM_OPCODE_SURVEY_V1',binary_sha256=sha(a.bin),
        level=6,threads=a.threads,codecs='step2/interleaved4,BWT shared min16,GAP-XOR,tolerance100,adaptive92,wave90/97,endpoint-select65',
        correction_selection=False,verification='every compressed BLOB round-tripped against its optimized original vector',
        index_base=0,arguments=vars(a).copy())
    config['arguments']={k:str(v) if isinstance(v,Path) else v for k,v in config['arguments'].items()}
    (a.out/'configuration.json').write_text(json.dumps(config,indent=2)+'\n')
    for ds in ('sorted','unsorted'):
        source=(a.input_dir/('gov2.'+ds)).resolve();stamp=source.stat()
        d=a.out/ds;d.mkdir()
        cmd=[str(a.bin.resolve()),'-t','-s','-level','6','-u32in',str(source),'-bvout','/dev/null',
            '-block-select-step2','-periodic-step1-interleaved','4','-bwt-min-gain','16','-gap-xor',
            '-sblock-gap-tolerance','100','-sblock-adaptive','-sblock-wave-adaptive',
            '-block-endpoint-select','-bwt-shared','-threads',str(a.threads),'-batch-vectors','512',
            '-progress-vectors',str(a.progress),'-opcode-survey',str(d/'opcodes.csv')]
        print(f'Survey Gov2.{ds}: {source}',flush=True)
        (d/'command.json').write_text(json.dumps(cmd,indent=2)+'\n')
        seconds=logged(cmd,d/'survey.log')
        after=source.stat()
        if (stamp.st_size,stamp.st_mtime_ns)!=(after.st_size,after.st_mtime_ns):raise ValueError('Source changed during survey')
        (d/'COMPLETE.json').write_text(json.dumps(dict(source=str(source),source_bytes=stamp.st_size,
            source_mtime_ns=stamp.st_mtime_ns,seconds=seconds),indent=2)+'\n')
    cmd=['python3',str(Path(__file__).with_name('playlist.py')),'select','--survey',str(a.out),
        '--tool',str(a.tool.resolve()),'--top',str(a.top),'--opcode',str(a.opcode),'--min-bytes',str(a.min_bytes)]
    report=subprocess.check_output(cmd,text=True)
    (a.out/'report.txt').write_text(report);print(report,end='')
    print(f'Survey and playlists: {a.out}')
except BaseException as e:
    (a.out/'FAILURE.txt').write_text(str(e)+'\n');raise
