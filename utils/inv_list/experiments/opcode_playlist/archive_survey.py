#!/usr/bin/env python3
"""Move a completed USB survey into the persistent survey collection."""
import argparse,json,os,subprocess
from datetime import datetime,timezone
from pathlib import Path
from playlist import sha


def archive(source,library):
    source=source.resolve();library=library.resolve()
    if source==library or source in library.parents:
        raise ValueError('Library must be outside the source survey')
    if (source/'FAILURE.txt').exists():raise ValueError('Survey has a failure marker')
    for name in ('configuration.json','selection-summary.json','report.txt','sorted/COMPLETE.json','unsorted/COMPLETE.json'):
        if not (source/name).is_file():raise ValueError(f'Survey is not complete: missing {name}')
    summaries=json.loads((source/'selection-summary.json').read_text())
    if set(summaries)!= {'sorted','unsorted'}:raise ValueError('Incomplete paired selection summary')
    playlists=source/'playlists'
    bv_paths=sorted(p for p in playlists.glob('*.bv') if not p.name.startswith('.'))
    if len(bv_paths)<2:raise ValueError('Missing paired playlists')
    for bv in bv_paths:
        meta=json.loads(bv.with_suffix('.json').read_text())
        if meta['format']!='BM_OPCODE_PLAYLIST_V1' or sha(bv)!=meta['playlist_sha256']:
            raise ValueError(f'Invalid playlist: {bv}')
        if sha(playlists/meta['ranked_csv'])!=meta['ranked_sha256']:
            raise ValueError(f'Changed ranked CSV: {bv}')
    library.mkdir(parents=True,exist_ok=True)
    destination=library/source.name
    if destination.exists():raise ValueError('Destination already exists; no overwrite')
    if source.stat().st_dev!=library.stat().st_dev:
        raise ValueError('Archive requires the same filesystem (rename only)')
    os.rename(source,destination)
    (destination/'relocation.json').write_text(json.dumps(dict(
        original_directory=str(source),survey_directory=str(destination),
        moved_utc=datetime.now(timezone.utc).isoformat(),
        method='same-filesystem rename; original artifacts preserved'),indent=2)+'\n')
    print(f'Archived survey: {destination}')
    return destination

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--survey',type=Path,required=True)
    p.add_argument('--library',type=Path,default=Path('/Volumes/USB_1TB/BitMagic-experiments/surveys'))
    a=p.parse_args()
    if not str(a.survey.resolve()).startswith('/Volumes/USB_1TB/') or not str(a.library.resolve()).startswith('/Volumes/USB_1TB/'):
        p.error('Both paths must be on USB_1TB')
    if ' on /Volumes/USB_1TB (' not in subprocess.check_output(['mount'],text=True):p.error('USB_1TB is not mounted')
    archive(a.survey,a.library)
