"""Shared output validation and logged progress for survey/replay harnesses."""
import os,shutil,subprocess,time
from pathlib import Path


def usb_output(path):
    if not str(path.resolve()).startswith('/Volumes/USB_1TB/'):
        raise ValueError('Output must be under /Volumes/USB_1TB')
    if ' on /Volumes/USB_1TB (' not in subprocess.check_output(['mount'],text=True):
        raise ValueError('USB_1TB is not mounted')
    if path.exists():raise ValueError('Output already exists')
    if shutil.disk_usage('/Volumes/USB_1TB').free<1024**3:
        raise OSError('Need at least 1 GiB free before survey/replay')
    path.mkdir(parents=True)


def logged(command,log):
    console=Path(__file__).resolve().parent.parent/'current/progress_console.py'
    start=time.monotonic()
    with log.open('wb') as out:
        progress=subprocess.Popen(['python3',str(console),'--hide-opcode-profile'],stdin=subprocess.PIPE)
        proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        try:
            for line in proc.stdout:
                out.write(line);out.flush();progress.stdin.write(line);progress.stdin.flush()
            rc=proc.wait()
        finally:
            if proc.poll() is None:proc.terminate();proc.wait()
            progress.stdin.close();progress_rc=progress.wait()
    if rc or progress_rc:raise RuntimeError(f'Command failed ({rc}/{progress_rc}); see {log}')
    return time.monotonic()-start
