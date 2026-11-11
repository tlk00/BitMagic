#!/usr/bin/env python3
"""Build one optimized reference generator, separately from benchmark execution."""
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import shlex
import subprocess
import time
from gov2_reference import git_info

APP = Path(__file__).resolve().parents[1]


def main():
    """Force a release build and record its exact command, compiler and checksum."""
    machine = platform.machine().lower()
    if machine in ('arm64', 'aarch64'):
        backend, define = 'NEON', '-DBMNEONOPT'
    elif machine in ('x86_64', 'amd64'):
        flags = ''
        if Path('/proc/cpuinfo').exists():
            flags = Path('/proc/cpuinfo').read_text().lower()
        elif platform.system() == 'Darwin':
            flags = subprocess.check_output(['sysctl', '-n', 'machdep.cpu.leaf7_features'], text=True).lower()
        if 'avx2' not in flags:
            raise SystemExit('AVX2 capability could not be confirmed on this x86 host')
        backend, define = 'AVX2', '-DBMAVX2OPT'
    else:
        raise SystemExit(f'Unsupported optimized reference platform: {machine}')
    provenance = git_info()
    output = APP / ('bminv_reference_' + backend.lower())
    compiler_command = os.environ.get('CXX', 'c++')
    command = ['make', '-B', 'DEBUG=NO', f'CXX={compiler_command}', f'BMOPTFLAGS={define}', 'bminv']
    log = output.with_suffix('.build.log')
    with log.open('w') as stream:
        subprocess.run(command, cwd=APP, stdout=stream, stderr=subprocess.STDOUT, check=True)
    if git_info()['source_sha256'] != provenance['source_sha256']:
        raise SystemExit('Source changed during build; rerun build_reference.py')
    shutil.copy2(APP / 'bminv', output)
    compiler = subprocess.check_output([*shlex.split(compiler_command), '--version'], text=True)
    metadata = dict(backend=backend, architecture=machine, platform=platform.platform(),
                    command=command, compiler=compiler, build_flags=define,
                    binary_sha256=hashlib.sha256(output.read_bytes()).hexdigest(), built_unix=time.time(),
                    provenance=provenance, build_log_sha256=hashlib.sha256(log.read_bytes()).hexdigest())
    output.with_suffix('.build.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(f'Built {output}\nBuild log: {log}')


if __name__ == '__main__':
    main()
