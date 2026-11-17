#!/usr/bin/env python3
"""Run serialization level and XOR comparisons on two prepared NED vectors."""
import argparse
import csv
import datetime
import pathlib
import re
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--program', required=True, type=pathlib.Path)
    parser.add_argument('--plain', required=True, type=pathlib.Path)
    parser.add_argument('--remap', required=True, type=pathlib.Path)
    parser.add_argument('--output', required=True, type=pathlib.Path)
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('--repeats must be positive')
    program = args.program.resolve(strict=True)
    sources = [args.plain.resolve(strict=True), args.remap.resolve(strict=True)]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    pattern = re.compile(r'Serialization measurement bytes=(\d+) serialize_ms=([\d.eE+-]+) deserialize_ms=([\d.eE+-]+)')
    fields = ['remap', 'xor', 'level', 'repeat', 'bytes', 'serialize_ms', 'deserialize_ms']
    with (output / 'results.csv').open('w', newline='') as result:
        writer = csv.DictWriter(result, fieldnames=fields)
        writer.writeheader()
        for remap, source in enumerate(sources):
            for xor in (0, 1):
                for level in range(6, 11):
                    label = f'remap{remap}_xor{xor}_level{level}'
                    sizes = set()
                    for repeat in range(1, args.repeats + 1):
                        command = [str(program), '-svin', str(source), '-svout',
                                   str(output / (label + '.sv')), '-level', str(level)]
                        if xor:
                            command.append('-xor')
                        if repeat == 1:
                            command.append('-profile-opcodes')
                        print(f'{datetime.datetime.now().isoformat(timespec="seconds")} {label} repeat={repeat}', flush=True)
                        completed = subprocess.run(command, stdout=subprocess.PIPE,
                                                   stderr=subprocess.STDOUT, text=True)
                        (output / f'{label}_repeat{repeat}.log').write_text(completed.stdout)
                        if completed.returncode:
                            raise RuntimeError(f'{label} repeat {repeat} failed: {completed.returncode}')
                        match = pattern.search(completed.stdout)
                        if not match:
                            raise RuntimeError(f'{label} repeat {repeat}: missing measurement')
                        size, encode, decode = match.groups()
                        sizes.add(size)
                        if len(sizes) != 1:
                            raise RuntimeError(f'{label}: serialized size changed between repetitions')
                        writer.writerow(dict(remap=remap, xor=xor, level=level, repeat=repeat,
                                             bytes=size, serialize_ms=encode, deserialize_ms=decode))
                        result.flush()
    print(f'Complete: {output}', flush=True)


if __name__ == '__main__':
    try:
        main()
    except (OSError, RuntimeError) as exc:
        print(f'Error: {exc}', file=sys.stderr)
        sys.exit(1)
