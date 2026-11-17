#!/usr/bin/env python3
"""Compare incumbent and four-recipe XOR selection on a prepared sparse vector."""
import argparse
import csv
import pathlib
import re
import statistics
import subprocess


def main():
    """Run paired repetitions, verify baseline bytes and write logs and summaries."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--program', required=True, type=pathlib.Path)
    parser.add_argument('--input', required=True, type=pathlib.Path)
    parser.add_argument('--output', required=True, type=pathlib.Path)
    parser.add_argument('--reference', type=pathlib.Path,
                        help='Previous NED matrix directory for baseline byte comparisons')
    parser.add_argument('--levels', type=int, nargs='+', default=[6, 7, 8])
    parser.add_argument('--repeats', type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1 or any(level < 0 or level > 10 for level in args.levels):
        parser.error('Positive repeats and levels from 0 to 10 are required')
    program = args.program.resolve(strict=True)
    source = args.input.resolve(strict=True)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    measurement = re.compile(r'Serialization measurement bytes=(\d+) serialize_ms=([\d.eE+-]+) deserialize_ms=([\d.eE+-]+)')
    counter = re.compile(r'XOR trials blocks=(\d+) alternatives=(\d+) changed=(\d+) saved_bytes=(\d+)')
    fields = ['level', 'mode', 'repeat', 'bytes', 'serialize_ms', 'deserialize_ms',
              'blocks', 'alternatives', 'changed', 'saved_bytes']
    results = []
    with (output / 'results.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for level in args.levels:
            for repeat in range(1, args.repeats + 1):
                # Alternate pair order to reduce consistent order bias.
                modes = ['current', 'four'] if repeat % 2 else ['four', 'current']
                for mode in modes:
                    label = f'level{level}_{mode}'
                    blob = output / f'{label}.sv'
                    command = [str(program), '-svin', str(source), '-svout', str(blob),
                               '-level', str(level), '-xor']
                    if mode == 'four':
                        command.append('-xor-candidates4')
                    if repeat == 1:
                        command.append('-profile-opcodes')
                    print(f'{label} repeat={repeat}', flush=True)
                    run = subprocess.run(command, stdout=subprocess.PIPE,
                                         stderr=subprocess.STDOUT, text=True)
                    (output / f'{label}_repeat{repeat}.log').write_text(run.stdout)
                    match = measurement.search(run.stdout)
                    if run.returncode or not match:
                        raise RuntimeError(f'{label} repeat {repeat} failed')
                    if repeat == 1 and 'reconciled=1' not in run.stdout:
                        raise RuntimeError(f'{label}: missing reconciled opcode profile')
                    if mode == 'current' and args.reference:
                        reference = args.reference / f'remap0_xor1_level{level}.sv'
                        if blob.read_bytes() != reference.read_bytes():
                            raise RuntimeError(f'{label}: default differs from prior matrix')
                    row = dict(zip(fields[:6], [level, mode, repeat, *match.groups()]))
                    stats = counter.search(run.stdout)
                    if mode == 'four' and not stats:
                        raise RuntimeError(f'{label}: missing experiment counters')
                    row.update(zip(fields[6:], stats.groups() if stats else ['0'] * 4))
                    writer.writerow(row)
                    stream.flush()
                    results.append(row)
    with (output / 'summary.csv').open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['level', 'mode', 'bytes', 'median_serialize_ms', 'median_deserialize_ms',
                         'blocks', 'alternatives', 'changed', 'saved_bytes'])
        for level in args.levels:
            for mode in ['current', 'four']:
                rows = [r for r in results if r['level'] == level and r['mode'] == mode]
                if len({r['bytes'] for r in rows}) != 1:
                    raise RuntimeError('Serialized sizes changed between repetitions')
                writer.writerow([level, mode, rows[0]['bytes'],
                                 statistics.median(float(r['serialize_ms']) for r in rows),
                                 statistics.median(float(r['deserialize_ms']) for r in rows),
                                 *[rows[0][key] for key in fields[6:]]])
    print(f'Complete: {output}', flush=True)


if __name__ == '__main__':
    main()
