#!/usr/bin/env python3
"""Report complete L8 tree versus exact stream sizes, times and opcode counts."""
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
for dataset in ('sorted', 'unsorted'):
    print('Gov2.' + dataset)
    data = []
    for mode in ('exact', 'tree'):
        run = root / (dataset + '-' + mode)
        assert not (run / 'FAILURE.txt').exists()
        assert (run / 'report.txt').is_file()
        row = json.loads((run / dataset / 'profile' / 'measurement.json').read_text())
        assert row['vectors'] == 1107205 and row['ids'] == 5979715441
        data.append(row)
        print(mode, 'bytes', row['size'], 'bpi', round(8*row['size']/row['ids'], 6),
              'encode/verify/decode s', row['timings'])
    print('regret bytes', data[1]['size']-data[0]['size'])
    for code in (68, 92, 98, 108, 109):
        counts = []
        for row in data:
            match = next((r for r in row['profile']['rows'] if r['code'] == code), None)
            counts.append((match['count'], match['bytes']) if match else (0, 0))
        print('opcode', code, 'exact/tree count and bytes', *counts)
    print((root / (dataset + '-leaves.tsv')).read_text())
