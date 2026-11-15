#!/usr/bin/env python3
"""Summarize validated full L9 preset, exact and tree Gov2 collections."""
import json
import sys
from pathlib import Path
root = Path(sys.argv[1])
for dataset in ('sorted', 'unsorted'):
    print('Gov2.' + dataset)
    rows = {}
    for mode in ('preset', 'exact', 'tree'):
        run = root / (dataset + '-' + mode)
        assert not (run / 'FAILURE.txt').exists()
        assert (run / 'report.txt').is_file()
        text = (run / 'report.txt').read_text()
        assert 'Accounting and raw-ID verification: PASS' in text
        row = json.loads((run / dataset / 'profile' / 'measurement.json').read_text())
        assert row['vectors'] == 1107205 and row['ids'] == 5979715441
        rows[mode] = row
        print(mode, 'bytes', row['size'], 'bits/ID', round(8*row['size']/row['ids'], 6),
              'encode/verify/decode s', row['timings'])
    for mode in ('preset', 'tree'):
        print(mode, 'size minus exact', rows[mode]['size']-rows['exact']['size'])
    for code in (68, 92, 98, 108, 109):
        counts = []
        for mode in ('preset', 'exact', 'tree'):
            match = next((r for r in rows[mode]['profile']['rows'] if r['code'] == code), None)
            counts.append((match['count'], match['bytes']) if match else (0, 0))
        print('opcode', code, 'preset/exact/tree count and bytes', *counts)
    print('Per-leaf regret at the wave selection stage')
    print((root / (dataset + '-leaves.tsv')).read_text())
