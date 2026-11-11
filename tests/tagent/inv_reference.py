"""Discover and validate an immutable inv_list Gov2 reference collection."""
import csv
from pathlib import Path


def inspect_reference(path, full=False, required_levels=None):
    """Check matrix identity and file sizes; return a human-readable result."""
    root = Path(path).expanduser().resolve()
    try:
        if (root / '.reference-lock').exists():
            return False, 'reference generation lock exists'
        with (root / 'references.tsv').open(newline='') as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        expected = {(f'Gov2.{d}', n) for d in ('sorted', 'unsorted') for n in range(1, 11)}
        if len(rows) != 20 or {(r['dataset'], int(r['compression_level'])) for r in rows} != expected:
            return False, 'expected 20 sorted/unsorted reference cases'
        if full and not (root / 'matrix.json').is_file():
            return False, 'matrix.json is required for pedantic mode'
        for row in rows:
            dataset = row['dataset'].split('.')[-1]
            expected_name = f'gov2.{dataset}.level-{int(row["compression_level"]):02d}.bv'
            if row['stream_file'] != expected_name:
                return False, 'unexpected reference filename'
            if required_levels is not None and int(row['compression_level']) not in required_levels:
                continue
            stream = root / expected_name
            for file, field in ((stream, 'compressed_bytes'), (stream.with_suffix('.idx'), 'index_bytes'),
                                (root / 'Gov2Flat' / ('gov2.' + dataset), 'input_bytes')):
                if not file.is_file() or file.stat().st_size != int(row[field]):
                    return False, 'missing or incomplete file: ' + str(file)
        return True, 'reference metadata and required streams, indexes and raw inputs verified'
    except (OSError, ValueError, KeyError, csv.Error) as error:
        return False, str(error)
