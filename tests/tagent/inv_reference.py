"""Discover and validate an immutable inv_list Gov2 reference collection."""
import csv
from pathlib import Path
import stat


def inspect_reference(path, full=False, required_levels=None, required_datasets=None):
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
        problems = []
        for row in rows:
            dataset = row['dataset'].split('.')[-1]
            expected_name = f'gov2.{dataset}.level-{int(row["compression_level"]):02d}.bv'
            if row['stream_file'] != expected_name:
                return False, ('unexpected stream_file for {} level {} in {}: expected {}, got {}'.format(
                    row['dataset'], row['compression_level'], root / 'references.tsv', expected_name, row['stream_file']))
            if required_levels is not None and int(row['compression_level']) not in required_levels:
                continue
            if required_datasets is not None and dataset not in required_datasets:
                continue
            stream = root / expected_name
            for file, field in ((stream, 'compressed_bytes'), (stream.with_suffix('.idx'), 'index_bytes'),
                                (root / 'Gov2Flat' / ('gov2.' + dataset), 'input_bytes')):
                context = '{} level {} | {} in {}'.format(row['dataset'], row['compression_level'],
                                                          field, root / 'references.tsv')
                try:
                    expected_size = int(row[field])
                except (ValueError, KeyError) as error:
                    problems.append('invalid expected size | {} | {} | {}'.format(file, context, error))
                    continue
                try:
                    info = file.stat()
                except FileNotFoundError:
                    problems.append('missing file | {} | {} | expected {} bytes'.format(file, context, expected_size))
                    continue
                except OSError as error:
                    problems.append('cannot inspect file | {} | {} | {}'.format(file, context, error))
                    continue
                if not stat.S_ISREG(info.st_mode):
                    problems.append('not a regular file | {} | {}'.format(file, context))
                elif info.st_size != expected_size:
                    problems.append('size mismatch | {} | actual {} bytes | expected {} bytes | {}'.format(
                        file, info.st_size, expected_size, context))
        if problems:
            return False, '{} reference problem(s):\n{}'.format(len(problems), '\n'.join(problems))
        return True, 'reference metadata and required streams, indexes and raw inputs verified'
    except (OSError, ValueError, KeyError, csv.Error) as error:
        return False, str(error)
