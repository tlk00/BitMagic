#!/usr/bin/env python3
"""Render compact, single-row progress; tee retains full diagnostics in the log."""
import re
import shutil
import sys


def compact_count(value):
    count = int(value)
    if count >= 1_000_000:
        return f'{count / 1_000_000:.2f}M'
    if count >= 10_000:
        return f'{count / 1_000:.0f}k'
    return f'{count:,}'


def compact_bytes(value):
    count = int(value)
    for scale, label in ((1 << 30, 'GiB'), (1 << 20, 'MiB'), (1 << 10, 'KiB')):
        if count >= scale:
            return f'{count / scale:.1f}{label}'
    return f'{count}B'


def compact_progress(line, columns):
    match = re.match(r'^Progress (compress|encode|verify|decode)\s', line)
    if not match:
        return line[:max(columns - 1, 1)]
    phase = match.group(1)
    fields = dict(re.findall(r'([A-Za-z_][A-Za-z_0-9]*)=([^\s]+)', line))
    parts = [phase]
    if 'vectors' in fields:
        parts.append(f'{int(fields["vectors"]):,}v')
    if 'completed' in fields:
        parts.append(f'done{compact_count(fields["completed"])}')
    if 'prepared' in fields:
        parts.append(f'prep{compact_count(fields["prepared"])}')
    if 'elapsed_s' in fields:
        parts.append(f'{float(fields["elapsed_s"]):.0f}s')
    if 'output_bytes' in fields:
        parts.append(f'out{compact_bytes(fields["output_bytes"])}')
    elif 'encoded_bytes' in fields:
        parts.append(f'read{compact_bytes(fields["encoded_bytes"])}')
    options = []
    if 'bwt_rejected' in fields:
        options.append(f'BWT-reject{compact_count(fields["bwt_rejected"])}')
    if 'evaluated' in fields:
        options.append(f'eval{compact_count(fields["evaluated"])}')
    elif 'evaluations' in fields:
        options.append(f'eval{compact_count(fields["evaluations"])}')
    if 'saved_bytes' in fields:
        options.append(f'saved{compact_bytes(fields["saved_bytes"])}')
    if 'vectors_per_s' in fields:
        options.append(f'{float(fields["vectors_per_s"]):.1f}v/s')
    if 'accepted' in fields:
        options.append(f'accepted{compact_count(fields["accepted"])}')
    if 'trials' in fields:
        options.append(f'trials{compact_count(fields["trials"])}')
    if 'input_bytes' in fields:
        options.append(f'in{compact_bytes(fields["input_bytes"])}')
    limit = max(columns - 1, 1)  # stay one column short to prevent terminal wrap
    result = ' '.join(parts)
    for item in options:
        candidate = result + ' ' + item
        if len(candidate) <= limit:
            result = candidate
    return result[:limit]


hide_opcode_profile = '--hide-opcode-profile' in sys.argv[1:]
in_opcode_profile = False
interactive = sys.stdout.isatty()
pending = None
active = False
for raw in sys.stdin.buffer:
    line = raw.decode('utf-8', errors='replace').rstrip('\r\n')
    if hide_opcode_profile and line.startswith(('Opcode profile blobs=', 'Opcode cardinality ids=')):
        in_opcode_profile = True
        continue
    if in_opcode_profile:
        if line.startswith('Opcode profile scope='):
            in_opcode_profile = False
            continue
        if (line.startswith(('Opcode group ', 'Opcode ID ', 'Opcode efficiency:')) or
            re.match(r'^.+\s+\d+\s+\d+\s+\d+\s+\d+\s+[\d.]+$', line) or
            re.match(r'^\s*\d+\s+\d+\s+(?:[\d.]+|n/a)\s+(?:[\d.]+|n/a)$', line)):
            continue
        # Keep unexpected messages/errors visible even if a table is interrupted.
        in_opcode_profile = False
    if line.startswith(('Selection block ', 'BWT block ', 'Periodic block ',
                        'Horizontal block ')):
        continue
    if line.startswith('Progress '):
        width = shutil.get_terminal_size(fallback=(80, 24)).columns if interactive else 80
        pending = compact_progress(line, width)
        if interactive:
            sys.stdout.write('\r\x1b[2K' + pending)
            sys.stdout.flush()
            active = True
        else:
            sys.stdout.write(pending + '\n')
            sys.stdout.flush()
            pending = None
        continue
    if interactive and active:
        sys.stdout.write('\r\x1b[2K')
        active = False
    pending = None
    sys.stdout.write(line + '\n')
    sys.stdout.flush()
if interactive and active:
    sys.stdout.write('\r\x1b[2K')
if pending:
    sys.stdout.write(pending + '\n')
sys.stdout.flush()
