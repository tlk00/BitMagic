#!/usr/bin/env python3
"""Summarize exact 108/109 search trials and measure simple gate regret."""
import argparse
from collections import Counter, defaultdict
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("scores", type=Path)
args = parser.parse_args()
thresholds = (1, 2, 4, 8, 16, 32, 64, 128, 192, 224, 240, 248, 252, 255, 256)
counts = Counter()
widths = defaultdict(Counter)
blocks = defaultdict(Counter)
costs = Counter()
by_outcome = defaultdict(Counter)
gates = {t: Counter() for t in thresholds}
map_modes = Counter()
groups = 0
with args.scores.open() as source:
    header = next(source).rstrip().split("\t")
    index = {name: i for i, name in enumerate(header)}
    for line in source:
        fields = line.rstrip().split("\t")
        if len(fields) != len(header):
            raise SystemExit("Incomplete score row")
        scope = int(fields[1])
        flat, hier = int(fields[6]), int(fields[7])
        result = "unavailable108" if not flat else "108" if flat < hier else "109" if hier < flat else "tie"
        widths[scope][result] += 1
        if scope != 6:
            continue
        row = {name: int(fields[i]) for name, i in index.items()}
        if not flat:
            raise SystemExit("Overall flat score unexpectedly unavailable")
        groups += 1
        if row["group"] != groups - 1:
            raise SystemExit("Missing or reordered observation groups")
        rival, occupied = row["rival"], row["blocks"]
        chosen = min(flat, hier, rival)
        useful = min(flat, hier) < rival
        counts[result] += 1
        counts["wave_useful"] += useful
        counts["useful_" + result] += useful
        blocks[occupied][result] += 1
        blocks[occupied]["useful_" + result] += useful
        for key in ("flat_map_bits", "hier_map_bits", "flat_header_bits", "hier_header_bits",
                    "flat_payload_bits", "hier_payload_bits", "flat_padding_bits", "hier_padding_bits"):
            costs[key] += row[key]
            by_outcome[result][key] += row[key]
        for prefix, size in (("flat", flat), ("hier", hier)):
            total = sum(row[prefix + "_" + part + "_bits"] for part in ("map", "header", "payload", "padding"))
            if total + 8 != 8 * size:
                raise SystemExit("Component costs do not reconcile")
        map_modes[row["hier_map_mode"]] += 1
        counts["dictionary_best_" + result] += row["hier_mode"] == 3
        counts["same_width_payload_" + result] += (
            row["flat_width"] == row["hier_width"] and
            row["flat_payload_bits"] == row["hier_payload_bits"])
        for threshold, data in gates.items():
            if occupied >= threshold:
                data["skip109_calls"] += 1
                loss = min(rival, flat) - chosen
                data["skip109_lost_bytes"] += loss
                data["skip109_errors"] += loss > 0
                data["skip109_max"] = max(data["skip109_max"], loss)
            if occupied <= threshold:
                data["skip108_calls"] += 1
                loss = min(rival, hier) - chosen
                data["skip108_lost_bytes"] += loss
                data["skip108_errors"] += loss > 0
                data["skip108_max"] = max(data["skip108_max"], loss)
if any(sum(widths[scope].values()) != groups for scope in range(7)):
    raise SystemExit("Missing width rows")
print("Exact search trials; costs are pre-emission and higher-count opcode 98 may later replace records.")
print("Strict comparisons use complete padded bytes. Ties are reported separately.")
print("A wave trial is useful when its best byte size strictly beats the ordinary/adaptive incumbent.")
print(f"Overall groups={groups:,}; counts={dict(counts)}")
print("Selected best-width component sums in bits, before incumbent fallback")
for key, value in costs.items():
    print(key, value)
print("Selected best-width costs by outcome (bits)")
for outcome, data in sorted(by_outcome.items()):
    print(outcome, dict(data))
print("Comparison per width; scope 0..5 means 128..4096 bits; scope 6 selects each format's best width")
for scope, data in sorted(widths.items()):
    print(scope, dict(data))
print("Occupied source blocks: strict 108, strict 109, ties, useful 108, useful 109, useful ties")
for occupied, data in sorted(blocks.items()):
    print(occupied, *(data[key] for key in ("108", "109", "tie", "useful_108", "useful_109", "useful_tie")))
print("Gate regret against exact selection with ordinary/adaptive fallback, summed in bytes")
print("threshold skip109_if_blocks_ge_calls errors lost_bytes max_loss skip108_if_blocks_le_calls errors lost_bytes max_loss")
for threshold, data in gates.items():
    print(threshold, data["skip109_calls"], data["skip109_errors"], data["skip109_lost_bytes"], data["skip109_max"],
          data["skip108_calls"], data["skip108_errors"], data["skip108_lost_bytes"], data["skip108_max"])
print("109 map coders: root*4+wave, each 0=all, 1=gamma, 2=BIC occupied, 3=BIC empty")
for code, count in sorted(map_modes.items()):
    print(code, count)
