#!/usr/bin/env python3
"""Profile actual isolated level-10 groups against constituent block records on Gov2.sorted."""
import argparse
import csv
from collections import Counter, defaultdict
from datetime import datetime
import mmap
from pathlib import Path
import random
import struct
import subprocess
import sys

from run_posting_xor import offsets


def expand(original, total, sample):
    """Retain each original draw and sample the other subdivisions of its stratum."""
    n = len(original)
    if not n or sample < n or sample > total or sample % n:
        raise ValueError("sample must be a positive multiple of the original sample within collection size")
    rng = random.Random(20261111)
    selected = []
    for i, anchor in enumerate(original):
        lo, hi = i * total // n, (i + 1) * total // n
        if not lo <= anchor < hi:
            raise ValueError(f"original ordinal {anchor} outside expected stratum {i}")
        subdivisions = sample // n
        for j in range(subdivisions):
            a = lo + j * (hi - lo) // subdivisions
            b = lo + (j + 1) * (hi - lo) // subdivisions
            selected.append(anchor if a <= anchor < b else rng.randrange(a, b))
    return selected


def report(out, selected, original, previous):
    """Reconcile the replay and summarize group outcomes by route and locality."""
    stats = defaultdict(Counter)
    codecs = defaultdict(Counter)
    replay = Counter()
    seen = set()
    gains = []
    with (out / "groups.tsv").open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            vec = int(row["ordinal"])
            seen.add(vec)
            cohort = "original" if vec in original else "additional"
            base, alt = int(row["baseline_B"]), int(row["constituent_B"])
            gain = base - alt
            bit, gap = int(row["bit"]), int(row["gap"])
            ids, runs = int(row["ids"]), int(row["one_runs"])
            code = int(row["group_opcode"])
            if cohort == "original":
                replay["ids"] += ids
                replay["groups"] += 1
                replay["isolated_level10_bytes"] += base
                if bit >= 2:
                    replay["ungrouping_only_saved_bytes"] += max(0, gain)
            if code:
                outcome = "win" if gain > 0 else "tie" if gain == 0 else "loss"
                for record in filter(None, row["alt_records"].split(",")):
                    _, name, count, size = record.split(":")
                    for label in (cohort, "combined"):
                        codecs[label, code, outcome, name].update(records=int(count), bytes=int(size))
            if gain > 0:
                gains.append((gain, cohort, row))
            # Keep direct routes visible in totals. Breakdown rows characterize
            # the groups that actually emitted a superblock opcode.
            breakdown = [("all", "all"), ("opcode", str(code) if code else "individual")]
            if code:
                breakdown += [
                    ("block_types", "mixed" if bit and gap else "has_BIT" if bit else "no_BIT"),
                    ("mean_one_run", "1" if ids == runs else "<4" if ids < 4*runs else "4+"),
                    ("largest_slot_share", ">=50%" if 2*int(row["max_slot_ids"]) >= ids else "<50%"),
                ]
            for population in ("combined", cohort):
                for dimension, category in breakdown:
                    s = stats[population, dimension, category]
                    s.update(groups=1, baseline_B=base, wins=int(gain > 0),
                             ties=int(gain == 0), losses=int(gain < 0),
                             saved_if_selective_B=max(0, gain), cost_if_forced_B=max(0, -gain))
    if seen != set(selected):
        raise ValueError("profiled ordinals differ from selected ordinals")
    for key in ("groups", "ids", "isolated_level10_bytes", "ungrouping_only_saved_bytes"):
        expected = previous.get(key, 0) if key == "ungrouping_only_saved_bytes" else previous[key]
        if replay[key] != expected:
            raise ValueError(f"10K replay disagrees on {key}: {replay[key]} vs {expected}")
    with (out / "winners.tsv").open("w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["saved_B", "cohort", *reader.fieldnames])
        for gain, cohort, row in sorted(gains, key=lambda x: -x[0]):
            writer.writerow([gain, cohort, *row.values()])
    with (out / "constituent_opcodes.tsv").open("w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["cohort", "group_opcode", "individual_outcome", "constituent_opcode",
                         "record_count", "record_bytes"])
        for key, s in sorted(codecs.items()):
            writer.writerow([*key, s["records"], s["bytes"]])
    columns = "groups wins ties losses baseline_B saved_if_selective_B cost_if_forced_B".split()
    with (out / "report.tsv").open("w") as f:
        writer = csv.writer(f, delimiter="\t", lineterminator="\n")
        writer.writerow(["cohort", "dimension", "category", *columns, "weighted_saved_pct"])
        for key, s in sorted(stats.items()):
            writer.writerow([*key, *(s[k] for k in columns),
                             f'{100*s["saved_if_selective_B"]/s["baseline_B"]:.6f}'])
    print("cohort        vectors   groups    wins     saved_B   baseline_B")
    for label, n in (("original", len(original)), ("additional", len(selected)-len(original)),
                     ("combined", len(selected))):
        s = stats[label, "all", "all"]
        if n:
            print(f'{label:12} {n:8} {s["groups"]:8} {s["wins"]:7} '
                  f'{s["saved_if_selective_B"]:11} {s["baseline_B"]:12}')
    print(f'Prior cohort reconciled: {replay["ungrouping_only_saved_bytes"]:,} B in groups with >=2 BIT blocks.')


def main():
    """Replay saved ordinals and pipe sampled raw records to the in-memory profiler."""
    usb = Path("/Volumes/USB_1TB")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=usb / "int_set_bench/Gov2Flat/gov2.sorted")
    parser.add_argument("--prior", type=Path, required=True, help="directory of the completed XOR survey")
    parser.add_argument("--bin", type=Path, default=Path("./gov2_posting_xor"))
    parser.add_argument("--sample", type=int, default=50000)
    parser.add_argument("--out", type=Path, default=usb / "BitMagic-experiments" /
                        ("gov2-sorted-group-profile-" + datetime.now().strftime("%Y%m%d-%H%M%S")))
    args = parser.parse_args()
    if sys.byteorder != "little" or (usb in args.out.resolve().parents and not usb.is_mount()):
        parser.error("requires a little-endian host and a mounted output volume")
    original = list(map(int, (args.prior / "sorted.ordinals.txt").read_text().splitlines()))
    with (args.prior / "sorted.tsv").open() as f:
        previous = {r["metric"]: int(r["value"]) for r in csv.DictReader(f, delimiter="\t")}
    if len(original) != previous["vectors"]:
        raise ValueError("prior report and ordinals disagree")
    with args.input.open("rb") as f, mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as data:
        starts = offsets(data)
        selected = expand(original, len(starts)-1, args.sample)
        args.out.mkdir(parents=True, exist_ok=False)
        print(f"Profiling {len(selected):,} vectors. Reports: {args.out}", flush=True)
        (args.out / "ordinals.txt").write_text("".join(f"{i}\n" for i in selected))
        with (args.out / "groups.tsv").open("wb") as rows, subprocess.Popen(
                [args.bin.resolve(), "--group-profile"], stdin=subprocess.PIPE, stdout=rows) as child:
            try:
                with child.stdin as stream:
                    for i in selected:
                        stream.write(struct.pack("<I", i))
                        stream.write(data[starts[i]:starts[i+1]])
            except BrokenPipeError:
                pass
            if child.wait():
                raise ValueError(f"group profiler exited with {child.returncode}")
    report(args.out, selected, set(original), previous)
    print(f"Reports: {args.out}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError) as error:
        sys.exit(f"Profile failed: {error}")
