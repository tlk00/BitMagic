#!/usr/bin/env python3
"""Sample paired Gov2 posting-list ordinals and report horizontal XOR trials."""
import argparse
from array import array
from collections import Counter
from datetime import datetime
import mmap
from pathlib import Path
import random
import struct
import subprocess
import sys


def offsets(data):
    """Index native little-endian length-prefixed records and validate boundaries."""
    starts = array("Q")
    pos = 0
    while pos < len(data):
        if pos + 4 > len(data):
            raise ValueError(f"truncated record header at {pos}")
        count, = struct.unpack_from("<I", data, pos)
        end = pos + 4 + count * 4
        if not count or end > len(data):
            raise ValueError(f"invalid record at {pos}, count={count}")
        starts.append(pos)
        pos = end
    starts.append(pos)
    return starts


def run_dataset(args, name, expected_total):
    """Scan one collection and feed sampled records directly to the survey."""
    source = args.input_dir / f"gov2.{name}"
    print(f"Indexing {source}", flush=True)
    with source.open("rb") as raw, mmap.mmap(raw.fileno(), 0, access=mmap.ACCESS_READ) as data:
        starts = offsets(data)
        total = len(starts) - 1
        if expected_total is not None and total != expected_total:
            raise ValueError("paired collections have different vector counts")
        if args.sample > total:
            raise ValueError(f"sample {args.sample} exceeds collection size {total}")
        rng = random.Random(args.seed)
        selected = [rng.randrange(i * total // args.sample, (i + 1) * total // args.sample)
                    for i in range(args.sample)]
        (args.out / f"{name}.ordinals.txt").write_text("".join(f"{i}\n" for i in selected))
        print(f"{name}: sampling {args.sample:,} of {total:,} vectors", flush=True)
        metrics = args.out / f"{name}.tsv"
        with metrics.open("w") as results:
            results.write("metric\tvalue\n")
            results.flush()
            with subprocess.Popen([args.bin.resolve()], stdin=subprocess.PIPE, stdout=results) as child:
                try:
                    with child.stdin as stream:
                        for i in selected:
                            stream.write(struct.pack("<I", i))
                            stream.write(data[starts[i]:starts[i + 1]])
                except BrokenPipeError:
                    pass
                if child.wait():
                    raise ValueError(f"{name} survey exited with {child.returncode}")
    counts = Counter()
    for line in metrics.read_text().splitlines()[1:]:
        key, value = line.split("\t")
        counts[key] = int(value)
    if counts["vectors"] != args.sample:
        raise ValueError(f"incomplete results for {name}")
    return total, counts


def main():
    """Run both datasets and render a weighted summary from measured byte totals."""
    usb = Path("/Volumes/USB_1TB")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=usb / "int_set_bench/Gov2Flat")
    parser.add_argument("--bin", type=Path, default=Path("./gov2_posting_xor"))
    parser.add_argument("--sample", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261111)
    parser.add_argument("--out", type=Path, default=usb / "BitMagic-experiments" /
                        ("gov2-horizontal-survey-" + datetime.now().strftime("%Y%m%d-%H%M%S")))
    args = parser.parse_args()
    if args.sample < 1 or sys.byteorder != "little":
        parser.error("a positive sample and little-endian host are required")
    if usb in args.out.resolve().parents and not usb.is_mount():
        parser.error(f"{usb} is not mounted")
    for name in ("sorted", "unsorted"):
        if not (args.input_dir / f"gov2.{name}").is_file():
            parser.error(f"missing gov2.{name} in {args.input_dir}")
    if not args.bin.is_file():
        parser.error(f"build {args.bin} with make gov2-posting-xor")
    args.out.mkdir(parents=True, exist_ok=False)
    print(f"Reports: {args.out}", flush=True)
    total = None
    rows = []
    for name in ("sorted", "unsorted"):
        total, c = run_dataset(args, name, total)
        for scorer in ("pop", "run", "union"):
            gain = c[scorer + "_group_saved_bytes"]
            rows.append([name, scorer, c["vectors"], c["pairs"], c["isolated_level10_bytes"],
                         c[scorer + "_block_saved_bytes"], c[scorer + "_winning_groups"], gain,
                         f'{100 * gain / c["isolated_level10_bytes"]:.5f}'])
    header = "dataset scorer vectors pairs baseline_bytes block_saved_B winning_groups net_saved_B net_pct"
    report = "\n".join(["\t".join(header.split())] + ["\t".join(map(str, r)) for r in rows])
    (args.out / "report.tsv").write_text(report + "\n")
    print(report)
    print(f"Reports: {args.out}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        sys.exit(f"Survey failed: {exc}")
