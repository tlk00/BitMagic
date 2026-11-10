#!/usr/bin/env python3
"""Summarize only complete, emitted and verified Gov2 collection streams."""
from pathlib import Path
import re
import sys


def require(condition, message):
    if not condition:
        raise SystemExit(message)


out, bv_stream, input_name, requested = map(str, sys.argv[1:5])
out = Path(out)
bv_stream = Path(bv_stream)
require(bv_stream.is_file(), f"Missing BV reference: {bv_stream}")
bv_bytes = bv_stream.stat().st_size
bv_index = bv_stream.with_suffix(".idx")
require(bv_index.is_file(), f"Missing BV reference index: {bv_index}")
bv_index_bytes = bv_index.stat().st_size
rows = {}
for repr_name in ("raw", "gaps"):
    for codec in ("plain", "xor", "select"):
        name = f"{repr_name}.{codec}"
        stream = out / f"{name}.sv"
        index = out / f"{name}.idx"
        require(stream.is_file() and index.is_file(), f"Missing {name} stream/index")
        results = {}
        for phase in ("compress", "verify", "decode"):
            path = out / f"{name}.{phase}.log"
            require(path.is_file(), f"Missing phase log: {path}")
            lines = re.findall(r"^Result (.+)$", path.read_text(errors="replace"), re.M)
            require(len(lines) == 1, f"Missing or repeated result in {path}")
            fields = dict(part.split("=", 1) for part in lines[0].split() if "=" in part)
            require(fields.get("action") == {"compress": "encode", "verify": "verify", "decode": "decode"}[phase],
                    f"Wrong action in {path}")
            results[phase] = fields
        size = stream.stat().st_size
        counts = {int(result["vectors"]) for result in results.values()}
        require(len(counts) == 1 and next(iter(counts)) > 0, f"Vector count mismatch for {name}")
        if requested != "0":
            require(next(iter(counts)) == int(requested), f"Requested vector count mismatch for {name}")
        require(all(int(result["stream_bytes"]) == size for result in results.values()),
                f"Complete stream size mismatch for {name}")
        rows[name] = (size, results)

counts = {int(results["compress"]["vectors"]) for _, results in rows.values()}
ints = {int(results["compress"]["original_ints"]) for _, results in rows.values()}
require(len(counts) == 1 and len(ints) == 1, "Input counts differ between variants")
count = next(iter(counts))
original_ints = next(iter(ints))
require(count > 0 and original_ints > 0, "Empty comparison")

if bv_stream.parent == out:
    for phase, needle in (("compress", "Total vectors="),
                          ("verify", "Verification complete."),
                          ("decode", "Decode complete.")):
        path = out / f"bv.{phase}.log"
        require(path.is_file() and needle in path.read_text(errors="replace"),
                f"BV {phase} did not finish")
    source = out / "bv.compress.log"
    match = re.search(r"^Total vectors=(\d+)$", source.read_text(errors="replace"), re.M)
    require(match and int(match.group(1)) == count, "BV reference count differs")

bv_prefix = "bv" if bv_stream.parent == out else "auto"
bv_times = []
for phase in ("compress", "verify", "decode"):
    path = bv_stream.parent / f"{bv_prefix}.{phase}.time"
    if path.is_file():
        values = path.read_text(errors="replace").split()
        bv_times.append(float(values[0]) if values else None)
    else:
        bv_times.append(None)
def shown_time(value):
    return f"{value:8.2f}" if value is not None else f"{'—':>8}"

print(f"Correctness: PASS; all six SV streams verified against {input_name} and decoded ({count:,} vectors).")
print(f"Original IDs: {original_ints:,}; BV reference: {bv_bytes:,} stream bytes, {bv_index_bytes:,} index bytes.")
print(f"BV reference stream: {bv_stream}")
print("Bits/ID uses all original IDs as the denominator; stored bits/ID includes the index.")
print("Variant       Stream bytes    Index bytes   Bits/ID   Stored/ID  Save vs BV  Save vs SV plain  Encode s  Verify s  Decode s")
print("------------  --------------  -----------  --------  ----------  ----------  ----------------  --------  --------  --------")
print(f"{'BV reference':12}  {bv_bytes:14,}  {bv_index_bytes:11,}  "
      f"{8*bv_bytes/original_ints:8.4f}  {8*(bv_bytes+bv_index_bytes)/original_ints:10.4f}  "
      f"{'—':>10}  {'—':>16}  " + "  ".join(shown_time(t) for t in bv_times))
for repr_name in ("raw", "gaps"):
    for codec in ("plain", "xor", "select"):
        name = f"{repr_name}.{codec}"
        size, phases = rows[name]
        times = [float(phases[p]["elapsed_s"]) for p in ("compress", "verify", "decode")]
        plain = rows[f"{repr_name}.plain"][0]
        index_size = (out / f"{name}.idx").stat().st_size
        print(f"{name:12}  {size:14,}  {index_size:11,}  "
              f"{8*size/original_ints:8.4f}  {8*(size+index_size)/original_ints:10.4f}  "
              f"{bv_bytes-size:+10,}  {plain-size:+16,}  "
              f"{times[0]:8.2f}  {times[1]:8.2f}  {times[2]:8.2f}")
    fields = rows[f"{repr_name}.select"][1]["compress"]
    xor_fields = rows[f"{repr_name}.xor"][1]["compress"]
    if "xor_blocks" in xor_fields:
        print(f"{repr_name} XOR-selected blocks: {int(xor_fields['xor_blocks']):,}.")
    print(f"{repr_name} selected BIT-block winners: " + ", ".join(
        f"{key}={fields.get(key, '?')}" for key in
        ("ordinary", "periodic", "bwt", "composed", "step2")))
best_name, (best_size, _) = min(rows.items(), key=lambda item: item[1][0])
print(f"Smallest SV stream: {best_name}, {best_size:,} bytes; "
      f"difference from BV reference: {best_size-bv_bytes:+,} bytes "
      f"({100*(best_size/bv_bytes-1):+.4f}%).")
print("Savings columns are bytes; negative values mean the SV stream is larger.")
