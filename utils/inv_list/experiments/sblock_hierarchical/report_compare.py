#!/usr/bin/env python3
"""Report validated 108-only, both-format and 109-only full-collection runs."""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("control", type=Path)
parser.add_argument("both", type=Path)
parser.add_argument("hierarchical", type=Path)
args = parser.parse_args()
runs = [args.control, args.both, args.hierarchical]
for run in runs:
    if not (run / "verification-report.txt").is_file():
        raise SystemExit(f"Missing completed validation report in {run}")
    if (run / "FAILURE.txt").exists():
        raise SystemExit(f"Run failed in {run}")

for dataset in ("sorted", "unsorted"):
    print(f"Gov2.{dataset}")
    print()
    print("| Configuration | Bytes | Bits/ID | Encode s | Verify s | Decode s | Saved vs 108-only |")
    print("|---|---:|---:|---:|---:|---:|---:|")
    records = []
    for label, run in zip(("108 only", "108 + 109", "109 only"), runs):
        data = json.loads((run / dataset / "profile" / "measurement.json").read_text())
        if data["vectors"] != 1107205 or data["ids"] != 5979715441:
            raise SystemExit(f"Incomplete collection in {run}/{dataset}")
        if records and (data["ids"], data["vectors"]) != (records[0]["ids"], records[0]["vectors"]):
            raise SystemExit("Validation counts differ")
        records.append(data)
        baseline = records[0]["size"]
        times = data["timings"]
        print(f"| {label} | {data['size']:,} | {8*data['size']/data['ids']:.6f} | "
              f"{times[0]:.3f} | {times[1]:.3f} | {times[2]:.3f} | {baseline-data['size']:+,} |")
    print(f"Retiring 108 changes size by {records[2]['size']-records[1]['size']:+,} bytes "
          "relative to allowing both formats.")
    print()
    print("| Opcode | Name | 108-only bytes/count | Both bytes/count | 109-only bytes/count |")
    print("|---:|---|---:|---:|---:|")
    profiles = [{r["code"]: r for r in d["profile"]["rows"]} for d in records]
    for code in sorted(set().union(*(set(p) for p in profiles))):
        rows = [p.get(code, {}) for p in profiles]
        cells = [f"{r.get('bytes', 0):,} / {r.get('count', 0):,}" for r in rows]
        name = next(r["name"] for r in rows if r)
        print(f"| {code} | {name} | " + " | ".join(cells) + " |")
    if profiles[0].get(109, {}).get("count", 0) or profiles[2].get(108, {}).get("count", 0):
        raise SystemExit("Excluded opcode was emitted")
    print()
