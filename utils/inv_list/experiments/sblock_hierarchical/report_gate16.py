#!/usr/bin/env python3
"""Compare validated full L8 exact and <=16 occupied-block gate collections."""
import argparse
import json
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("root", type=Path)
args = parser.parse_args()
for mode in ("exact", "gate16"):
    run = args.root / mode
    if (run / "FAILURE.txt").exists() or not (run / "verification-report.txt").is_file():
        raise SystemExit(f"Missing validated run: {run}")
for dataset in ("sorted", "unsorted"):
    data = []
    print(f"Gov2.{dataset}\n")
    print("| Selection | BV bytes | Bits/ID | Encode s | Verify s | Decode s |")
    print("|---|---:|---:|---:|---:|---:|")
    for mode in ("exact", "gate16"):
        path = args.root / mode / dataset / "profile" / "measurement.json"
        row = json.loads(path.read_text())
        if row["vectors"] != 1107205 or row["ids"] != 5979715441:
            raise SystemExit(f"Incomplete collection: {path}")
        data.append(row)
        print(f"| {mode} | {row['size']:,} | {8*row['size']/row['ids']:.6f} | "
              + " | ".join(f"{t:.3f}" for t in row["timings"]) + " |")
    exact, gate = data
    delta = gate["size"] - exact["size"]
    print(f"\nGate size loss relative to exact selection: {delta:+,} bytes "
          f"({100*delta/exact['size']:+.6f}%).")
    print(f"Observed encode change: {gate['timings'][0]-exact['timings'][0]:+.3f} seconds "
          f"({100*(gate['timings'][0]/exact['timings'][0]-1):+.2f}%).\n")
    print("| Opcode | Name | Exact records | Gate records | Exact bytes | Gate bytes |")
    print("|---:|---|---:|---:|---:|---:|")
    profiles = [{r["code"]: r for r in d["profile"]["rows"]} for d in data]
    for code in sorted(set(profiles[0]) | set(profiles[1])):
        rows = [p.get(code, {}) for p in profiles]
        name = next(r["name"] for r in rows if r)
        print(f"| {code} | {name} | {rows[0].get('count',0):,} | {rows[1].get('count',0):,} | "
              f"{rows[0].get('bytes',0):,} | {rows[1].get('bytes',0):,} |")
    print()
print("Both configurations use L8 with all other preset competitors retained, ten workers and native ARM NEON.")
print("The gate skips 108 only for <=16 occupied source blocks and retains ordinary/adaptive fallback.")
print("Larger groups use exact joint selection. Timings are single-pass wall observations.")
