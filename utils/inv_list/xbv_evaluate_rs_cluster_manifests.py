#!/usr/bin/env python3
"""Score C++ RS cluster manifests against a completed all-pairs RS benchmark."""

import argparse
import csv
from pathlib import Path
import sys
import time


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("benchmark_prefix", help="completed -rs-bench output prefix")
    ap.add_argument("manifests", nargs="+", help="C++ -rs-cluster-out CSV files")
    ap.add_argument("--out", help="summary CSV (default: PREFIX.cluster-replay.csv)")
    args = ap.parse_args()
    prefix = args.benchmark_prefix
    if "Status: complete" not in Path(prefix + ".report.txt").read_text():
        raise ValueError("RS benchmark is incomplete")
    with open(prefix + ".vectors.csv", newline="") as src:
        ids = {r["id"] for r in csv.DictReader(src)}
    if len(ids) < 2:
        raise ValueError("Too few benchmark vectors")
    runs = []
    for manifest in args.manifests:
        report = Path(manifest + ".report.txt").read_text()
        if "Status: complete" not in report:
            raise ValueError(f"Cluster manifest incomplete: {manifest}")
        with open(manifest, newline="") as src:
            rows = list(csv.DictReader(src))
        groups = {}
        for row in rows:
            vector_id = row["vector_id"]
            if vector_id not in ids or vector_id in groups:
                raise ValueError(f"Manifest IDs do not match benchmark: {manifest}")
            groups[vector_id] = row["cluster_id"]
        if groups.keys() != ids:
            raise ValueError(f"Manifest must cover every benchmark vector exactly once: {manifest}")
        feature = next((line.split("=", 1)[1] for line in report.splitlines()
                        if line.startswith("Feature=")), "unknown")
        seed = next((line.split("=", 1)[1] for line in report.splitlines()
                     if line.startswith("Seed=")), "unknown")
        runs.append({"manifest": manifest, "feature": feature, "seed": seed,
                     "groups": groups, "pairs": 0, "positives": 0, "best": {}})
    expected_pairs = len(ids) * (len(ids) - 1) // 2
    all_best = {}
    positive_count = 0
    pair_count = 0
    began = time.perf_counter()
    with open(prefix + ".pairs.csv", newline="") as src:
        reader = csv.reader(src)
        header = next(reader)
        columns = {name: header.index(name) for name in
                   ("a_id", "b_id", "target_id", "a_cost", "b_cost", "xor_cost", "qualified")}
        for row in reader:
            pair_count += 1
            a, b = row[columns["a_id"]], row[columns["b_id"]]
            qualified = row[columns["qualified"]] == "1"
            if qualified:
                positive_count += 1
                target = row[columns["target_id"]]
                gain = max(int(row[columns["a_cost"]]), int(row[columns["b_cost"]])) - int(row[columns["xor_cost"]])
                all_best[target] = max(gain, all_best.get(target, 0))
            for run in runs:
                if run["groups"][a] != run["groups"][b]:
                    continue
                run["pairs"] += 1
                if qualified:
                    run["positives"] += 1
                    best = run["best"]
                    best[target] = max(gain, best.get(target, 0))
            if pair_count % 4_000_000 == 0:
                print(f"Scored {pair_count:,} labeled pairs", flush=True)
    if pair_count != expected_pairs:
        raise ValueError(f"Expected {expected_pairs} pairs, got {pair_count}")
    gain_sum = sum(all_best.values())
    out = args.out or prefix + ".cluster-replay.csv"
    with open(out, "w", newline="") as dst:
        writer = csv.writer(dst)
        writer.writerow(("manifest", "feature", "seed", "vectors", "pairs_evaluated",
                         "pair_fraction", "qualifying_pairs_retained", "pair_recall",
                         "targets_recovered", "target_recall", "best_gain_recovery"))
        for run in runs:
            writer.writerow((run["manifest"], run["feature"], run["seed"], len(ids),
                             run["pairs"], run["pairs"] / pair_count,
                             run["positives"],
                             run["positives"] / positive_count if positive_count else "",
                             len(run["best"]),
                             len(run["best"]) / len(all_best) if all_best else "",
                             sum(run["best"].values()) / gain_sum if gain_sum else ""))
    print(f"Scored {len(runs)} C++ plans against {pair_count:,} pairs in {time.perf_counter()-began:.2f}s")
    print(f"Results: {out}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, IndexError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
