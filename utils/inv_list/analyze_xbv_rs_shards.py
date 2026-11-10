#!/usr/bin/env python3
"""Compare balanced RS and serialized-size shards using completed RS pair labels.

This is an offline candidate-coverage experiment. It does not run XOR scoring.
"""

import argparse
from array import array
import csv
import math
import random
import statistics
import sys
import time
from pathlib import Path


def options():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("prefix", help="completed -rs-bench output prefix")
    ap.add_argument("--sizes", default="500,1000,2000",
                    help="maximum primary vectors per shard (default: 500,1000,2000)")
    ap.add_argument("--overlap", default="0,0.1",
                    help="fraction of vectors also placed in a second shard (default: 0,0.1)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", help="output prefix (default: PREFIX.shards)")
    args = ap.parse_args()
    try:
        args.sizes = sorted(set(int(x) for x in args.sizes.split(",")))
        args.overlap = sorted(set(float(x) for x in args.overlap.split(",")))
    except ValueError as exc:
        ap.error(str(exc))
    if not args.sizes or any(x < 2 for x in args.sizes):
        ap.error("--sizes must contain integers >= 2")
    if not args.overlap or any(not math.isfinite(x) or x < 0 or x >= 1 for x in args.overlap):
        ap.error("--overlap fractions must be in [0,1)")
    return args


def read_features(prefix):
    report = Path(prefix + ".report.txt").read_text()
    if "Status: complete" not in report:
        raise ValueError("RS benchmark report is not complete")
    with open(prefix + ".vectors.csv", newline="") as src:
        vectors = list(csv.DictReader(src))
    ids = [v["id"] for v in vectors]
    if len(ids) != len(set(ids)) or len(ids) < 2:
        raise ValueError("Need at least two distinct vector IDs")
    pos = {v: i for i, v in enumerate(ids)}
    regions = [array("H") for _ in ids]
    with open(prefix + ".profiles.csv", newline="") as src:
        reader = csv.DictReader(src)
        for row in reader:
            i = pos[row["id"]]
            if int(row["block"]) != len(regions[i]) // 3:
                raise ValueError("Missing or out-of-order RS block")
            regions[i].extend((int(row["first_count"]), int(row["second_count"]),
                               int(row["third_count"])))
    width = max(map(len, regions))
    if width % 3:
        raise ValueError("Incomplete three-region profile")
    for i, row in enumerate(regions):
        if sum(row) != int(vectors[i]["popcount"]):
            raise ValueError("RS population does not match vector popcount")
        row.extend([0] * (width - len(row)))

    def aggregate(blocks_per_group):
        span = 3 * blocks_per_group
        return [array("I", (sum(row[j:j + span]) for j in range(0, width, span)))
                for row in regions]

    features = {"superblock": aggregate(256), "blocks16": aggregate(16),
                "blocks8": aggregate(8), "regions3": regions}
    sizes = [int(v["serialized_bytes"]) for v in vectors]
    memory = [int(v["memory_bytes"]) for v in vectors]
    return ids, sizes, memory, features


def distance(a, b, a_total, b_total):
    denom = a_total + b_total
    if not denom:
        return 0.0
    total = 0
    for x, y in zip(a, b):
        total += abs(x - y)
    return total / denom


def capacities(n, k):
    return [n // k + (g < n % k) for g in range(k)]


def size_clusters(sizes, cap):
    n = len(sizes)
    k = math.ceil(n / cap)
    limits = capacities(n, k)
    labels = [-1] * n
    sorted_ids = sorted(range(n), key=lambda i: (sizes[i], i))
    cursor = 0
    centers = []
    for group, count in enumerate(limits):
        members = sorted_ids[cursor:cursor + count]
        for i in members:
            labels[i] = group
        centers.append(statistics.median(math.log2(max(1, sizes[i])) for i in members))
        cursor += count
    distances = [[abs(math.log2(max(1, s)) - center) for center in centers]
                 for s in sizes]
    return labels, distances


def choose_seeds(rows, totals, k, rng):
    n = len(rows)
    seeds = [rng.randrange(n)]
    while len(seeds) < k:
        available = [i for i in rng.sample(range(n), min(n, 64)) if i not in seeds]
        if not available:
            available = [i for i in range(n) if i not in seeds]
        next_seed = max(available,
                        key=lambda i: min(distance(rows[i], rows[s], totals[i], totals[s])
                                          for s in seeds))
        seeds.append(next_seed)
    return [rows[i] for i in seeds]


def assign(rows, totals, centers, limits):
    center_totals = [sum(center) for center in centers]
    edges = []
    distances = []
    for i, row in enumerate(rows):
        ds = [distance(row, center, totals[i], center_totals[g])
              for g, center in enumerate(centers)]
        distances.append(ds)
        edges.extend((d, i, g) for g, d in enumerate(ds))
    edges.sort()
    labels = [-1] * len(rows)
    remaining = limits[:]
    for _, i, g in edges:
        if labels[i] < 0 and remaining[g]:
            labels[i] = g
            remaining[g] -= 1
    if any(g < 0 for g in labels) or any(remaining):
        raise ValueError("Balanced assignment failed")
    return labels, distances


def rs_clusters(rows, cap, seed):
    n = len(rows)
    k = math.ceil(n / cap)
    limits = capacities(n, k)
    totals = [sum(row) for row in rows]
    rng = random.Random(seed)
    centers = choose_seeds(rows, totals, k, rng)
    labels, _ = assign(rows, totals, centers, limits)
    # One centroid update followed by a fresh capacity-constrained assignment.
    sums = [[0] * len(rows[0]) for _ in range(k)]
    for i, row in enumerate(rows):
        dst = sums[labels[i]]
        for j, value in enumerate(row):
            dst[j] += value
    centers = [[(x + limits[g] // 2) // limits[g] for x in sums[g]]
               for g in range(k)]
    return assign(rows, totals, centers, limits)


def membership(labels, distances, fraction):
    masks = [1 << group for group in labels]
    extra = round(len(labels) * fraction)
    if extra:
        boundary = []
        for i, group in enumerate(labels):
            alternate = min((g for g in range(len(distances[i])) if g != group),
                            key=lambda g: (distances[i][g], g))
            margin = distances[i][alternate] - distances[i][group]
            boundary.append((margin, i, alternate))
        boundary.sort()
        for _, i, alternate in boundary[:extra]:
            masks[i] |= 1 << alternate
    return masks


def run():
    args = options()
    prefix = args.prefix
    out = args.out or prefix + ".shards"
    began = time.perf_counter()
    ids, sizes, memory, features = read_features(prefix)
    n = len(ids)
    print(f"Loaded {n} vector profiles; full labels: {n * (n - 1) // 2:,} pairs", flush=True)

    strategies = []
    for cap in args.sizes:
        for name in ("size", "superblock", "blocks16", "blocks8", "regions3"):
            if cap >= n:
                continue
            print(f"Clustering {name}, cap={cap}...", flush=True)
            start = time.perf_counter()
            if name == "size":
                labels, distances = size_clusters(sizes, cap)
            else:
                labels, distances = rs_clusters(features[name], cap,
                                                args.seed + cap * 13 +
                                                ("superblock", "blocks16", "blocks8", "regions3").index(name))
            cluster_seconds = time.perf_counter() - start
            for fraction in args.overlap:
                masks = membership(labels, distances, fraction)
                counts = [0] * (max(labels) + 1)
                group_memory = [0] * len(counts)
                for i, mask in enumerate(masks):
                    for g in range(len(counts)):
                        if mask & (1 << g):
                            counts[g] += 1
                            group_memory[g] += memory[i]
                strategies.append({"name": name, "cap": cap, "overlap": fraction,
                                   "groups": len(counts), "min_members": min(counts),
                                   "max_members": max(counts), "max_memory": max(group_memory),
                                   "seconds": cluster_seconds,
                                   "masks": masks, "pairs": 0, "positive": 0, "best": {}})
            print(f"  {cluster_seconds:.2f}s; {len(labels)} primary assignments", flush=True)
    if not strategies:
        raise ValueError("Every requested shard cap is at least the number of vectors")

    id_pos = {v: i for i, v in enumerate(ids)}
    strategy_masks = [s["masks"] for s in strategies]
    all_best = {}
    all_positive = 0
    pair_rows = 0
    scan_start = time.perf_counter()
    with open(prefix + ".pairs.csv", newline="") as src:
        reader = csv.reader(src)
        header = next(reader)
        needed = ("a_id", "b_id", "target_id", "a_cost", "b_cost", "xor_cost", "qualified")
        columns = {name: header.index(name) for name in needed}
        for row in reader:
            pair_rows += 1
            a = id_pos[row[columns["a_id"]]]
            b = id_pos[row[columns["b_id"]]]
            qualified = row[columns["qualified"]] == "1"
            if qualified:
                all_positive += 1
                target = row[columns["target_id"]]
                gain = max(int(row[columns["a_cost"]]),
                           int(row[columns["b_cost"]])) - int(row[columns["xor_cost"]])
                all_best[target] = max(gain, all_best.get(target, 0))
            for s, masks in zip(strategies, strategy_masks):
                if masks[a] & masks[b]:
                    s["pairs"] += 1
                    if qualified:
                        s["positive"] += 1
                        best = s["best"]
                        best[target] = max(gain, best.get(target, 0))
            if pair_rows % 4_000_000 == 0:
                print(f"Labeled pairs inspected: {pair_rows:,}", flush=True)
    if pair_rows != n * (n - 1) // 2:
        raise ValueError(f"Incomplete pair table: {pair_rows} rows")
    scan_seconds = time.perf_counter() - scan_start
    total_gain = sum(all_best.values())
    target_count = len(all_best)
    path = Path(out + ".csv")
    with path.open("w", newline="") as dst:
        writer = csv.writer(dst)
        writer.writerow(("method", "primary_cap", "overlap_fraction", "groups",
                         "smallest_group", "largest_group", "largest_group_memory_bytes",
                         "pairs_evaluated",
                         "pair_fraction", "qualifying_pairs_retained", "pair_recall",
                         "targets_recovered", "target_recall", "best_gain_recovery",
                         "clustering_seconds"))
        for s in strategies:
            writer.writerow((s["name"], s["cap"], s["overlap"], s["groups"],
                             s["min_members"], s["max_members"], s["max_memory"], s["pairs"],
                             s["pairs"] / pair_rows, s["positive"],
                             s["positive"] / all_positive if all_positive else "",
                             len(s["best"]), len(s["best"]) / target_count if target_count else "",
                             sum(s["best"].values()) / total_gain if total_gain else "",
                             s["seconds"]))
    summary = (f"Vectors={n}; labeled pairs={pair_rows}; qualifying={all_positive}; "
               f"candidate targets={target_count}\n"
               f"Profile and clustering seconds={scan_start - began:.2f}; "
               f"pair-label scan seconds={scan_seconds:.2f}; total seconds={time.perf_counter() - began:.2f}\n"
               "All methods have identical primary group sizes at each cap. "
               "Overlap adds one membership to the indicated fraction of vectors.\n"
               "Pair budgets and coverage are exact for this sample; clustering is unsupervised. "
               "Gains are structural heuristic units, not serialized bytes.\n")
    Path(out + ".report.txt").write_text(summary)
    print(summary, end="")
    print(f"Results: {path}")


if __name__ == "__main__":
    try:
        run()
    except (OSError, ValueError, KeyError, IndexError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
