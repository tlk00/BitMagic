#!/usr/bin/env python3
"""Combine completed 7K/3K pair searches and quantify cross-cluster coverage."""

import argparse
from collections import Counter
import csv
from pathlib import Path
import sys


def require_complete(path):
    if "Status: complete" not in path.read_text(encoding="utf-8"):
        raise ValueError(f"Incomplete report: {path}")


def read_manifest(path):
    require_complete(Path(str(path) + ".report.txt"))
    result = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            vector = int(row["vector_id"])
            if vector in result:
                raise ValueError(f"Duplicate vector {vector} in {path}")
            result[vector] = (int(row["cluster_id"]), int(row["source_offset"]),
                              int(row["source_bytes"]), int(row["indexed_input_bytes"]))
    if not result:
        raise ValueError(f"Empty manifest: {path}")
    return result


def read_pairs(path):
    require_complete(Path(str(path) + ".report.txt"))
    result = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames
        for row in reader:
            key = (int(row["lead_id"]), int(row["target_id"]))
            if key in result:
                raise ValueError(f"Duplicate pair {key} in {path}")
            result[key] = row
    if not result:
        raise ValueError(f"Empty pair CSV: {path}")
    return columns, result


def choose2(n):
    return n * (n - 1) // 2


def structural_best(rows):
    best = {}
    for (_, target), row in rows.items():
        gain = int(row["gain_cost"])
        if gain > best.get(target, -1):
            best[target] = gain
    return best


def builder_top4(rows):
    candidates = {}
    for (lead, target), row in rows.items():
        lead_bytes = int(row["lead_bytes"])
        target_bytes = int(row["target_bytes"])
        if (lead_bytes, lead) < (target_bytes, target):
            dep, storage_target, cost = lead, target, int(row["target_cost"])
        else:
            dep, storage_target, cost = target, lead, int(row["lead_cost"])
        gain = cost - int(row["xor_cost"])
        candidates.setdefault(storage_target, {})[dep] = gain
    return {target: tuple(sorted(deps.items(), key=lambda item: (-item[1], item[0]))[:4])
            for target, deps in candidates.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--large", default="gov2.rs.cluster.regions3")
    parser.add_argument("--small", default="gov2.rs.cluster.regions3.3k")
    parser.add_argument("--pairs-out", default="gov2.rs.cluster.union.pairs.csv")
    parser.add_argument("--report-out", default="gov2.rs.cluster.cross-eval.txt")
    args = parser.parse_args()
    pairs_out = Path(args.pairs_out)
    report_out = Path(args.report_out)
    for output in (pairs_out, report_out):
        if output.exists():
            raise ValueError(f"Output already exists: {output}")

    large_manifest = read_manifest(Path(args.large + ".manifest.csv"))
    small_manifest = read_manifest(Path(args.small + ".manifest.csv"))
    if large_manifest.keys() != small_manifest.keys():
        raise ValueError("Cluster manifests select different vector sets")
    for vector, large in large_manifest.items():
        if large[1:] != small_manifest[vector][1:]:
            raise ValueError(f"Cluster manifests disagree about vector {vector}")

    large_columns, large = read_pairs(Path(args.large + ".pairs.csv"))
    small_columns, small = read_pairs(Path(args.small + ".pairs.csv"))
    if large_columns != small_columns:
        raise ValueError("Pair CSV columns differ")
    for key in large.keys() & small.keys():
        if large[key] != small[key]:
            raise ValueError(f"Pair outcome differs between searches: {key}")
    for key in large.keys() | small.keys():
        if key[0] not in large_manifest or key[1] not in large_manifest:
            raise ValueError(f"Pair outside both manifests: {key}")

    large_sizes = Counter(entry[0] for entry in large_manifest.values())
    small_sizes = Counter(entry[0] for entry in small_manifest.values())
    intersections = Counter((large_manifest[v][0], small_manifest[v][0])
                            for v in large_manifest)
    large_space = sum(choose2(n) for n in large_sizes.values())
    small_space = sum(choose2(n) for n in small_sizes.values())
    common_space = sum(choose2(n) for n in intersections.values())
    large_cross = {key for key in large
                   if small_manifest[key[0]][0] != small_manifest[key[1]][0]}
    small_cross = {key for key in small
                   if large_manifest[key[0]][0] != large_manifest[key[1]][0]}
    if len(large_cross) != len(large.keys() - small.keys()) or \
            len(small_cross) != len(small.keys() - large.keys()):
        raise ValueError("Pair-space and qualifying-pair intersection disagree")
    union = dict(large)
    union.update(small)

    large_best, small_best, union_best = (structural_best(rows)
                                          for rows in (large, small, union))
    large_top, small_top, union_top = (builder_top4(rows)
                                      for rows in (large, small, union))
    lines = [
        "Cross-cluster evaluation from completed 7K and 3K Gov2 searches",
        f"Large input={args.large}; small input={args.small}",
        f"Eligible vectors={len(large_manifest)}",
        f"Evaluated pair spaces: large={large_space}; small={small_space}; "
        f"overlap={common_space}; union={large_space + small_space - common_space}",
        f"Additional evaluated cross-cluster pairs vs large={small_space - common_space}; "
        f"vs small={large_space - common_space}",
        f"Qualifying pairs: large={len(large)}; small={len(small)}; "
        f"overlap={len(large.keys() & small.keys())}; union={len(union)}",
        f"Large-run qualifying pairs crossing small-cluster boundaries={len(large_cross)}; "
        f"targets={len({target for _, target in large_cross})}",
        f"Small-run qualifying pairs crossing large-cluster boundaries={len(small_cross)}; "
        f"targets={len({target for _, target in small_cross})}",
    ]
    for label, old_best, old_top in (("large", large_best, large_top),
                                      ("small", small_best, small_top)):
        lines.extend([
            f"Union vs {label}: new structural targets={len(union_best.keys() - old_best.keys())}; "
            f"improved existing targets={sum(union_best[t] > old_best[t] for t in old_best)}; "
            f"best structural gain sum delta={sum(union_best.values()) - sum(old_best.values())}",
            f"Union vs {label}: new builder targets={len(union_top.keys() - old_top.keys())}; "
            f"changed top-4 candidate lists={sum(union_top[t] != old_top[t] for t in old_top)}",
        ])
    lines.extend([
        f"Union pair CSV={pairs_out}; rows={len(union)}",
        "Structural gains are heuristic. Builder candidate lists use the smaller "
        "(standalone_bytes,id) as dependency; measured store savings require -xor-build.",
        "Status: complete",
    ])

    with pairs_out.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=large_columns, lineterminator="\n")
        writer.writeheader()
        for key in sorted(union, key=lambda item: (item[1], item[0])):
            writer.writerow(union[key])
    report = "\n".join(lines) + "\n"
    report_out.write_text(report, encoding="utf-8")
    print(report, end="")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
