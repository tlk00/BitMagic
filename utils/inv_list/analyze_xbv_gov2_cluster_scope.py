#!/usr/bin/env python3
"""Compare completed clustered searches and their verified XOR stores."""

import argparse
import csv
from pathlib import Path
import re
import sys


def number(pattern, text, path):
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"Missing statistic in {path}: {pattern}")
    return int(match.group(1))


def seconds(text, path):
    match = re.search(r"elapsed seconds=([0-9.]+)", text, re.IGNORECASE)
    if not match:
        raise ValueError(f"Missing elapsed time in {path}")
    return float(match.group(1))


def read_run(prefix):
    plan_path = Path(prefix + ".manifest.csv.report.txt")
    search_path = Path(prefix + ".pairs.csv.report.txt")
    build_path = Path(prefix + ".xor.bv.report.txt")
    verify_path = Path(prefix + ".xor.verify.log")
    plan = plan_path.read_text()
    search = search_path.read_text()
    build = build_path.read_text()
    verify = verify_path.read_text()
    if any("Status: complete" not in text for text in (plan, search, build)):
        raise ValueError(f"Plan, search, or build is incomplete: {prefix}")
    if "XOR verification complete." not in verify:
        raise ValueError(f"Verification is incomplete: {prefix}")
    with open(prefix + ".pairs.csv.best.csv", newline="") as stream:
        best = {int(row["target_id"]): int(row["gain_cost"])
                for row in csv.DictReader(stream)}
    return {
        "prefix": prefix,
        "vectors": number(r"searched vectors=(\d+)", search, search_path),
        "store_vectors": number(r"XOR store build complete[^\n]*Vectors=(\d+)", build, build_path),
        "clusters": number(r"\nClusters=(\d+)", plan, plan_path),
        "peak_memory": number(r"Maximum cluster decoded memory=(\d+)", plan, plan_path),
        "pairs": number(r"internal-shard pairs=(\d+)", search, search_path),
        "qualifying": number(r"Qualifying pairs=(\d+)", search, search_path),
        "targets": number(r"candidate targets=(\d+)", search, search_path),
        "source": number(r"Input bytes consumed=(\d+)", build, build_path),
        "standalone": number(r"standalone level-6 bytes=(\d+)", build, build_path),
        "payload": number(r"Payload file bytes=(\d+)", build, build_path),
        "metadata": number(r"Reference metadata bytes=(\d+)", build, build_path),
        "offsets": number(r"Offset index bytes=(\d+)", build, build_path),
        "store": number(r"Total XOR store bytes=(\d+)", build, build_path),
        "plain": number(r"plain XOR=(\d+) vectors", build, build_path),
        "reference": number(r"reference group=(\d+) vectors", build, build_path),
        "search_seconds": seconds(search, search_path),
        "plan_seconds": float(re.search(r"Total planning seconds=([0-9.]+)", plan).group(1)),
        "build_seconds": seconds(build, build_path),
        "verify_seconds": seconds(verify, verify_path),
        "best": best,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline")
    parser.add_argument("trial")
    parser.add_argument("--pair-index", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    baseline = read_run(args.baseline)
    trial = read_run(args.trial)
    if baseline["standalone"] != trial["standalone"]:
        raise ValueError("Standalone baselines differ; stores are not comparable")
    if baseline["source"] != trial["source"]:
        raise ValueError("CM6 source byte totals differ; stores are not comparable")
    if baseline["vectors"] != trial["vectors"]:
        raise ValueError("Search vector sets differ")
    with open(args.pair_index, newline="") as stream:
        rows = csv.DictReader(stream)
        first = next(rows, None)
        if first is None:
            raise ValueError("Empty pair index")
        total_ints = int(first["count"])
        indexed_bytes = int(first["indexed_input_bytes"])
        indexed_vectors = 1
        for row in rows:
            total_ints += int(row["count"])
            indexed_vectors += 1
    if indexed_vectors != baseline["store_vectors"] or indexed_vectors != trial["store_vectors"]:
        raise ValueError("Pair index vector count does not match XOR stores")
    if indexed_bytes != baseline["source"] or total_ints == 0:
        raise ValueError("Pair index does not match XOR source or has no integers")
    b_targets = set(baseline["best"])
    t_targets = set(trial["best"])
    lines = [
        "Cluster-scope experiment: verified whole-collection XOR stores",
        f"Baseline: {args.baseline}",
        f"Trial: {args.trial}",
        "",
        "Metric                              Baseline             Trial",
    ]
    for label, key in (
        ("Clusters", "clusters"),
        ("Internal pairs", "pairs"),
        ("Qualifying pairs", "qualifying"),
        ("Structural candidate targets", "targets"),
        ("Standalone level-6 bytes", "standalone"),
        ("Input CM6 bytes", "source"),
        ("Total integers", "total_ints"),
        ("Plain XOR vectors", "plain"),
        ("Reference-group vectors", "reference"),
        ("Payload bytes", "payload"),
        ("Reference metadata bytes", "metadata"),
        ("Offset index bytes", "offsets"),
        ("Total store bytes", "store"),
        ("Peak cluster memory estimate", "peak_memory"),
    ):
        left = total_ints if key == "total_ints" else baseline[key]
        right = total_ints if key == "total_ints" else trial[key]
        lines.append(f"{label:<35} {left:>18,} {right:>18,}")
    for run in (baseline, trial):
        run["saving"] = run["standalone"] - run["store"]
        run["saving_pct"] = 100 * run["saving"] / run["standalone"]
        run["source_saving_pct"] = 100 * (run["source"] - run["store"]) / run["source"]
        run["best_gain_sum"] = sum(run["best"].values())
    lines.extend([
        f"{'Net saving bytes':<35} {baseline['saving']:>18,} {trial['saving']:>18,}",
        f"{'Net saving percent':<35} {baseline['saving_pct']:>17.4f}% {trial['saving_pct']:>17.4f}%",
        f"{'Saving vs input CM6':<35} {baseline['source_saving_pct']:>17.4f}% {trial['source_saving_pct']:>17.4f}%",
        f"{'Input CM6 bits/int':<35} {8 * baseline['source'] / total_ints:>18.6f} {8 * trial['source'] / total_ints:>18.6f}",
        f"{'Standalone bits/int':<35} {8 * baseline['standalone'] / total_ints:>18.6f} {8 * trial['standalone'] / total_ints:>18.6f}",
        f"{'XOR store bits/int':<35} {8 * baseline['store'] / total_ints:>18.6f} {8 * trial['store'] / total_ints:>18.6f}",
        f"{'Search seconds':<35} {baseline['search_seconds']:>18.2f} {trial['search_seconds']:>18.2f}",
        f"{'Planning seconds':<35} {baseline['plan_seconds']:>18.2f} {trial['plan_seconds']:>18.2f}",
        f"{'Build seconds':<35} {baseline['build_seconds']:>18.2f} {trial['build_seconds']:>18.2f}",
        f"{'Verification seconds':<35} {baseline['verify_seconds']:>18.2f} {trial['verify_seconds']:>18.2f}",
        f"{'Best structural gain sum':<35} {baseline['best_gain_sum']:>18,} {trial['best_gain_sum']:>18,}",
        "",
        f"Trial minus baseline store bytes: {trial['store'] - baseline['store']:,}",
        f"Trial minus baseline saving bytes: {trial['saving'] - baseline['saving']:,}",
        f"Best-pair target overlap: {len(b_targets & t_targets):,}",
        f"Baseline-only best-pair targets: {len(b_targets - t_targets):,}",
        f"Trial-only best-pair targets: {len(t_targets - b_targets):,}",
        f"Trial pair-budget reduction: {100 * (baseline['pairs'] - trial['pairs']) / baseline['pairs']:.2f}%",
        f"Shared targets with larger trial gain: {sum(trial['best'][i] > baseline['best'][i] for i in b_targets & t_targets):,}",
        f"Shared targets with larger baseline gain: {sum(trial['best'][i] < baseline['best'][i] for i in b_targets & t_targets):,}",
        "Structural gains are heuristic; actual store bytes include all payload and metadata files.",
    ])
    report = "\n".join(lines) + "\n"
    Path(args.out).write_text(report)
    print(report, end="")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
