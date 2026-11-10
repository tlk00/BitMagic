#!/usr/bin/env python3
"""Stream bminv RS benchmark CSVs; no third-party packages required.

Reports structural-label recall, target coverage, best structural-gain recovery,
feature threshold curves, and measured timing speedups. Does not equate the
structural score with serialized bytes. Run on completed benchmark prefixes.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import time


def divide(a, b):
    return a / b if b else None


def fmt(value):
    return "n/a" if value is None else f"{value:.6g}"


def integer_bin(numerator, denominator):
    return min(100, (numerator * 100 + denominator - 1) // denominator) if denominator else 0


def analyze(prefix):
    started = time.perf_counter()
    prefix = str(prefix)
    report = Path(prefix + ".report.txt").read_text()
    if not report.rstrip().endswith("Status: complete"):
        raise ValueError(f"Incomplete benchmark: {prefix}")
    with open(prefix + ".timings.csv", newline="") as source:
        timing = list(csv.DictReader(source))
    with open(prefix + ".vectors.csv", newline="") as source:
        vectors = list(csv.DictReader(source))
    if not timing or len(vectors) < 2:
        raise ValueError("Missing timings or vector manifest")
    ids = {int(v['id']) for v in vectors}
    thresholds = sorted({float(row['threshold']) for row in timing if row['mode'] == 'rs'})
    accepted = [0] * len(thresholds)
    true_positive = [0] * len(thresholds)
    checksum_cost = [0] * len(thresholds)
    checksum_count = [0] * len(thresholds)
    best_retained = [{} for _ in thresholds]
    best_all = {}
    # Binned feature curves include negatives; no whole pair table is retained.
    feature_names = ['population_distance', 'block_distance', 'rs_distance', 'shape_distance', 'size_mismatch']
    histograms = {name: [[0, 0] for _ in range(101)] for name in feature_names}
    best_by_bin = {}  # exact rational RS distance thresholds in increments of .01
    pairs = qualified = total_cost = total_count = 0
    last = started
    with open(prefix + ".pairs.csv", newline="") as source:
        reader = csv.reader(source)
        header = next(reader)
        columns = {name: header.index(name) for name in (
            'a_id', 'b_id', 'target_id', 'a_bytes', 'b_bytes', 'a_count', 'b_count',
            'a_cost', 'b_cost', 'xor_cost', 'xor_count', 'qualified', 'block_l1',
            'region_l1', 'shape_distance')}
        for row in reader:
            a_id, b_id = int(row[columns['a_id']]), int(row[columns['b_id']])
            if a_id not in ids or b_id not in ids or a_id >= b_id:
                raise ValueError("Invalid pair IDs")
            ac, bc = int(row[columns['a_count']]), int(row[columns['b_count']])
            ab, bb = int(row[columns['a_bytes']]), int(row[columns['b_bytes']])
            cost, count = int(row[columns['xor_cost']]), int(row[columns['xor_count']])
            block, region = int(row[columns['block_l1']]), int(row[columns['region_l1']])
            q = int(row[columns['qualified']])
            target = int(row[columns['target_id']])
            gain = max(int(row[columns['a_cost']]), int(row[columns['b_cost']])) - cost
            if not 0 <= block <= region <= count or q not in (0, 1):
                raise ValueError("Invalid RS bound or pair label")
            denominator = ac + bc
            rs_bin = integer_bin(region, denominator)
            shape = float(row[columns['shape_distance']])
            bins = [integer_bin(abs(ac-bc), denominator), integer_bin(block, denominator), rs_bin,
                    min(100, max(0, math.ceil(shape*100))) if shape >= 0 else None,
                    integer_bin(abs(ab-bb), max(ab, bb))]
            for name, bucket in zip(feature_names, bins):
                if bucket is not None:
                    histograms[name][bucket][0] += 1
                    histograms[name][bucket][1] += q
            pairs += 1; qualified += q; total_cost += cost; total_count += count
            if q:
                best_all[target] = max(best_all.get(target, 0), gain)
                bins_for_target = best_by_bin.setdefault(target, [0] * 101)
                bins_for_target[rs_bin] = max(bins_for_target[rs_bin], gain)
            for k, threshold in enumerate(thresholds):
                # Match the C++ screening boundary, including equality.
                if region <= threshold * denominator:
                    accepted[k] += 1; true_positive[k] += q
                    checksum_cost[k] += cost; checksum_count[k] += count
                    if q:
                        best_retained[k][target] = max(best_retained[k].get(target, 0), gain)
            if pairs % 1000000 == 0 and time.perf_counter() - last >= 10:
                print(f"{prefix}: analyzed {pairs:,} pairs", flush=True)
                last = time.perf_counter()
    expected_pairs = len(vectors) * (len(vectors)-1)//2
    if pairs != expected_pairs:
        raise ValueError(f"Expected {expected_pairs} pairs, found {pairs}")
    # Cross-check timed runs with the independently collected exact pair table.
    for row in timing:
        if row['mode'] == 'baseline':
            expected = (pairs, qualified, total_cost, total_count)
        else:
            k = thresholds.index(float(row['threshold']))
            expected = (accepted[k], true_positive[k], checksum_cost[k], checksum_count[k])
        actual = tuple(int(row[key]) for key in ('evaluated', 'qualified', 'checksum_cost', 'checksum_count'))
        if actual != expected:
            raise ValueError(f"Timing/label mismatch: {row['mode']} {row['threshold']}")
    baseline_times = [float(row['seconds']) for row in timing if row['mode'] == 'baseline']
    baseline_s = statistics.median(baseline_times)
    load_s = float(timing[0]['load_seconds']); profile_s = float(timing[0]['profile_seconds'])
    gain_total = sum(best_all.values())
    summary_rows = []
    for k, threshold in enumerate(thresholds):
        runs = [float(row['seconds']) for row in timing if row['mode'] == 'rs' and float(row['threshold']) == threshold]
        elapsed = statistics.median(runs)
        row = dict(threshold=threshold, pairs=pairs, accepted=accepted[k], rejected=pairs-accepted[k],
                   true_positive=true_positive[k], false_negative=qualified-true_positive[k],
                   false_positive=accepted[k]-true_positive[k],
                   true_negative=pairs-qualified-accepted[k]+true_positive[k],
                   precision=divide(true_positive[k], accepted[k]), recall=divide(true_positive[k], qualified),
                   f1=divide(2*true_positive[k], accepted[k]+qualified),
                   evaluation_reduction=divide(pairs-accepted[k], pairs),
                   recovered_targets=len(best_retained[k]), target_recall=divide(len(best_retained[k]), len(best_all)),
                   best_structural_gain_recovery=divide(sum(best_retained[k].values()), gain_total),
                   baseline_seconds=baseline_s, screened_seconds=elapsed,
                   screened_min_seconds=min(runs), screened_max_seconds=max(runs),
                   search_speedup=divide(baseline_s, elapsed),
                   speedup_including_rs_build=divide(baseline_s, profile_s+elapsed),
                   speedup_including_load_and_rs_build=divide(load_s+baseline_s, load_s+profile_s+elapsed))
        summary_rows.append(row)
    with open(prefix + '.analysis.csv', 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(summary_rows[0]))
        writer.writeheader(); writer.writerows(summary_rows)
    target_counts = [0]*101; target_gains = [0]*101
    for bins in best_by_bin.values():
        best = 0
        for b in range(101):
            best = max(best, bins[b]); target_counts[b] += best > 0; target_gains[b] += best
    with open(prefix + '.feature-curves.csv', 'w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['feature', 'threshold', 'accepted', 'true_positive', 'false_negative', 'precision',
                         'recall', 'evaluation_reduction', 'target_recall', 'best_structural_gain_recovery'])
        for name, histogram in histograms.items():
            feature_pairs = sum(x[0] for x in histogram)
            feature_qualified = sum(x[1] for x in histogram)
            kept = positives = 0
            for b, (num, pos) in enumerate(histogram):
                kept += num; positives += pos
                writer.writerow([name, b/100, kept, positives, feature_qualified-positives,
                                 divide(positives, kept), divide(positives, feature_qualified),
                                 divide(feature_pairs-kept, feature_pairs),
                                 divide(target_counts[b], len(best_all)) if name == 'rs_distance' else '',
                                 divide(target_gains[b], gain_total) if name == 'rs_distance' else ''])
    text = [f"RS analysis: {prefix}", f"Vectors={len(vectors)}; pairs={pairs}; qualifying={qualified}; targets={len(best_all)}",
            f"Baseline search median={baseline_s:.6f}s; range={min(baseline_times):.6f}..{max(baseline_times):.6f}s",
            f"Load={load_s:.6f}s; RS build/extraction={profile_s:.6f}s",
            "All timing counts and checksums match exact pair labels.",
            "Structural gains are heuristic units, not serialized-byte savings.",
            "Including-load/build speedups are derived from separately measured stages; detailed CSV I/O is excluded.",
            "Threshold results are descriptive for these sampled sets, not held-out validation."]
    for row in summary_rows:
        text.append(f"cutoff={row['threshold']:g}: rejected={row['evaluation_reduction']:.2%}, "
                    f"pair recall={fmt(row['recall'])}, target recall={fmt(row['target_recall'])}, "
                    f"best gain recovery={fmt(row['best_structural_gain_recovery'])}, "
                    f"search speedup={fmt(row['search_speedup'])}x, "
                    f"including RS build={fmt(row['speedup_including_rs_build'])}x, "
                    f"including load={fmt(row['speedup_including_load_and_rs_build'])}x")
    text.append(f"Python analysis seconds={time.perf_counter()-started:.6f}")
    Path(prefix + '.analysis.txt').write_text('\n'.join(text) + '\n')
    print('\n'.join(text), flush=True)
    return dict(prefix=prefix, vectors=len(vectors), pairs=pairs, qualifying=qualified, results=summary_rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', nargs='+', help='Benchmark output prefixes (without .pairs.csv)')
    parser.add_argument('--summary-json', help='Optional combined analysis JSON file')
    args = parser.parse_args()
    results = [analyze(prefix) for prefix in args.prefix]
    if args.summary_json:
        Path(args.summary_json).write_text(json.dumps(results, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
