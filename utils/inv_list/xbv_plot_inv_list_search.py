#!/usr/bin/env python3
"""Render separate static dot plots of inv_list XOR search targets and arguments.

The target plot has one dot per target ID, sized by its best candidate score.
The argument plot has one dot per CSV pair: repeated argument IDs therefore
represent multiple candidate pairs.
X is the zero-based vector ID; Y is its standalone serialized size in bytes.
Argument colors indicate the sign of the selected gain column, not whether
the builder ultimately selected that argument.
"""

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("pixel dimensions must be positive")
    return number


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="inv_list qualifying-pair CSV")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--gain-column", help="signed column for argument colors; default: gain_bytes if present, otherwise gain_cost")
    parser.add_argument("--target-score-column", help="projected gain for target dot size; default: gain_cost if present, otherwise the argument gain column")
    parser.add_argument("--width", type=positive_int, default=8000, help="PNG width in pixels (default: 8000)")
    parser.add_argument("--height", type=positive_int, default=4000, help="PNG height in pixels (default: 4000)")
    parser.add_argument("--log-bytes", action="store_true", help="use a log scale for serialized sizes")
    return parser.parse_args()


def load_pairs(path, requested_gain, requested_target_score):
    target_sizes, target_scores = {}, {}
    argument_ids, argument_sizes, gains = [], [], []
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        names = set(reader.fieldnames or ())
        required = {"lead_id", "target_id", "lead_bytes", "target_bytes"}
        missing = required - names
        if missing:
            raise ValueError(f"missing CSV columns: {', '.join(sorted(missing))}")
        gain_column = requested_gain or next(
            (name for name in ("gain_bytes", "gain_cost") if name in names), None
        )
        if not gain_column or gain_column not in names:
            raise ValueError("no gain column found; pass --gain-column COLUMN")
        score_column = requested_target_score or ("gain_cost" if "gain_cost" in names else gain_column)
        if score_column not in names:
            raise ValueError(f"target score column not found: {score_column}")
        for line_no, row in enumerate(reader, start=2):
            try:
                argument_id = int(row["lead_id"])
                target_id = int(row["target_id"])
                argument_bytes = int(row["lead_bytes"])
                target_bytes = int(row["target_bytes"])
                gain = float(row[gain_column])
                score = float(row[score_column])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid CSV row at line {line_no}: {exc}") from exc
            if min(argument_id, target_id) < 0 or min(argument_bytes, target_bytes) <= 0:
                raise ValueError(f"invalid ID or byte size at line {line_no}")
            if not math.isfinite(gain) or not math.isfinite(score):
                raise ValueError(f"non-finite gain or target score at line {line_no}")
            previous = target_sizes.setdefault(target_id, target_bytes)
            if previous != target_bytes:
                raise ValueError(f"target {target_id} has conflicting byte sizes")
            target_scores[target_id] = max(score, target_scores.get(target_id, -math.inf))
            argument_ids.append(argument_id)
            argument_sizes.append(argument_bytes)
            gains.append(gain)
    if not gains:
        raise ValueError("CSV contains no pair rows")
    return target_sizes, target_scores, np.asarray(argument_ids), np.asarray(argument_sizes), np.asarray(gains), gain_column, score_column


def marker_areas(scores):
    """Map nonnegative scores to marker areas (points squared) using log1p."""
    transformed = np.log1p(np.maximum(scores, 0))
    low, high = transformed.min(), transformed.max()
    if high == low:
        return np.full(scores.shape, 32.0)
    return 5.0 + 75.0 * (transformed - low) / (high - low)


def make_figure(width, height):
    dpi = 100
    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.12, top=0.86)
    ax.grid(True, color="#e6e6e6", linewidth=0.5, zorder=0)
    ax.set_xlabel("Vector ID (zero-based input record)")
    ax.set_ylabel("Standalone serialized bytes")
    ax.ticklabel_format(axis="x", style="plain", useOffset=False)
    return fig, ax, dpi


def finish(fig, ax, path, dpi, log_bytes):
    if log_bytes:
        ax.set_yscale("log")
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)


def main():
    args = parse_args()
    target_sizes, target_scores, argument_ids, argument_sizes, gains, gain_column, score_column = load_pairs(
        args.csv, args.gain_column, args.target_score_column
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    status_file = Path(str(args.csv) + ".report.txt")
    status = ""
    if status_file.exists():
        report = status_file.read_text(encoding="utf-8", errors="replace")
        if "Status: incomplete" in report:
            status = " · partial search"

    ids = np.fromiter(target_sizes.keys(), dtype=np.int64)
    sizes = np.fromiter(target_sizes.values(), dtype=np.int64)
    scores = np.fromiter((target_scores[id_] for id_ in ids), dtype=np.float64)
    areas = marker_areas(scores)
    fig, ax, dpi = make_figure(args.width, args.height)
    ax.scatter(ids, sizes, s=areas, c="#2166d1", alpha=0.65, edgecolors="none")
    ax.set_title("Compression targets", loc="left", pad=25)
    ax.text(0, 1.015, f"{len(ids):,} unique vectors{status} · dot area scales with log(1 + best {score_column})",
            transform=ax.transAxes, fontsize=9, color="#555555")
    for score in sorted(set(np.quantile(scores, [0, 0.5, 1]))):
        area = marker_areas(np.asarray([score, scores.min(), scores.max()]))[0]
        ax.scatter([], [], s=area, c="#2166d1", alpha=0.65,
                   label=f"{score:,.0f}")
    ax.legend(title=f"Best {score_column}", loc="upper left", frameon=True,
              labelspacing=1.5, borderpad=1)
    target_path = args.output_dir / "inv_list_targets.png"
    finish(fig, ax, target_path, dpi, args.log_bytes)

    fig, ax, dpi = make_figure(args.width, args.height)
    classes = (
        (gains < 0, "Negative gain", "#c93c3c"),
        (gains == 0, "Zero gain", "#888888"),
        (gains > 0, "Positive gain", "#168a55"),
    )
    for mask, label, color in classes:
        if mask.any():
            ax.scatter(argument_ids[mask], argument_sizes[mask], s=4, c=color,
                       alpha=0.45, edgecolors="none", label=f"{label}: {mask.sum():,} pairs")
    ax.set_title("Compression arguments", loc="left", pad=25)
    ax.text(0, 1.015, f"{len(gains):,} candidate pairs{status} · color = sign of {gain_column}",
            transform=ax.transAxes, fontsize=9, color="#555555")
    ax.legend(loc="upper right", frameon=True, markerscale=4)
    argument_path = args.output_dir / "inv_list_arguments.png"
    finish(fig, ax, argument_path, dpi, args.log_bytes)

    print(f"Gain column: {gain_column}")
    print(f"Target score column: {score_column}; per-target maximum, log-scaled marker area")
    print(f"Pairs: {len(gains):,}; unique targets: {len(ids):,}; unique arguments: {len(set(argument_ids)):,}")
    print(f"Positive/zero/negative pairs: {(gains > 0).sum():,}/{(gains == 0).sum():,}/{(gains < 0).sum():,}")
    if gain_column == "gain_cost":
        print("Note: gain_cost is a structural heuristic, not measured byte savings.")
    if status:
        print("Note: companion search report marks this CSV as incomplete.")
    print(target_path)
    print(argument_path)


if __name__ == "__main__":
    main()
