#!/usr/bin/env python3
"""Draw a search-derived XOR dependency forest from an inv_list pair CSV.

Each pair is oriented as the XOR builder does: the vector with the smaller
(standalone serialized bytes, ID) is the argument. For each possible target,
the script retains only its highest positive projected structural gain. This
one-parent choice makes a forest. It is a *search projection*, not a plot of
the references actually selected by -xor-build.
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
import shutil
import subprocess


def positive_int(value):
    number = int(value)
    if number <= 0:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return number


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="inv_list qualifying-pair CSV")
    parser.add_argument("--output-dir", type=Path, default=Path("."))
    parser.add_argument("--top-components", type=positive_int, default=6,
                        help="number of largest trees to draw, one PNG per tree (default: 6)")
    parser.add_argument("--root", type=int, help="instead draw the tree containing this vector ID")
    parser.add_argument("--all-components", action="store_true", help="draw every tree; may create many files")
    parser.add_argument("--min-gain", type=positive_int, default=1,
                        help="minimum projected structural gain after builder orientation (default: 1)")
    parser.add_argument("--dpi", type=positive_int, default=180,
                        help="Graphviz PNG resolution (default: 180)")
    args = parser.parse_args()
    if args.root is not None and args.root < 0:
        parser.error("--root must be nonnegative")
    if args.root is not None and args.all_components:
        parser.error("choose --root or --all-components, not both")
    return args


def read_forest(path, min_gain):
    best = {}  # target -> (projected gain, argument ID)
    sizes = {}
    rows = 0
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"lead_id", "target_id", "lead_cost", "target_cost",
                    "xor_cost", "lead_bytes", "target_bytes"}
        missing = required - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"missing CSV columns: {', '.join(sorted(missing))}")
        for line_no, row in enumerate(reader, start=2):
            try:
                a, b = int(row["lead_id"]), int(row["target_id"])
                ac, bc, xc = int(row["lead_cost"]), int(row["target_cost"]), int(row["xor_cost"])
                ab, bb = int(row["lead_bytes"]), int(row["target_bytes"])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid CSV row at line {line_no}: {exc}") from exc
            if min(a, b) < 0 or a == b or min(ab, bb) <= 0:
                raise ValueError(f"invalid ID or byte size at line {line_no}")
            for vector_id, size in ((a, ab), (b, bb)):
                previous = sizes.setdefault(vector_id, size)
                if previous != size:
                    raise ValueError(f"vector {vector_id} has conflicting byte sizes")
            rows += 1
            if (ab, a) < (bb, b):
                argument, target, projected_gain = a, b, bc - xc
            else:
                argument, target, projected_gain = b, a, ac - xc
            if projected_gain < min_gain:
                continue
            previous = best.get(target)
            if previous is None or projected_gain > previous[0] or (
                projected_gain == previous[0] and argument < previous[1]
            ):
                best[target] = (projected_gain, argument)
    if not rows:
        raise ValueError("CSV contains no pair rows")
    if not best:
        raise ValueError("no positive projected dependencies pass --min-gain")
    children = defaultdict(list)
    for target, (_, argument) in best.items():
        if (sizes[argument], argument) >= (sizes[target], target):
            raise ValueError("dependency order is not acyclic")
        children[argument].append(target)
    for followers in children.values():
        followers.sort()
    return rows, sizes, best, children


def components(best, children):
    roots = sorted(set(children) - set(best))
    result = []
    node_to_root = {}
    for root in roots:
        nodes = []
        max_depth = 0
        stack = [(root, 0)]
        while stack:
            node, depth = stack.pop()
            if node in node_to_root:
                raise ValueError("cycle or shared target in forest")
            node_to_root[node] = root
            nodes.append(node)
            max_depth = max(max_depth, depth)
            stack.extend((child, depth + 1) for child in reversed(children[node]))
        result.append((root, nodes, max_depth))
    if len(node_to_root) != len(set(children) | set(best)):
        raise ValueError("forest contains an unreachable node or a cycle")
    result.sort(key=lambda item: (-len(item[1]), -item[2], item[0]))
    return result, node_to_root


def write_summary(path, all_components, best):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("rank", "root_id", "nodes", "edges", "max_depth", "sum_projected_gain_cost"))
        for rank, (root, nodes, depth) in enumerate(all_components, start=1):
            writer.writerow((rank, root, len(nodes), len(nodes) - 1, depth,
                             sum(best[node][0] for node in nodes if node in best)))


def dot_text(selected, sizes, best, children, rows, total_components, status):
    selected_nodes = {node for _, nodes, _ in selected for node in nodes}
    max_gain = max(best[node][0] for node in selected_nodes if node in best)
    lines = [
        "digraph xor_forest {",
        '  graph [rankdir=TB, bgcolor="white", fontname="Helvetica", fontsize=18,',
        '         labeljust=l, labelloc=t, nodesep=0.15, ranksep=0.55, splines=polyline];',
        '  node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=9,',
        '        color="#596477", penwidth=0.6, margin="0.07,0.04"];',
        '  edge [color="#598d71", arrowsize=0.5, fontname="Helvetica", fontsize=8];',
        f'  label="Search-derived XOR dependency forest\\n{rows:,} input pairs; '
        f'{len(selected)} of {total_components:,} trees shown{status}";',
    ]
    for index, (root, nodes, depth) in enumerate(selected, start=1):
        lines.append(f"  subgraph cluster_{index} {{")
        lines.append('    color="#d7dce5"; style=rounded; margin=18;')
        lines.append(f'    label="Tree {index}: root {root}; {len(nodes)} vectors; depth {depth}";')
        for node in nodes:
            if node == root:
                fill = "#58ad76"  # standalone argument/root
            elif children[node]:
                fill = "#a18acb"  # both target and argument
            else:
                fill = "#6c9ddd"  # target only
            lines.append(f'    n{node} [label="{node}\\n{sizes[node]} B", fillcolor="{fill}"];')
        for target in nodes:
            if target not in best:
                continue
            gain, argument = best[target]
            width = 1.0 + 2.0 * math.log1p(gain) / math.log1p(max_gain)
            lines.append(f"    n{argument} -> n{target} [penwidth={width:.2f}];")
        lines.append("  }")
    lines.append("}")
    return "\n".join(lines) + "\n"


def main():
    args = parse_args()
    dot = shutil.which("dot")
    if not dot:
        raise SystemExit("Graphviz 'dot' is required; install Graphviz or add dot to PATH")
    rows, sizes, best, children = read_forest(args.csv, args.min_gain)
    all_components, node_to_root = components(best, children)
    if args.root is not None:
        root = node_to_root.get(args.root)
        if root is None:
            raise ValueError(f"vector {args.root} has no retained dependency in this forest")
        selected = [next(component for component in all_components if component[0] == root)]
    elif args.all_components:
        selected = all_components
    else:
        selected = all_components[:args.top_components]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = args.output_dir / "inv_list_forest_components.csv"
    write_summary(summary_path, all_components, best)
    report_path = Path(str(args.csv) + ".report.txt")
    status = ""
    if report_path.exists() and "Status: incomplete" in report_path.read_text(encoding="utf-8", errors="replace"):
        status = " (partial search)"
    rank_by_root = {root: rank for rank, (root, _, _) in enumerate(all_components, start=1)}
    png_paths = []
    for component in selected:
        root = component[0]
        basename = f"inv_list_forest_rank_{rank_by_root[root]:04d}_root_{root}"
        dot_path = args.output_dir / f"{basename}.dot"
        png_path = args.output_dir / f"{basename}.png"
        dot_path.write_text(dot_text([component], sizes, best, children, rows,
                                     len(all_components), status), encoding="utf-8")
        subprocess.run([dot, f"-Gdpi={args.dpi}", "-Tpng", str(dot_path), "-o", str(png_path)], check=True)
        png_paths.append(png_path)
    print(f"Input pairs: {rows:,}; retained dependencies: {len(best):,}; trees: {len(all_components):,}")
    print(f"Shown: {len(selected)} trees, {sum(len(nodes) for _, nodes, _ in selected):,} vectors")
    print("Green roots; purple vectors that are both targets and arguments; blue target leaves.")
    print("Edges use the best positive projected structural gain per target; this is not the built store.")
    if status:
        print("Note: companion search report marks this CSV as incomplete.")
    print(summary_path)
    for path in png_paths:
        print(path)


if __name__ == "__main__":
    main()
