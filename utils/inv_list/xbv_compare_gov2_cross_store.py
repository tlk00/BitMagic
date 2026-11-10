#!/usr/bin/env python3
"""Compare three verified whole-collection Gov2 XOR stores."""

import argparse
from pathlib import Path
import re
import sys


def extract(text, pattern, path):
    match = re.search(pattern, text)
    if not match:
        raise ValueError(f"Missing statistic in {path}: {pattern}")
    return int(match.group(1))


def read_store(prefix):
    build_path = Path(prefix + ".xor.bv.report.txt")
    verify_path = Path(prefix + ".xor.verify.log")
    build = build_path.read_text(encoding="utf-8")
    verify = verify_path.read_text(encoding="utf-8")
    if not re.search(r"^Status: complete$", build, re.MULTILINE):
        raise ValueError(f"Incomplete build: {build_path}")
    matches = re.search(r"XOR verification complete\. Vectors=(\d+); exact matches=(\d+)", verify)
    if not matches or matches.group(1) != matches.group(2):
        raise ValueError(f"Incomplete or inexact verification: {verify_path}")
    return {
        "vectors": int(matches.group(1)),
        "source": extract(build, r"Input bytes consumed=(\d+)", build_path),
        "standalone": extract(build, r"standalone level-6 bytes=(\d+)", build_path),
        "payload": extract(build, r"Payload file bytes=(\d+)", build_path),
        "metadata": extract(build, r"Reference metadata bytes=(\d+)", build_path),
        "offsets": extract(build, r"Offset index bytes=(\d+)", build_path),
        "store": extract(build, r"Total XOR store bytes=(\d+)", build_path),
        "plain": extract(build, r"plain XOR=(\d+) vectors", build_path),
        "group": extract(build, r"reference group=(\d+) vectors", build_path),
        "build_seconds": float(re.search(r"Elapsed seconds=([0-9.]+)", build).group(1)),
        "verify_seconds": float(re.search(r"elapsed seconds=([0-9.]+)", verify).group(1)),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--large", default="gov2.rs.cluster.regions3")
    parser.add_argument("--small", default="gov2.rs.cluster.regions3.3k")
    parser.add_argument("--union", default="gov2.rs.cluster.union")
    parser.add_argument("--out", default="gov2.rs.cluster.cross-store.txt")
    args = parser.parse_args()
    stores = {label: read_store(getattr(args, label))
              for label in ("large", "small", "union")}
    if len({store["vectors"] for store in stores.values()}) != 1 or \
            len({store["source"] for store in stores.values()}) != 1 or \
            len({store["standalone"] for store in stores.values()}) != 1:
        raise ValueError("Stores do not share the same vector count and baselines")
    large, small, union = (stores[label] for label in ("large", "small", "union"))
    lines = [
        "Verified whole-collection Gov2 XOR stores",
        f"Large={args.large}; small={args.small}; union={args.union}",
        f"Vectors={union['vectors']}; source CM6 bytes={union['source']}; "
        f"standalone level-6 bytes={union['standalone']}",
        "",
        f"{'Metric':<28} {'7K':>15} {'3K':>15} {'Union':>15}",
    ]
    for label, key in (("Total store bytes", "store"), ("Payload bytes", "payload"),
                       ("Reference metadata bytes", "metadata"), ("Offset index bytes", "offsets"),
                       ("Plain XOR vectors", "plain"), ("Reference-group vectors", "group"),
                       ("Build seconds", "build_seconds"), ("Verify seconds", "verify_seconds")):
        value = lambda store: f"{store[key]:,.2f}" if isinstance(store[key], float) else f"{store[key]:,}"
        lines.append(f"{label:<28} {value(large):>15} {value(small):>15} {value(union):>15}")
    lines.extend([
        "",
        f"Union saving vs standalone bytes={union['standalone'] - union['store']:,}; "
        f"pct={100 * (union['standalone'] - union['store']) / union['standalone']:.4f}",
        f"Union saving vs input CM6 bytes={union['source'] - union['store']:,}; "
        f"pct={100 * (union['source'] - union['store']) / union['source']:.4f}",
        f"Union minus 7K store bytes={union['store'] - large['store']:,}",
        f"Union minus 3K store bytes={union['store'] - small['store']:,}",
        "All three stores were fully decoded and exactly matched the CM6 oracle.",
    ])
    report = "\n".join(lines) + "\n"
    Path(args.out).write_text(report, encoding="utf-8")
    print(report, end="")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, AttributeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)
