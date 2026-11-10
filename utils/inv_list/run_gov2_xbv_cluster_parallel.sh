#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$script_dir"

bminv=${BMINV:-./bminv}
source_file=${SOURCE:-gov2.unsorted.bv.cm6}
pair_index=${PAIR_INDEX:-gov2.unsorted.pair-shards.csv}
manifest=${MANIFEST:-gov2.rs.cluster.regions3.3k.manifest.csv}
control=${CONTROL_CSV:-gov2.rs.cluster.regions3.3k.pairs.csv}
parallel=${PARALLEL_CSV:-gov2.rs.cluster.regions3.3k.parallel-24g-9t.pairs.csv}
threads=${PAIR_THREADS:-9}
memory=${PAIR_CLUSTER_MEMORY_BYTES:-25769803776}
comparison=${COMPARISON_REPORT:-gov2.rs.cluster.regions3.3k.parallel-24g-9t.compare.txt}

stage()
{
    printf '\n========================================================================\n'
    printf 'Stage %s/3: %s\n' "$1" "$2"
    printf '========================================================================\n\n'
}

require_complete()
{
    local path=$1
    if [[ ! -s "$path" ]] || ! grep -qx 'Status: complete' "$path"; then
        echo "ERROR: missing or incomplete report: $path" >&2
        exit 1
    fi
}

stage 1 'Validate the collection, 3K cluster manifest, and completed control run'
command -v "$bminv" >/dev/null || { echo "ERROR: bminv missing: $bminv" >&2; exit 127; }
[[ -r "$source_file" && -r "$pair_index" && -r "$manifest" ]] || {
    echo 'ERROR: source, pair index, or cluster manifest missing' >&2; exit 1;
}
require_complete "$manifest.report.txt"
require_complete "$control.report.txt"
[[ -s "$control" && -s "$control.best.csv" ]] || {
    echo 'ERROR: completed control pair CSVs missing' >&2; exit 1;
}
[[ "$parallel" != "$control" ]] || { echo 'ERROR: parallel output must differ from control' >&2; exit 1; }

stage 2 'Run one-cluster-per-task pair search with a total active-memory budget'
if [[ -e "$parallel" || -e "$parallel.best.csv" || -e "$parallel.report.txt" ]]; then
    require_complete "$parallel.report.txt"
    [[ -s "$parallel" && -s "$parallel.best.csv" ]] || {
        echo 'ERROR: parallel output is incomplete; inspect it before retrying' >&2; exit 1;
    }
    echo "Completed parallel search already exists: $parallel"
else
    if ! "$bminv" -pair-search -bvin "$source_file" -pair-index-in "$pair_index" \
        -pair-cluster-in "$manifest" -pair-cluster-parallel \
        -pair-cluster-memory "$memory" -csv "$parallel" \
        -pair-max-ratio 0 -pair-min-pct 5 -pair-min-gain 256 \
        -pair-threads "$threads" -level 6; then
        echo 'ERROR: cluster-parallel search failed; partial outputs retained for inspection' >&2
        exit 1
    fi
    require_complete "$parallel.report.txt"
fi

stage 3 'Compare exact CSV results and elapsed search time with the control'
if ! cmp -s "$control" "$parallel"; then
    echo "ERROR: pair CSV differs from control: $parallel" >&2
    exit 1
fi
if ! cmp -s "$control.best.csv" "$parallel.best.csv"; then
    echo "ERROR: best-pair CSV differs from control: $parallel.best.csv" >&2
    exit 1
fi

python3 - "$control.report.txt" "$parallel.report.txt" "$comparison" <<'PY'
import re
import sys

def read_report(path):
    with open(path, encoding='utf-8') as stream:
        text = stream.read()
    match = re.search(r'^Elapsed seconds=([0-9.]+)$', text, re.MULTILINE)
    if not match:
        raise SystemExit(f'ERROR: elapsed time missing from {path}')
    return float(match.group(1))

control = read_report(sys.argv[1])
parallel = read_report(sys.argv[2])
lines = [
    'Pair and best-pair CSVs: byte-identical',
    f'Control elapsed seconds={control:.6f}',
    f'Cluster-parallel elapsed seconds={parallel:.6f}',
    f'Search speedup={control / parallel:.6f}x',
]
with open(sys.argv[3], 'w', encoding='utf-8') as stream:
    stream.write('\n'.join(lines) + '\n')
print('\n'.join(lines))
print(f'Comparison report: {sys.argv[3]}')
PY
