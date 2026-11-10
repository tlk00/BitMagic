#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$script_dir"

bminv=${BMINV:-./bminv}
source_file=${SOURCE:-gov2.unsorted.bv.cm6}
pair_index=${PAIR_INDEX:-gov2.unsorted.pair-shards.csv}
baseline=${BASELINE_PREFIX:-gov2.rs.cluster.regions3}
trial=${TRIAL_PREFIX:-gov2.rs.cluster.regions3.3k}
trial_vectors=${TRIAL_VECTORS:-3000}
trial_memory=${RS_CLUSTER_MEMORY_BYTES:-8589934592}
trial_min_bytes=${RS_MIN_BYTES:-4096}
trial_coverage=${RS_BYTE_COVERAGE:-90}
trial_seed=${RS_SEED:-42}
pair_threads=${PAIR_THREADS:-8}
xor_memory=${XOR_MEMORY_BYTES:-8589934592}
ref_top=${REF_TOP:-4}
report_out=${SCOPE_REPORT:-gov2.rs.cluster.scope.txt}

stage()
{
    printf '\n========================================================================\n'
    printf 'Stage %s/7: %s\n' "$1" "$2"
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

build_store()
{
    local prefix=$1
    local pair_csv=$prefix.pairs.csv
    local out=$prefix.xor.bv
    local refs=$prefix.xor.refs
    local offsets=$prefix.xor.offsets
    local report=$out.report.txt
    if [[ -e "$out" || -e "$refs" || -e "$offsets" || -e "$report" ]]; then
        require_complete "$report"
        grep -Fqx "Pair CSV input=$pair_csv" "$report" || {
            echo "ERROR: existing XOR store used a different pair CSV: $out" >&2; exit 1;
        }
        [[ -s "$out" && -s "$refs" && -s "$offsets" ]] || {
            echo "ERROR: incomplete XOR store at $prefix" >&2; exit 1;
        }
        echo "Completed XOR store already exists: $out"
        return
    fi
    if ! "$bminv" -xor-build -bvin "$source_file" -pair-index-in "$pair_index" \
        -pairs-in "$pair_csv" -xorout "$out" -xor-refs "$refs" \
        -xor-index "$offsets" -ref-top "$ref_top" -xor-memory "$xor_memory" \
        -level 6; then
        echo "ERROR: indexed XOR build failed for $prefix" >&2
        exit 1
    fi
    require_complete "$report"
}

verify_store()
{
    local prefix=$1
    local out=$prefix.xor.bv
    local refs=$prefix.xor.refs
    local offsets=$prefix.xor.offsets
    local log=$prefix.xor.verify.log
    if [[ -s "$log" ]] && grep -q '^XOR verification complete\.' "$log"; then
        echo "Completed verification already exists: $log"
        return
    fi
    if ! "$bminv" -xor-verify -bvin "$source_file" -pair-index-in "$pair_index" \
        -xorin "$out" -xor-refs "$refs" -xor-index "$offsets" \
        -xor-memory "$xor_memory" -level 6 | tee "$log"; then
        echo "ERROR: indexed XOR verification failed for $prefix" >&2
        exit 1
    fi
    grep -q '^XOR verification complete\.' "$log" || {
        echo "ERROR: verification did not report completion for $prefix" >&2; exit 1;
    }
}

stage 1 "Check the completed baseline search and source index"
command -v "$bminv" >/dev/null || { echo "ERROR: bminv missing: $bminv" >&2; exit 127; }
[[ -r "$source_file" && -r "$pair_index" ]] || {
    echo "ERROR: source or pair index missing" >&2; exit 1;
}
require_complete "$baseline.manifest.csv.report.txt"
require_complete "$baseline.pairs.csv.report.txt"
[[ -s "$baseline.pairs.csv" ]] || { echo "ERROR: baseline pair CSV missing" >&2; exit 1; }

stage 2 "Plan and search the smaller-cluster trial if needed"
if [[ -e "$trial.pairs.csv.report.txt" ]]; then
    require_complete "$trial.manifest.csv.report.txt"
    require_complete "$trial.pairs.csv.report.txt"
    [[ -s "$trial.pairs.csv" ]] || { echo "ERROR: trial pair CSV missing" >&2; exit 1; }
    echo "Completed trial search already exists: $trial.pairs.csv"
elif [[ -e "$trial.manifest.csv" || -e "$trial.manifest.csv.report.txt" ]]; then
    require_complete "$trial.manifest.csv.report.txt"
    [[ -s "$trial.manifest.csv" ]] || { echo "ERROR: trial manifest missing" >&2; exit 1; }
    [[ ! -e "$trial.pairs.csv" && ! -e "$trial.pairs.csv.best.csv" ]] || {
        echo "ERROR: partial trial search exists; remove only its incomplete outputs before retrying" >&2; exit 1;
    }
    "$bminv" -pair-search -bvin "$source_file" -pair-index-in "$pair_index" \
        -pair-cluster-in "$trial.manifest.csv" -csv "$trial.pairs.csv" \
        -pair-max-ratio 0 -pair-min-pct 5 -pair-min-gain 256 \
        -pair-threads "$pair_threads" -level 6 || {
        echo "ERROR: trial pair search failed" >&2; exit 1;
    }
else
    BMINV="$bminv" SOURCE="$source_file" PAIR_INDEX="$pair_index" \
        RS_CLUSTER_PREFIX="$trial" RS_CLUSTER_VECTORS="$trial_vectors" \
        RS_CLUSTER_MEMORY_BYTES="$trial_memory" RS_MIN_BYTES="$trial_min_bytes" \
        RS_BYTE_COVERAGE="$trial_coverage" RS_SEED="$trial_seed" \
        PAIR_THREADS="$pair_threads" sh ./run_gov2_xbv_rs_cluster_search.sh || {
        echo "ERROR: trial cluster plan/search failed" >&2; exit 1;
    }
fi
require_complete "$trial.pairs.csv.report.txt"

stage 3 "Build baseline XOR store from the completed 7K-cap search"
build_store "$baseline"
stage 4 "Decode and verify every baseline vector against CM6"
verify_store "$baseline"
stage 5 "Build trial XOR store from the smaller-cluster search"
build_store "$trial"
stage 6 "Decode and verify every trial vector against CM6"
verify_store "$trial"
stage 7 "Compare search coverage and measured whole-store compression"
python3 analyze_xbv_gov2_cluster_scope.py "$baseline" "$trial" \
    --pair-index "$pair_index" --out "$report_out"
echo "Large-scope experiment complete: $report_out"
