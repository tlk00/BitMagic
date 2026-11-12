#!/bin/bash
#
# Bounded parallel 64-bit stress-test runner for Unix-like systems.
# Compatible with the Bash 3.2 version included with macOS.

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

JOBS=""
RUN_NAME=""
LOG_DIR="logs/parallel"
ORDER_FILE=""
FAIL_FAST=0
PEDANTIC_ARGS=()
FAST=0
ONLY_STRESS=0
NO_STRESS=0
TESTS=()

usage()
{
    cat <<EOF
Usage: $0 [runner options] [test selectors]

  -fast, --fast        Reduce stress checks; keep fixture sizes
                      Incompatible with -pedantic / --pedantic.
  -pedantic, --pedantic
                      Sweep bvector and sparse-vector serialization levels 0-10
                      Bvector pairs: 11 same-level + 1 mixed case.
                      Sparse vectors: 11 levels; full float bookmark variants at level 6.
                      Sampled interval-16 float bookmarks at levels 7-10.
                      Default: rotate explicit levels, non-exhaustive.

Runs the selected 64-bit stress-test binary. A build name is required so the same
script can run platform-specific builds such as stress64_release_avx2 on Linux
or stress64_release_neon on Apple Silicon. With no test selectors, all
independent test groups are run in parallel.

Runner options:
  -j N, --jobs N        Run at most N test processes concurrently
  --build NAME          Use ./NAME as the test binary (required unless
                        --binary is used)
  -b BIN, --binary BIN  Use BIN as the test binary path/name (required unless
                        --build is used)
  -l DIR, --log-dir DIR Write logs below DIR (default: logs/parallel)
  --fail-fast           Stop launching new groups after the first failure
  --heavy-first         Use learned longest-first order (also the default).
                        New selectors run first; missing selectors are ignored.
                        History: ./logs/parallel/order.tsv, shared by all builds.
                        Completed non-fast jobs update history, including failed runs.
  --order-file FILE     Use a separate learned timing history at FILE
  -onlystress,
    --only-stress       Pass -onlystress to every t64.cpp test process
  -nostress,
    --no-stress         Pass -nostress to every t64.cpp test process
  -silent, --silent     Accepted for t64.cpp compatibility; workers are silent
  -h, --help            Display this help and exit

Test selectors verified against t64.cpp:
  -bvb,  -bvbasic       Bit-vector basic tests (includes bvector serialization)
  -bvb0, -bvb1a/b/c     Independent bit-vector basic subgroups
  -complexity           Vector/XOR complexity all-pairs validation
  -bvser, -bvset        Bit-vector serialization tests
  -fileser              File serialization tests
  -svser                Sparse-vector level/opcode/index/gather tests
  -svfileser            Sparse-vector file and stringstream tests
  -strsvstream          Disk-free string streaming and gather tests
  -svfstream            Float file/stringstream serialization and gather tests
  -svfsteam             Alias for -svfstream
  -svindex              Sparse-vector index persistence tests
  -bvl,  -bvo, -bvops   All bit-vector logical operations
  -bvl0, -bvops0        Logical-operation correctness tests
  -bvl0a/b              AND/OR correctness / XOR/SUB correctness
  -bvl1, -bvops1        OR stress test
  -bvl2, -bvops2        AND stress test
  -bvl3, -bvops3        SUB stress test
  -bvl4, -bvops4        XOR stress test
  -bvl1a/b ... -bvl4a/b Two trial batches per logical operation
  -bvs,  -bvshift       Bit-vector shift/insert/erase tests
  -bvs0/1/2/3           Direct shifts / shift aggregator / insert / erase
  -rc,   -rankc         Rank-compression tests
  -agg,  -aggregator    Aggregator tests
  -sv                   Sparse-vector tests
  -sv0, -sv1a/b/c       Core / range scan / signed scan / sorting groups
  -csv                  Compressed sparse-vector tests
  -csv0a0               Compressed sparse-vector basic checks
  -csv0a3, -csv0a4      Full set/null stress / load stress
  -csv0a30/31/32/33    Disjoint randomized set/null stress shards
  -csv0b, -csv1a/b      Signed and gather / GT scan / serial and scan groups
  -strsv, -svstr        String sparse-vector tests
  -cc                   Compressed-collection tests
  -svf                  Floating-point sparse-vector tests
  -svf0                 Float basic/scanner tests
  -svf0a/b              Float basic / unbounded scanner tests
  -svf0c1, -svf0c3     Linear / skewed float scanner tests
  -svf0c20/21/22/23    Random float scanner query partitions
  -svf1                 Float/RSC index/gather/scanner aggregate
  -svf1a0/1/2          Float index/gather levels 0-5 / 6 / 7-10
  -svf1b0/1/2          RSC float index/gather levels 0-5 / 6 / 7-10
  -svf1c               RSC float scanner tests
  -svf1a, -svf1b       Runner groups for plain / RSC float index and gather
  -svserquick           Runner group for sparse serialization, index, and streams

Parsed by t64.cpp but currently without an independently dispatched test block:
  -ll, -llevel, -s, -support
These compatibility selectors are accepted here but omitted from the default
partition. Use -svf for current float coverage.

Examples:
  ./run_all_pl.sh -b stress64_release_neon -j 7 --pedantic -svf1a0 -svf1a1 -svf1a2 -svf1b0 -svf1b1 -svf1b2 -svf1c
  ./run_all_pl.sh -b stress64_release_neon --pedantic -svser -svfileser -svfstream -strsvstream
  ./run_all_pl.sh -b stress64_release_neon -j 1 --pedantic -bvl1
  ./run_all_pl.sh -b bmtest64 -j 7 --fail-fast --heavy-first
  ./run_all_pl.sh --build stress64_release
  ./run_all_pl.sh --build stress64_release --heavy-first
  ./run_all_pl.sh --build stress64_release_avx2 --jobs 6
  ./run_all_pl.sh --binary ./stress64_debug --jobs 3
  ./run_all_pl.sh --build stress64_release --no-stress
EOF
}

is_valid_test()
{
    case "$1" in
        -ll|-llevel|-s|-support|-bvb|-bvbasic|-bvb0|-bvb1[a-c]|-complexity|-bvser|-bvset|\
        -fileser|-svfileser|-svser|-svserquick|-strsvstream|-svfstream|-svindex|\
        -bvo|-bvops|-bvl|-bvl0|-bvops0|-bvl0[a-b]|-bvl1|-bvops1|-bvl1[a-b]|-bvl2|-bvops2|-bvl2[a-b]|\
        -bvl3|-bvops3|-bvl3[a-b]|-bvl4|-bvops4|-bvl4[a-b]|-bvs|-bvshift|-bvs[0-3]|-rc|-rankc|-agg|-aggregator|\
        -sv|-sv0|-sv1[a-c]|-csv|-csv0a0|-csv0a3|-csv0a4|-csv0a3[0-3]|-csv0b|-csv1[a-b]|-strsv|-svstr|-cc|\
        -svf|-svf0|-svf0a|-svf0b|-svf0c1|-svf0c3|-svf0c2[0-3]|-svf1|-svf1a|-svf1b|-svf1a0|-svf1a1|-svf1a2|-svf1b0|-svf1b1|-svf1b2|-svf1c)
            return 0
            ;;
        *)
            return 1
            ;;
    esac
}

while [ "$#" -gt 0 ]; do
    case "$1" in
        -j|--jobs)
            [ "$#" -ge 2 ] || {
                echo "Missing value for $1" >&2
                exit 2
            }
            JOBS="$2"
            shift 2
            ;;
        --build)
            [ "$#" -ge 2 ] || {
                echo "Missing value for $1" >&2
                exit 2
            }
            RUN_NAME="$2"
            shift 2
            ;;
        -b|--binary)
            [ "$#" -ge 2 ] || {
                echo "Missing value for $1" >&2
                exit 2
            }
            RUN_NAME="$2"
            shift 2
            ;;
        -l|--log-dir)
            [ "$#" -ge 2 ] || {
                echo "Missing value for $1" >&2
                exit 2
            }
            LOG_DIR="$2"
            shift 2
            ;;
        -fast|--fast)
            FAST=1
            shift
            ;;
        -pedantic|--pedantic)
            PEDANTIC_ARGS=(--pedantic)
            shift
            ;;
        --fail-fast)
            FAIL_FAST=1
            shift
            ;;
        --heavy-first)
            # Compatibility option: learned scheduling is always enabled.
            shift
            ;;
        --order-file)
            [ "$#" -ge 2 ] || { echo "Missing value for $1" >&2; exit 2; }
            ORDER_FILE="$2"
            shift 2
            ;;
        -onlystress|--only-stress)
            ONLY_STRESS=1
            shift
            ;;
        -nostress|--no-stress)
            NO_STRESS=1
            shift
            ;;
        -silent|--silent)
            shift
            ;;
        -svfsteam)
            TESTS[${#TESTS[@]}]="-svfstream"
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            if is_valid_test "$1"; then
                TESTS[${#TESTS[@]}]="$1"
                shift
            else
                echo "Unknown option or test selector: $1" >&2
                usage >&2
                exit 2
            fi
            ;;
    esac
done

if [ "$ONLY_STRESS" -ne 0 ] && [ "$NO_STRESS" -ne 0 ]; then
    echo "--only-stress and --no-stress are mutually exclusive" >&2
    exit 2
fi

if [ "$FAST" -ne 0 ] && [ "${#PEDANTIC_ARGS[@]}" -ne 0 ]; then
    echo "--fast and --pedantic are mutually exclusive" >&2
    exit 2
fi

ARCH="$(uname -m)"

if [ -z "$RUN_NAME" ]; then
    echo "Missing required --build NAME or --binary BIN" >&2
    usage >&2
    exit 2
fi

case "$RUN_NAME" in
    */*) RUN_PATH="$RUN_NAME" ;;
    *)   RUN_PATH="./$RUN_NAME" ;;
esac

if [ ! -x "$RUN_PATH" ]; then
    echo "Test binary is missing or not executable: $RUN_PATH" >&2
    exit 2
fi

if [ -z "$JOBS" ]; then
    JOBS="$(sysctl -n hw.logicalcpu 2>/dev/null)"
    [ -n "$JOBS" ] || JOBS="$(nproc 2>/dev/null)"
    [ -n "$JOBS" ] || JOBS=4
    [ "$JOBS" -le 8 ] || JOBS=8
fi

case "$JOBS" in
    ''|*[!0-9]*|0)
        echo "Invalid job count: $JOBS" >&2
        exit 2
        ;;
esac

if [ "${#TESTS[@]}" -eq 0 ]; then
    # Initial longest-first guess for the untimed 64-bit partitions. Completed
    # normal runs replace this order with measured weights from order.tsv.
    TESTS=(
        -csv0a30 -csv0a31 -csv0a32 -csv0a33 -svf0c20 -svf0c21 -svf0c22 -svf0c23
        -bvl1a -bvl2a -bvl3a -bvl4a -bvl1b -bvl2b -bvl3b -bvl4b
        -sv1b -bvb1b -bvb1c -bvs0 -csv1b -csv0a4
        -bvb0 -bvb1a -bvser -svserquick -bvl0a -bvl0b -bvs1 -bvs2 -bvs3
        -rc -agg -sv0 -sv1a -sv1c -csv0a0 -csv0b -csv1a -strsv -cc
        -svf0a -svf0b -svf0c1 -svf0c3 -svf1a -svf1b
    )
fi

RUN_LOG_DIR="${LOG_DIR}/${RUN_NAME}"
mkdir -p "$RUN_LOG_DIR" || exit 1
RUN_LOG_DIR="$(cd "$RUN_LOG_DIR" && pwd)" || exit 1

# One learned order per suite, independent of binary and log-directory overrides.
# Rows contain selector and latest elapsed seconds, sorted longest first.
[ -n "$ORDER_FILE" ] || ORDER_FILE="${SCRIPT_DIR}/logs/parallel/order.tsv"
# Fast timings are not representative: use existing history read-only.
ORDER_INPUT="$ORDER_FILE"
if [ "$FAST" -eq 0 ]; then
    mkdir -p "$(dirname "$ORDER_FILE")" || exit 1
    touch "$ORDER_FILE" || exit 1
elif [ ! -f "$ORDER_FILE" ]; then
    ORDER_INPUT=/dev/null
fi

# Unknown selectors precede timed selectors; ties preserve selection order.
ORDERED_TESTS=()
while IFS= read -r test_opt; do
    ORDERED_TESTS[${#ORDERED_TESTS[@]}]="$test_opt"
done < <(
    printf '%s\n' "${TESTS[@]}" |
        awk 'FILENAME == ARGV[1] {
                 if (NF == 2 && $2 ~ /^[0-9]+$/) duration[$1] = $2
                 next
             }
             { known = ($1 in duration)
               print known, (known ? duration[$1] : 0), FNR, $1 }' \
            "$ORDER_INPUT" - |
        LC_ALL=C sort -k1,1n -k2,2nr -k3,3n |
        awk '{ print $4 }'
)
TESTS=("${ORDERED_TESTS[@]}")

# Persist each completed job so a later interruption retains earlier results.
# Keep other selectors for future targeted runs; scheduling filters them out.
record_timing()
{
    # Keep elapsed times in the report, but never learn weights from fast runs.
    [ "$FAST" -eq 0 ] || return 0
    local selector="$1" duration="$2" temporary
    temporary="$(mktemp "${ORDER_FILE}.XXXXXX")" || return 1
    if awk -v selector="$selector" -v duration="$duration" '
        NF == 2 && $2 ~ /^[0-9]+$/ && $1 != selector { print $1, $2 }
        END { print selector, duration }' "$ORDER_FILE" |
        LC_ALL=C sort -k2,2nr -k1,1 > "$temporary" &&
        mv "$temporary" "$ORDER_FILE"; then
        return 0
    fi
    rm -f "$temporary"
    echo "Warning: could not update timing history: $ORDER_FILE" >&2
    return 1
}

PIDS=()
PID_TESTS=()
PID_LOGS=()
PID_STARTS=()
PID_SECONDS=()
PID_STATUS=()
FAILED_TESTS=()
FAILED_LOGS=()
FAILED_STATUS=()
FAILURES=0
STOP_LAUNCHING=0
RUNNING_COUNT=0

format_duration()
{
    seconds="$1"
    hours=$((seconds / 3600))
    minutes=$(((seconds % 3600) / 60))
    secs=$((seconds % 60))
    printf "%02d:%02d:%02d" "$hours" "$minutes" "$secs"
}

print_timing_table()
{
    echo "============================================================"
    echo "Test timings (launch order)"
    echo "============================================================"
    printf "%-12s %-10s %-12s %s\n" "Test" "Status" "Elapsed" "Log"
    printf "%-12s %-10s %-12s %s\n" "----" "------" "-------" "---"

    index=0
    while [ "$index" -lt "${#PID_TESTS[@]}" ]; do
        elapsed="${PID_SECONDS[$index]}"
        [ -n "$elapsed" ] || elapsed=0

        status="${PID_STATUS[$index]}"
        [ -n "$status" ] || status="not-run"

        printf "%-12s %-10s %-12s %s\n" \
            "${PID_TESTS[$index]}" \
            "$status" \
            "$(format_duration "$elapsed")" \
            "${PID_LOGS[$index]}"

        index=$((index + 1))
    done
}

start_test()
{
    test_opt="$1"
    test_name="${test_opt#-}"
    test_name="${test_name#-}"
    log_file="${RUN_LOG_DIR}/${test_name}.log"

    echo "START  $test_opt -> $log_file"

    case "$test_opt" in
        -svserquick)
            worker_args=(-svfileser -svser -svindex -strsvstream -svfstream -silent)
            ;;
        -svf1a)
            worker_args=(-svf1a0 -svf1a1 -svf1a2 -silent)
            ;;
        -svf1b)
            worker_args=(-svf1b0 -svf1b1 -svf1b2 -svf1c -silent)
            ;;
        *)
            worker_args=("$test_opt" -silent)
            ;;
    esac
    if [ "${#PEDANTIC_ARGS[@]}" -ne 0 ]; then
        worker_args[${#worker_args[@]}]=--pedantic
    fi
    if [ "$FAST" -ne 0 ]; then
        worker_args[${#worker_args[@]}]=--fast
    fi

    if [ "$ONLY_STRESS" -ne 0 ]; then
        "$RUN_PATH" "${worker_args[@]}" -onlystress >"$log_file" 2>&1 &
    elif [ "$NO_STRESS" -ne 0 ]; then
        "$RUN_PATH" "${worker_args[@]}" -nostress >"$log_file" 2>&1 &
    else
        "$RUN_PATH" "${worker_args[@]}" >"$log_file" 2>&1 &
    fi

    index=${#PIDS[@]}
    PIDS[$index]=$!
    PID_TESTS[$index]="$test_opt"
    PID_LOGS[$index]="$log_file"
    PID_STARTS[$index]="$(date +%s)"
    PID_SECONDS[$index]=""
    PID_STATUS[$index]="running"
}

reap_finished()
{
    index=0

    while [ "$index" -lt "${#PIDS[@]}" ]; do
        pid="${PIDS[$index]}"

        if [ -n "$pid" ] && ! kill -0 "$pid" 2>/dev/null; then
            wait "$pid"
            status=$?
            elapsed=$(($(date +%s) - PID_STARTS[$index]))
            PID_SECONDS[$index]="$elapsed"
            record_timing "${PID_TESTS[$index]}" "$elapsed"

            if [ "$status" -eq 0 ]; then
                PID_STATUS[$index]="pass"
                echo "PASS   ${PID_TESTS[$index]} ($(format_duration "$elapsed"))"
            else
                PID_STATUS[$index]="fail:$status"
                echo "FAIL   ${PID_TESTS[$index]} (exit $status, $(format_duration "$elapsed"))"
                echo "       Details: ${PID_LOGS[$index]}"

                failure_index=${#FAILED_TESTS[@]}
                FAILED_TESTS[$failure_index]="${PID_TESTS[$index]}"
                FAILED_LOGS[$failure_index]="${PID_LOGS[$index]}"
                FAILED_STATUS[$failure_index]="$status"
                FAILURES=$((FAILURES + 1))

                if [ "$FAIL_FAST" -ne 0 ]; then
                    STOP_LAUNCHING=1
                fi
            fi

            PIDS[$index]=""
        fi

        index=$((index + 1))
    done
}

count_running()
{
    count=0
    index=0

    while [ "$index" -lt "${#PIDS[@]}" ]; do
        [ -n "${PIDS[$index]}" ] && count=$((count + 1))
        index=$((index + 1))
    done

    RUNNING_COUNT=$count
}

echo "Architecture : $ARCH"
echo "Binary       : $RUN_PATH"
echo "Parallel jobs: $JOBS"
echo "Schedule     : learned heavy-first (new jobs first)"
echo "Order file   : $ORDER_FILE"
if [ "$ONLY_STRESS" -ne 0 ]; then
    echo "Test mode    : only stress loops"
elif [ "$NO_STRESS" -ne 0 ]; then
    echo "Test mode    : no stress loops"
fi
if [ "$FAST" -ne 0 ]; then
    echo "Workload     : fast (reduced stress checks; timing weights not updated)"
fi
if [ "${#PEDANTIC_ARGS[@]}" -ne 0 ]; then
    echo "Serialization: pedantic bvector/sparse-vector levels 0-10 (bvector pairs add one mixed case)"
else
    echo "Warning: THIS IS NON EXHAUSTIVE TEST: use --pedantic for exhaustive compression-level coverage."
fi
echo "Test groups  : ${#TESTS[@]}"
echo "Logs         : $RUN_LOG_DIR"
echo

test_index=0

while [ "$test_index" -lt "${#TESTS[@]}" ]; do
    [ "$STOP_LAUNCHING" -ne 0 ] && break

    count_running
    if [ "$RUNNING_COUNT" -lt "$JOBS" ]; then
        start_test "${TESTS[$test_index]}"
        test_index=$((test_index + 1))
    else
        reap_finished
        sleep 0.1
    fi
done

while :; do
    count_running
    [ "$RUNNING_COUNT" -eq 0 ] && break
    reap_finished
    sleep 0.1
done

echo

print_timing_table
echo

if [ "$FAILURES" -ne 0 ]; then
    echo "============================================================"
    echo "FAILED: $FAILURES test group(s)"
    echo "============================================================"

    index=0
    while [ "$index" -lt "${#FAILED_TESTS[@]}" ]; do
        echo "${FAILED_TESTS[$index]} exited ${FAILED_STATUS[$index]}"
        echo "Log: ${FAILED_LOGS[$index]}"
        echo "-------------------- log tail --------------------"
        tail -20 "${FAILED_LOGS[$index]}"
        echo "--------------------------------------------------"
        index=$((index + 1))
    done
    exit 1
fi

if [ "$STOP_LAUNCHING" -ne 0 ]; then
    echo "Stopped before all groups were launched."
    exit 1
fi

echo "All ${#TESTS[@]} test groups passed."
exit 0
