#!/bin/bash
#
# Bounded parallel stress-test runner for Unix-like systems.
# Compatible with the Bash 3.2 version included with macOS.

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

JOBS=""
RUN_NAME=""
LOG_DIR="logs/parallel"
FAIL_FAST=0
HEAVY_FIRST=0
TESTS=()

usage()
{
    cat <<EOF
Usage: $0 [runner options] [test selectors]

Runs the selected stress-test binary. A build name is required so the same
script can run platform-specific builds such as stress_release_avx2 on Linux
or stress_release_neon on Apple Silicon. With no test selectors, all
independent test groups are run in parallel.

Runner options:
  -j N, --jobs N        Run at most N test processes concurrently
  --build NAME          Use ./NAME as the test binary (required unless
                        --binary is used)
  -b BIN, --binary BIN  Use BIN as the test binary path/name (required unless
                        --build is used)
  -l DIR, --log-dir DIR Write logs below DIR (default: logs/parallel)
  --fail-fast           Stop launching new groups after the first failure
  --heavy-first         With the default test partition, launch heavier
                        groups first
  -h, --help            Display this help and exit

Test selectors from t.cpp:
  -ll,  -llevel         Low-level tests
  -s,   -support        Support-container tests
  -bvb0                 Bit-vector basic tests, part 0
  -bvb1                 Bit-vector basic tests, part 1
  -bvser                Bit-vector serialization tests
  -bvl0, -bvops0        Bit-vector logical operations, part 0
  -bvl1, -bvops1        Bit-vector logical operations, part 1
  -bvl2, -bvops2        Bit-vector logical operations, part 2
  -bvs,  -bvshift       Bit-vector shift/insert/erase tests
  -rc,   -rankc         Rank-compression tests
  -agg,  -aggregator    Aggregator tests
  -sv0                  Sparse-vector tests, part 0
  -sv1                  Sparse-vector tests, part 1
  -sv1a                 Sparse-vector tests, part 1a
  -sv1b                 Sparse-vector tests, part 1b
  -sv1c                 Sparse-vector tests, part 1c
  -sort, --sort         Sparse-vector sorting tests
  -csv0                 Compressed sparse-vector tests, part 0
  -csv0a                Compressed sparse-vector tests, part 0a
  -csv0b                Compressed sparse-vector tests, part 0b
  -csv0c                Compressed sparse-vector tests, part 0c
  -csv1                 Compressed sparse-vector tests, part 1
  -csv1a                Compressed sparse-vector tests, part 1a
  -csv1b                Compressed sparse-vector tests, part 1b
  -strsv, -svstr        String sparse-vector tests
  -cc                   Compressed-collection tests
  -ser                  All serialization tests
  -svf                  All floating-point sparse-vector tests
  -svf0                 Floating-point sparse-vector tests, part 0
  -svf0a                Floating-point sparse-vector tests, part 0a
  -svf0b                Floating-point sparse-vector tests, part 0b
  -svf0c                Floating-point sparse-vector tests, part 0c
  -svf1                 Floating-point sparse-vector tests, part 1

The broader overlapping selectors -bvb/-bvbasic, -bvo/-bvops/-bvl, -sv,
-csv/-rsc, -allsvser, and -svf are accepted for targeted runs but are not
included in the default partition because they duplicate work covered by
groups above.

Examples:
  ./run_macos_pl.sh -b bmtest -j 7 --fail-fast --heavy-first
  $0 --build stress_release
  $0 --build stress_release --heavy-first
  $0 --build stress_release_avx2 --jobs 6
  $0 --build stress_release_neon --jobs 4 -ll -bvb0 -bvb1
  $0 --binary ./stress_debug --jobs 3
EOF
}

is_valid_test()
{
    case "$1" in
        -ll|-llevel|-s|-support|-bvb|-bvbasic|-bvb0|-bvb1|-bvser|\
        -bvo|-bvops|-bvl|-bvl0|-bvops0|-bvl1|-bvops1|-bvl2|-bvops2|\
        -bvs|-bvshift|-rc|-rankc|-agg|-aggregator|\
        -sv|-sv0|-sv1|-sv1a|-sv1b|-sv1c|-sort|--sort|\
        -csv|-rsc|-csv0|-csv0a|-csv0b|-csv0c|-csv1|-csv1a|-csv1b|\
        -strsv|-svstr|-cc|-ser|-allsvser|-svf|-svf0|-svf0a|-svf0b|-svf0c|-svf1)
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
        --fail-fast)
            FAIL_FAST=1
            shift
            ;;
        --heavy-first)
            HEAVY_FIRST=1
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
    if [ "$HEAVY_FIRST" -ne 0 ]; then
        TESTS=(
            -svf0c -csv1a -csv0a -sv1c
            -bvb1 -bvs -strsv -ser
            -svf0b -bvl2 -svf1 -bvl0 -bvl1
            -sv1a -sv1b -csv0b -csv0c
            -ll -sv0 -bvb0 -bvser
            -agg -s -csv1b -sort -cc -rc
        )
    else
        TESTS=(
            -ll -s
            -bvb0 -bvb1 -bvser
            -bvl0 -bvl1 -bvl2
            -bvs -rc -agg
            -sv0 -sv1a -sv1b -sv1c -sort
            -csv0a -csv0c -csv0b -csv1a -csv1b
            -strsv -cc -ser
            -svf0a -svf0c -svf0b -svf1
        )
    fi
fi

RUN_LOG_DIR="${LOG_DIR}/${RUN_NAME}"
mkdir -p "$RUN_LOG_DIR" || exit 1
RUN_LOG_DIR="$(cd "$RUN_LOG_DIR" && pwd)" || exit 1

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

    # Every background process uses -silent to keep its log small.
    "$RUN_PATH" "$test_opt" -silent >"$log_file" 2>&1 &

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
if [ "$HEAVY_FIRST" -ne 0 ]; then
    echo "Schedule     : heavy-first"
else
    echo "Schedule     : natural"
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
