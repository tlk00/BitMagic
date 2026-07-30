#!/bin/bash
#
# Bounded parallel stress-test runner for macOS.
# Compatible with the Bash 3.2 version included with macOS.

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

JOBS=""
RUN_NAME=""
LOG_DIR="logs/macos"
FAIL_FAST=0
TESTS=()

usage()
{
    cat <<EOF
Usage: $0 [runner options] [test selectors]

Runs the release NEON stress-test binary by default on Apple Silicon. Debug
builds are never selected automatically; use --binary stress_debug explicitly.
With no test selectors, all independent test groups are run in parallel.

Runner options:
  -j N, --jobs N        Run at most N test processes concurrently
  -b BIN, --binary BIN  Use BIN instead of the architecture-specific release
                        default (stress_release_neon on Apple Silicon)
  -l DIR, --log-dir DIR Write logs below DIR (default: logs/macos)
  --fail-fast           Stop launching new groups after the first failure
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
  -sort, --sort         Sparse-vector sorting tests
  -csv0                 Compressed sparse-vector tests, part 0
  -csv1                 Compressed sparse-vector tests, part 1
  -strsv, -svstr        String sparse-vector tests
  -cc                   Compressed-collection tests
  -ser                  All serialization tests
  -svf                  Floating-point sparse-vector tests

The broader overlapping selectors -bvb/-bvbasic, -bvo/-bvops/-bvl, -sv,
-csv/-rsc, and -allsvser are accepted for targeted runs but are not included
in the default partition because they duplicate work covered by groups above.

Examples:
  $0
  $0 --jobs 6
  $0 --jobs 4 -ll -bvb0 -bvb1
  $0 --binary stress_release --jobs 6
  $0 --binary stress_debug --jobs 3
EOF
}

is_valid_test()
{
    case "$1" in
        -ll|-llevel|-s|-support|-bvb|-bvbasic|-bvb0|-bvb1|-bvser|\
        -bvo|-bvops|-bvl|-bvl0|-bvops0|-bvl1|-bvops1|-bvl2|-bvops2|\
        -bvs|-bvshift|-rc|-rankc|-agg|-aggregator|\
        -sv|-sv0|-sv1|-sort|--sort|-csv|-rsc|-csv0|-csv1|\
        -strsv|-svstr|-cc|-ser|-allsvser|-svf)
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
    case "$ARCH" in
        arm64|aarch64)
            if [ -x ./stress_release_neon ]; then
                RUN_NAME="stress_release_neon"
            else
                RUN_NAME="stress_release"
            fi
            ;;
        x86_64)
            if [ -x ./stress_release_avx2 ]; then
                RUN_NAME="stress_release_avx2"
            else
                RUN_NAME="stress_release"
            fi
            ;;
        *)
            echo "Unsupported architecture: $ARCH" >&2
            exit 2
            ;;
    esac
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
    TESTS=(
        -ll -s
        -bvb0 -bvb1 -bvser
        -bvl0 -bvl1 -bvl2
        -bvs -rc -agg
        -sv0 -sv1 -sort
        -csv0 -csv1
        -strsv -cc -ser -svf
    )
fi

RUN_LOG_DIR="${LOG_DIR}/${RUN_NAME}"
mkdir -p "$RUN_LOG_DIR" || exit 1
RUN_LOG_DIR="$(cd "$RUN_LOG_DIR" && pwd)" || exit 1

PIDS=()
PID_TESTS=()
PID_LOGS=()
FAILED_TESTS=()
FAILED_LOGS=()
FAILED_STATUS=()
FAILURES=0
STOP_LAUNCHING=0
RUNNING_COUNT=0

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
}

reap_finished()
{
    index=0

    while [ "$index" -lt "${#PIDS[@]}" ]; do
        pid="${PIDS[$index]}"

        if [ -n "$pid" ] && ! kill -0 "$pid" 2>/dev/null; then
            wait "$pid"
            status=$?

            if [ "$status" -eq 0 ]; then
                echo "PASS   ${PID_TESTS[$index]}"
            else
                echo "FAIL   ${PID_TESTS[$index]} (exit $status)"
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
