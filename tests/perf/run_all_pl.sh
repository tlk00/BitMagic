#!/bin/bash
# Run independent perf groups in bounded parallel processes (Bash 3.2 compatible).

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)" || exit 1
cd "$SCRIPT_DIR" || exit 1

usage()
{
    cat <<EOF
Usage: $0 --build NAME [-j N] [options] [group ...]
       $0 --binary PATH [-j N] [options] [group ...]

With no group names, run every group reported by the perf binary.
Each group has its own process, working directory, and timing log.

  -j N, --jobs N        Maximum concurrent processes (default 2)
  --build NAME          Binary in this directory, e.g. perf_release_neon
  -b PATH, --binary PATH  Binary path, relative to this directory or absolute
  -l DIR, --log-dir DIR  Log root (default logs/parallel)
  --fast                Pass --fast to each group; do not update timing weights
  --fail-fast           Stop launching groups after a failure; drain active jobs
  --order-file FILE     Normal-mode timing history (default logs/parallel/order.tsv)
  -h, --help            Show this help

Examples:
  ./run_all_pl.sh --build perf_release -j 4
  ./run_all_pl.sh --build perf_release_neon -j 2 --fast
  ./run_all_pl.sh --build perf_release -j 2 codec-decode float-rsc mmap
EOF
}

JOBS=2
RUN_NAME=""
LOG_DIR="logs/parallel"
ORDER_FILE=""
FAST=0
FAIL_FAST=0
SELECTED=()

while [ "$#" -gt 0 ]; do
    case "$1" in
        -j|--jobs)
            [ "$#" -ge 2 ] || { echo "Missing value for $1" >&2; exit 2; }
            JOBS="$2"
            shift 2
            ;;
        --build|-b|--binary)
            [ "$#" -ge 2 ] || { echo "Missing value for $1" >&2; exit 2; }
            RUN_NAME="$2"
            shift 2
            ;;
        -l|--log-dir)
            [ "$#" -ge 2 ] || { echo "Missing value for $1" >&2; exit 2; }
            LOG_DIR="$2"
            shift 2
            ;;
        --order-file)
            [ "$#" -ge 2 ] || { echo "Missing value for $1" >&2; exit 2; }
            ORDER_FILE="$2"
            shift 2
            ;;
        --fast)
            FAST=1
            shift
            ;;
        --fail-fast)
            FAIL_FAST=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -* )
            echo "Unknown option: $1" >&2
            exit 2
            ;;
        *)
            SELECTED[${#SELECTED[@]}]="$1"
            shift
            ;;
    esac
done

case "$JOBS" in
    ''|*[!0-9]*|0) echo "Invalid job count: $JOBS" >&2; exit 2 ;;
esac
[ -n "$RUN_NAME" ] || { echo "Specify --build or --binary" >&2; exit 2; }

case "$RUN_NAME" in
    /*) RUN_PATH="$RUN_NAME" ;;
    *) RUN_PATH="$SCRIPT_DIR/${RUN_NAME#./}" ;;
esac
[ -x "$RUN_PATH" ] || { echo "Binary is missing or not executable: $RUN_PATH" >&2; exit 2; }

GROUP_LIST="$("$RUN_PATH" --list-groups)" || exit 2
[ -n "$GROUP_LIST" ] || { echo "Binary does not list perf groups: $RUN_PATH" >&2; exit 2; }
ALL_GROUPS=()
while IFS= read -r group; do
    [ -n "$group" ] && ALL_GROUPS[${#ALL_GROUPS[@]}]="$group"
done <<< "$GROUP_LIST"

if [ "${#SELECTED[@]}" -eq 0 ]; then
    SELECTED=("${ALL_GROUPS[@]}")
fi
for group in "${SELECTED[@]}"; do
    valid=0
    for known in "${ALL_GROUPS[@]}"; do
        [ "$group" != "$known" ] || valid=1
    done
    [ "$valid" -eq 1 ] || { echo "Unknown perf group: $group" >&2; exit 2; }
done
for ((i=0; i<${#SELECTED[@]}; ++i)); do
    for ((j=i+1; j<${#SELECTED[@]}; ++j)); do
        [ "${SELECTED[$i]}" != "${SELECTED[$j]}" ] || {
            echo "Duplicate perf group: ${SELECTED[$i]}" >&2
            exit 2
        }
    done
done

mkdir -p "$LOG_DIR" || exit 1
LOG_DIR="$(cd "$LOG_DIR" && pwd)" || exit 1
RUN_DIR="$(mktemp -d "$LOG_DIR/run.XXXXXX")" || exit 1
[ -n "$ORDER_FILE" ] || ORDER_FILE="$LOG_DIR/order.tsv"
case "$ORDER_FILE" in
    /*) ;;
    *) ORDER_FILE="$SCRIPT_DIR/$ORDER_FILE" ;;
esac
if [ "$FAST" -eq 0 ]; then
    mkdir -p "$(dirname "$ORDER_FILE")" || exit 1
    touch "$ORDER_FILE" || exit 1
fi

# Initial normal-mode weights are rough seconds from the x86 full-run log.
# Learned group times replace them after a normal run.
default_weight()
{
    case "$1" in
        codec-decode) echo 510 ;;
        codec-encode|float-rsc|float-compare|bit-enum) echo 420 ;;
        float-plain|sparse-scan|sparse-access|bit-logic) echo 360 ;;
        bit-core|bit-shift|sparse-serial|float-serial|rsc|string|mmap) echo 300 ;;
        *) echo 0 ;;
    esac
}

schedule_file="$(mktemp "$RUN_DIR/schedule.XXXXXX")" || exit 1
for group in "${SELECTED[@]}"; do
    weight=0
    if [ -f "$ORDER_FILE" ]; then
        weight="$(awk -v name="$group" '$1 == name && $2 ~ /^[0-9]+$/ { n=$2 } END { print n+0 }' "$ORDER_FILE")"
    fi
    [ "$weight" -gt 0 ] || weight="$(default_weight "$group")"
    printf '%s\t%s\n' "$weight" "$group"
done | LC_ALL=C sort -k1,1nr -k2,2 > "$schedule_file"
ORDERED=()
while IFS="$(printf '\t')" read -r weight group; do
    [ -n "$group" ] && ORDERED[${#ORDERED[@]}]="$group"
done < "$schedule_file"
[ "${#ORDERED[@]}" -eq "${#SELECTED[@]}" ] || { echo "Incomplete schedule" >&2; exit 1; }

record_timing()
{
    [ "$FAST" -eq 0 ] || return 0
    local name="$1" seconds="$2" temporary
    temporary="$(mktemp "${ORDER_FILE}.XXXXXX")" || return 1
    if awk -v name="$name" -v seconds="$seconds" '
        NF == 2 && $1 != name { print $1, $2 }
        END { print name, seconds }' "$ORDER_FILE" > "$temporary" &&
        mv "$temporary" "$ORDER_FILE"; then
        return 0
    fi
    rm -f "$temporary"
    return 1
}

PIDS=()
NAMES=()
STARTS=()
STATUSES=()
ELAPSED=()
NEXT=0
RUNNING=0
FAILURES=0
STOP_LAUNCHING=0

stop_workers()
{
    for pid in "${PIDS[@]}"; do
        [ -z "$pid" ] || kill -TERM "$pid" 2>/dev/null || true
    done
}
trap 'stop_workers; exit 130' INT
trap 'stop_workers; exit 143' TERM
trap 'stop_workers; exit 129' HUP

launch_group()
{
    local name="$1" index=${#PIDS[@]} work_dir="$RUN_DIR/$1"
    mkdir -p "$work_dir" || exit 1
    echo "START $name -> $RUN_DIR/$name.log"
    if [ "$FAST" -eq 1 ]; then
        (cd "$work_dir" && "$RUN_PATH" --group "$name" --fast) > "$RUN_DIR/$name.log" 2>&1 &
    else
        (cd "$work_dir" && "$RUN_PATH" --group "$name") > "$RUN_DIR/$name.log" 2>&1 &
    fi
    PIDS[$index]=$!
    NAMES[$index]="$name"
    STARTS[$index]="$(date +%s)"
    STATUSES[$index]="running"
    ELAPSED[$index]=0
    RUNNING=$((RUNNING + 1))
}

reap_finished()
{
    local index pid status seconds
    for ((index=0; index<${#PIDS[@]}; ++index)); do
        pid="${PIDS[$index]}"
        if [ -n "$pid" ] && ! kill -0 "$pid" 2>/dev/null; then
            wait "$pid"
            status=$?
            seconds=$(( $(date +%s) - STARTS[$index] ))
            ELAPSED[$index]=$seconds
            record_timing "${NAMES[$index]}" "$seconds" || true
            if [ "$status" -eq 0 ]; then
                STATUSES[$index]=pass
                echo "PASS  ${NAMES[$index]} (${seconds}s)"
            else
                STATUSES[$index]="fail:$status"
                FAILURES=$((FAILURES + 1))
                echo "FAIL  ${NAMES[$index]} (exit $status, ${seconds}s): $RUN_DIR/${NAMES[$index]}.log"
                [ "$FAIL_FAST" -eq 0 ] || STOP_LAUNCHING=1
            fi
            PIDS[$index]=""
            RUNNING=$((RUNNING - 1))
        fi
    done
}

echo "Binary: $RUN_PATH"
echo "Jobs: $JOBS"
echo "Groups: ${#ORDERED[@]}"
echo "Logs: $RUN_DIR"
echo "Schedule: ${ORDERED[*]}"

while [ "$NEXT" -lt "${#ORDERED[@]}" ] && [ "$STOP_LAUNCHING" -eq 0 ]; do
    while [ "$RUNNING" -lt "$JOBS" ] && [ "$NEXT" -lt "${#ORDERED[@]}" ]; do
        launch_group "${ORDERED[$NEXT]}"
        NEXT=$((NEXT + 1))
    done
    reap_finished
    [ "$RUNNING" -eq 0 ] || sleep 0.2
done
while [ "$RUNNING" -gt 0 ]; do
    reap_finished
    [ "$RUNNING" -eq 0 ] || sleep 0.2
done

echo
printf '%-18s %-10s %-9s %s\n' Group Status Seconds Log
for ((i=0; i<${#NAMES[@]}; ++i)); do
    printf '%-18s %-10s %-9s %s\n' "${NAMES[$i]}" "${STATUSES[$i]}" \
        "${ELAPSED[$i]}" "$RUN_DIR/${NAMES[$i]}.log"
done
if [ "$FAILURES" -ne 0 ] || [ "$STOP_LAUNCHING" -ne 0 ]; then
    echo "FAILED: $FAILURES group(s); $NEXT/${#ORDERED[@]} launched"
    exit 1
fi
echo "All ${#ORDERED[@]} perf groups passed."
