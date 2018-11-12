#!/bin/sh

#  run_all_bg.sh
#  
#  Run all tests in background
#  Created by Anatoliy Kuznetsov on 9/2/17.
#

./stress_debug 2>&1 > debug.log &
PID_debug=$!

./stress_debug_sse2 2>&1 > debug_sse2.log &
PID_debug_sse2=$!

./stress_debug_sse42 2>&1 > debug_sse42.log &
PID_debug_sse42=$!

./stress_debug_avx2 2>&1 > debug_avx2.log &
PID_debug_avx2=$!


# ---------------------------------------------------------

wait $PID_debug_avx2
RET_debug_avx2=$?
if test "$RET_debug_avx2" != "0"; then
    echo "Error: debug AVX2 failed! " $RET_debug_avx2
else
    echo "debug AVX2 finished OK!"
fi

# ---------------------------------------------------------

wait $PID_debug_sse42
RET_debug_sse42=$?
if test "$RET_debug_sse42" != "0"; then
    echo "Error: debug SSE42 failed! " $RET_debug_sse42
else
    echo "debug SSE42 finished OK!"
fi

# ---------------------------------------------------------

wait $PID_debug_sse2
RET_debug_sse2=$?
if test "$RET_debug_sse2" != "0"; then
    echo "Error: debug SSE2 failed! " $RET_debug_sse2
else
    echo "debug SSE2 finished OK!"
fi

# ---------------------------------------------------------

wait $PID_debug
RET_debug=$?

if test "$RET_debug" != "0"; then
    echo "Error: debug failed! " $RET_debug
else
    echo "debug finished OK!"
fi

echo
echo "DONE"
