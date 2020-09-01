#!/bin/sh

#  run_all_bg.sh
#  
#  Run all tests in background
#  Created by Anatoliy Kuznetsov on 9/2/17.
#

#echo "DEBUG"
#./stress_debug || exit 1

./stress64_release 2>&1 > release.log &
PID_release=$!

./stress64_release_sse42 2>&1 > release_sse42.log &
PID_release_sse42=$!

./stress64_release_avx2 2>&1 > release_avx2.log &
PID_release_avx2=$!


# ---------------------------------------------------------

wait $PID_release_avx2
RET_release_avx2=$?
if test "$RET_release_avx2" != "0"; then
    echo "Error: Release AVX2 failed! " $RET_release_avx2
else
    echo "Release AVX2 finished OK!"
fi

# ---------------------------------------------------------

wait $PID_release_sse42
RET_release_sse42=$?
if test "$RET_release_sse42" != "0"; then
    echo "Error: Release SSE42 failed! " $RET_release_sse42
else
    echo "Release SSE42 finished OK!"
fi

# ---------------------------------------------------------


wait $PID_release
RET_release=$?

if test "$RET_release" != "0"; then
    echo "Error: Release failed! " $RET_release
else
    echo "Release finished OK!"
fi

echo
echo "DONE"
