#!/bin/sh

#  run_all_ser_bg.sh
#  
#  Run all serialization tests in background
#  Created by Anatoliy Kuznetsov on 9/2/17.
#


./stress_release -sv -cvs -strsv 2>&1 > release.log &
PID_release=$!

./stress_release_gcc -sv -cvs -strsv 2>&1 > release_gcc.log &
PID_release_gcc=$!

./stress_release_sse2 -sv -cvs -strsv 2>&1 > release_sse2.log &
PID_release_sse2=$!

./stress_release_sse42 -sv -cvs -strsv 2>&1 > release_sse42.log &
PID_release_sse42=$!

./stress_release_avx2 -sv -cvs -strsv 2>&1 > release_avx2.log &
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

wait $PID_release_sse2
RET_release_sse2=$?
if test "$RET_release_sse2" != "0"; then
    echo "Error: Release SSE2 failed! " $RET_release_sse2
else
    echo "Release SSE2 finished OK!"
fi

# ---------------------------------------------------------

wait $PID_release_gcc
RET_release_gcc=$?

if test "$RET_release_gcc" != "0"; then
    echo "Error: Release GCC failed! " $RET_release
else
    echo "Release GCC finished OK!"
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
