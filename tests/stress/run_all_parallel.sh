#!/bin/bash

#  run_all_parallel.sh
#  

# read cmd-line to identify which build to run
# as:
#  run_all_parallel.sh sse42

#BUILD_CLASS="debug"

BUILD_CLASS="release"
BUILD_TYPE="$1"

case $BUILD_TYPE in
    sse2)
      RUN_NAME="stress_${BUILD_CLASS}_${BUILD_TYPE}"
      ;;
    sse42)
      RUN_NAME="stress_${BUILD_CLASS}_${BUILD_TYPE}"
      ;;
    gcc)
      RUN_NAME="stress_${BUILD_CLASS}_${BUILD_TYPE}"
      ;;
    avx2)
      RUN_NAME="stress_${BUILD_CLASS}_${BUILD_TYPE}"
      ;;
    avx512)
      RUN_NAME="stress_${BUILD_CLASS}_${BUILD_TYPE}"
      ;;
    neon)
      RUN_NAME="stress_${BUILD_CLASS}_${BUILD_TYPE}"
      ;;
    *)    # unknown option select "release" build by default
      BUILD_TYPE="release"
      RUN_NAME="stress_${BUILD_TYPE}"
      ;;
esac

if [ -f "$RUN_NAME" ]; then
    echo $RUN_NAME
else
    echo "$RUN_NAME does not exist."
    exit 1;
fi

rm -f *.log

PID_ARR=()
declare -a OPT_MAP

TEST_OPT="-ll -s -bvb0 -bvb1 -bvser -bvl0 -bvl1 -bvl2 -bvs -rc -agg -sv0 "
for OPT in $TEST_OPT;
do
    LOG_OPT=${RUN_NAME}${OPT}.log
    echo ${RUN_NAME} ${OPT} ">" ${LOG_OPT}

    ./${RUN_NAME} ${OPT} -silent &> ${LOG_OPT} &
    PID=$!
    PID_ARR+=($PID)
    OPT_ARR+=($LOG_OPT)
    echo "pid="${PID}
    OPT_MAP[$PID]=${LOG_OPT}
done

#echo ${OPT_MAP[*]}
echo "----------------------"
echo "waiting for test cases.."

ERR_FLAG="0"

for PID in ${PID_ARR[*]};
do
    LOG=${OPT_MAP[$PID]}
    echo ${PID} ${LOG}
    wait ${PID}
    RET=$?

    if test "$RET" != "0"; then
        ERR_FLAG="1"
        echo
        echo "Error: test case failed! log file:" $LOG
        echo "==============================================================="
        echo "-------------------------------------------------------------->8"
        tail -10 $LOG
        echo "==============================================================="
        echo "-------------------------------------------------------------->8"
        echo
    else
        echo "$LOG OK!"
    fi
    echo
    echo
done


if test "$ERR_FLAG" != "0"; then
    echo "======================"
    echo "BUILD=" ${BUILD_TYPE}
    echo "ERROR(s) detected!"
    exit 1
fi

unset PID_ARR
unset OPT_ARR
unset OPT_MAP

PID_ARR=()
declare -a OPT_MAP

TEST_OPT="-sv1 -sort -csv0 -csv1 -strsv -cc -ser -allsvser"
for OPT in $TEST_OPT;
do
    LOG_OPT=${RUN_NAME}${OPT}.log
    echo ${RUN_NAME} ${OPT} ">" ${LOG_OPT}

    ./${RUN_NAME} ${OPT} -silent &> ${LOG_OPT} &
    PID=$!
    PID_ARR+=($PID)
    OPT_ARR+=($LOG_OPT)
    echo "pid="${PID}
    OPT_MAP[$PID]=${LOG_OPT}
done

#echo ${OPT_MAP[*]}
echo "----------------------"
echo "waiting for test cases.."

ERR_FLAG="0"

for PID in ${PID_ARR[*]};
do
    LOG=${OPT_MAP[$PID]}
    echo ${PID} ${LOG}
    wait ${PID}
    RET=$?

    if test "$RET" != "0"; then
        ERR_FLAG="1"
        echo
        echo "Error: test case failed! log file:" $LOG
        echo "==============================================================="
        echo "-------------------------------------------------------------->8"
        tail -10 $LOG
        echo "==============================================================="
        echo "-------------------------------------------------------------->8"
        echo
    else
        echo "$LOG OK!"
    fi
    echo
    echo
done


if test "$ERR_FLAG" != "0"; then
    echo "======================"
    echo "BUILD=" ${BUILD_TYPE}
    echo "ERROR(s) detected!"
    exit 1
else
    echo "DONE     " ${RUN_NAME}
fi

