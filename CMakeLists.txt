cmake_minimum_required(VERSION 3.12)

project("bm" VERSION 7.4.0
             DESCRIPTION "BitMagic Library"
             HOMEPAGE_URL "https://github.com/tlk00/BitMagic")

set( CMAKE_VERBOSE_MAKEFILE on )

# SET(BMOPTFLAGS "BMSSE42OPT" CACHE STRING "User-specified optimization flag")

find_package(Threads REQUIRED)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}/src")

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src >
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

set(BM_HEADERS ${PROJECT_SOURCE_DIR}/src/bm.h 
               ${PROJECT_SOURCE_DIR}/src/bmblocks.h
               ${PROJECT_SOURCE_DIR}/src/bm3vl.h
               ${PROJECT_SOURCE_DIR}/src/bmaggregator.h
               ${PROJECT_SOURCE_DIR}/src/bmalgo.h
               ${PROJECT_SOURCE_DIR}/src/bmbvimport.h
               ${PROJECT_SOURCE_DIR}/src/bmalgo_impl.h
               ${PROJECT_SOURCE_DIR}/src/bmalgo_similarity.h
               ${PROJECT_SOURCE_DIR}/src/bmalloc.h
               ${PROJECT_SOURCE_DIR}/src/bmbmatrix.h
               ${PROJECT_SOURCE_DIR}/src/bmbuffer.h
               ${PROJECT_SOURCE_DIR}/src/bmconst.h
               ${PROJECT_SOURCE_DIR}/src/bmdbg.h
               ${PROJECT_SOURCE_DIR}/src/bmdef.h
               ${PROJECT_SOURCE_DIR}/src/bmfunc.h
               ${PROJECT_SOURCE_DIR}/src/bmfwd.h
               ${PROJECT_SOURCE_DIR}/src/bmgamma.h
               ${PROJECT_SOURCE_DIR}/src/bmintervals.h
               ${PROJECT_SOURCE_DIR}/src/bmrandom.h
               ${PROJECT_SOURCE_DIR}/src/bmrs.h
               ${PROJECT_SOURCE_DIR}/src/bmserial.h
               ${PROJECT_SOURCE_DIR}/src/bmsimd.h
               ${PROJECT_SOURCE_DIR}/src/bmsparsevec.h
               ${PROJECT_SOURCE_DIR}/src/bmsparsevec_algo.h
               ${PROJECT_SOURCE_DIR}/src/bmsparsevec_compr.h
               ${PROJECT_SOURCE_DIR}/src/bmsparsevec_parallel.h
               ${PROJECT_SOURCE_DIR}/src/bmsparsevec_serial.h
               ${PROJECT_SOURCE_DIR}/src/bmsparsevec_util.h
               ${PROJECT_SOURCE_DIR}/src/bmstrsparsevec.h
               ${PROJECT_SOURCE_DIR}/src/bmtask.h
               ${PROJECT_SOURCE_DIR}/src/bmthreadpool.h
               ${PROJECT_SOURCE_DIR}/src/bmtimer.h
               ${PROJECT_SOURCE_DIR}/src/bmtrans.h
               ${PROJECT_SOURCE_DIR}/src/bmundef.h
               ${PROJECT_SOURCE_DIR}/src/bmutil.h
               ${PROJECT_SOURCE_DIR}/src/bmxor.h
               ${PROJECT_SOURCE_DIR}/src/encoding.h
               ${PROJECT_SOURCE_DIR}/src/bmsse2.h
               ${PROJECT_SOURCE_DIR}/src/bmsse4.h
               ${PROJECT_SOURCE_DIR}/src/bmsse_util.h
               ${PROJECT_SOURCE_DIR}/src/bmavx2.h
               ${PROJECT_SOURCE_DIR}/src/bmavx512.h
               ${PROJECT_SOURCE_DIR}/src/bmbmi1.h
               ${PROJECT_SOURCE_DIR}/src/bmbmi2.h
)


set (CMAKE_CXX_STANDARD 17)


if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    set(warnopt "-Wno-psabi")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    if ("${BMOPTFLAGS}" STREQUAL "BMSSE42OPT")
	    set(bmoptf "-march=nehalem -O2 -msse4.2 -DBMSSE42OPT")
    elseif("${BMOPTFLAGS}" STREQUAL "BMAVX2OPT")
        set(bmoptf "-march=skylake -mavx2 -O2 -DBMAVX2OPT")
    else()
	if (CMAKE_BUILD_TYPE MATCHES "Release")
		set(bmoptf "-march=native -O2")
	else()
		set(bmoptf "-march=native")
	endif()
    endif()

    set(flags "-Wall -Wextra -Wno-ignored-qualifiers -fPIC -fomit-frame-pointer")
    set(optflags -g0 -O2)
    set(bmptflags -pthread)
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl,-lpthread")
    

elseif (CMAKE_CXX_COMPILER_ID MATCHES "Intel")
    set(flags, "-tpp7 -march=core2 -restrict -DBM_HASRESTRICT -fno-fnalias -Wall")
    set(optflags, "-g0 -O2 -opt_report_fileopt.txt -opt_report_levelmax")
elseif (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    if ("${BMOPTFLAGS}" STREQUAL "BMSSE42OPT")
	set(bmoptf "-DBMSSE42OPT")
    elseif("${BMOPTFLAGS}" STREQUAL "BMAVX2OPT")
        set(bmoptf "-DBMAVX2OPT") 
    endif()

    set(flags "/W4 /EHsc /F 5000000 ")
    set(optflags "-O2")
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} ")
endif()

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${flags} ${bmoptf} ${bmptflags}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flags} ${bmoptf} ${bmptflags} ${warnopt}")



add_executable(bmsvutil ${PROJECT_SOURCE_DIR}/utils/svutil/svutil.cpp)

add_executable(bmtest ${PROJECT_SOURCE_DIR}/tests/stress/t.cpp ${BM_HEADERS})
add_executable(bmptest ${PROJECT_SOURCE_DIR}/tests/test_threads/ptest.cpp ${BM_HEADERS})
add_executable(bmtest64 ${PROJECT_SOURCE_DIR}/tests/stress64/t64.cpp ${BM_HEADERS})
add_executable(bmperf ${PROJECT_SOURCE_DIR}/tests/perf/perf.cpp ${BM_HEADERS})
add_executable(bmperf64 ${PROJECT_SOURCE_DIR}/tests/perf64/perf64.cpp ${BM_HEADERS})
add_executable(bmlnkutil ${PROJECT_SOURCE_DIR}/utils/lnkutil/lnkutil.cpp ${BM_HEADERS})
add_executable(bminv ${PROJECT_SOURCE_DIR}/utils/inv_list/inv_list.cpp ${BM_HEADERS})

add_executable(bvsample01 ${PROJECT_SOURCE_DIR}/samples/bvsample01/sample1.cpp ${BM_HEADERS})
add_executable(bvsample02 ${PROJECT_SOURCE_DIR}/samples/bvsample02/sample2.cpp ${BM_HEADERS})
add_executable(bvsample03 ${PROJECT_SOURCE_DIR}/samples/bvsample03/sample3.cpp ${BM_HEADERS})
add_executable(bvsample04 ${PROJECT_SOURCE_DIR}/samples/bvsample04/sample4.cpp ${BM_HEADERS})
add_executable(bvsample05 ${PROJECT_SOURCE_DIR}/samples/bvsample05/sample5.cpp ${BM_HEADERS})
add_executable(bvsample06 ${PROJECT_SOURCE_DIR}/samples/bvsample06/sample6.cpp ${BM_HEADERS})
add_executable(bvsample07 ${PROJECT_SOURCE_DIR}/samples/bvsample07/sample7.cpp ${BM_HEADERS})
add_executable(bvsample08 ${PROJECT_SOURCE_DIR}/samples/bvsample08/sample8.cpp ${BM_HEADERS})
add_executable(bvsample09 ${PROJECT_SOURCE_DIR}/samples/bvsample09/sample9.cpp ${BM_HEADERS})
add_executable(bvsample10 ${PROJECT_SOURCE_DIR}/samples/bvsample10/sample10.cpp ${BM_HEADERS})
add_executable(bvsample11 ${PROJECT_SOURCE_DIR}/samples/bvsample11/sample11.cpp ${BM_HEADERS})
add_executable(bvsample12 ${PROJECT_SOURCE_DIR}/samples/bvsample12/sample12.cpp ${BM_HEADERS})
add_executable(bvsample14 ${PROJECT_SOURCE_DIR}/samples/bvsample14/sample14.cpp ${BM_HEADERS})
add_executable(bvsample15 ${PROJECT_SOURCE_DIR}/samples/bvsample15/sample15.cpp ${BM_HEADERS})
add_executable(bvsample16 ${PROJECT_SOURCE_DIR}/samples/bvsample16/sample16.cpp ${BM_HEADERS})
add_executable(bvsample17 ${PROJECT_SOURCE_DIR}/samples/bvsample17/sample17.cpp ${BM_HEADERS})
add_executable(bvsample18 ${PROJECT_SOURCE_DIR}/samples/bvsample18/sample18.cpp ${BM_HEADERS})
add_executable(bvsample18a ${PROJECT_SOURCE_DIR}/samples/bvsample18a/sample18a.cpp ${BM_HEADERS})
add_executable(bvsample19 ${PROJECT_SOURCE_DIR}/samples/bvsample19/sample19.cpp ${BM_HEADERS})
add_executable(bvsample20 ${PROJECT_SOURCE_DIR}/samples/bvsample20/sample20.cpp ${BM_HEADERS})
add_executable(bvsample21 ${PROJECT_SOURCE_DIR}/samples/bvsample21/sample21.cpp ${BM_HEADERS})
add_executable(bvsample22 ${PROJECT_SOURCE_DIR}/samples/bvsample22/sample22.cpp ${BM_HEADERS})
add_executable(bvsample23 ${PROJECT_SOURCE_DIR}/samples/bvsample23/sample23.cpp ${BM_HEADERS})
add_executable(bvsample24 ${PROJECT_SOURCE_DIR}/samples/bvsample24/sample24.cpp ${BM_HEADERS})
add_executable(bvsample25 ${PROJECT_SOURCE_DIR}/samples/bvsample25/sample25.cpp ${BM_HEADERS})
add_executable(bvsample26 ${PROJECT_SOURCE_DIR}/samples/bvsample26/sample26.cpp ${BM_HEADERS})

add_executable(bvsetalgebra ${PROJECT_SOURCE_DIR}/samples/bvsetalgebra/bvsetalgebra.cpp ${BM_HEADERS})
add_executable(bv3vlogic ${PROJECT_SOURCE_DIR}/samples/bv3vlogic/bv3vlogic.cpp ${BM_HEADERS})

add_executable(bvsample01_64 ${PROJECT_SOURCE_DIR}/samples/bvsample01_64/bvsample01_64.cpp ${BM_HEADERS})


add_executable(xsample01 ${PROJECT_SOURCE_DIR}/samples/xsample01/xsample01.cpp ${BM_HEADERS})
add_executable(xsample02 ${PROJECT_SOURCE_DIR}/samples/xsample02/xsample02.cpp ${BM_HEADERS})
add_executable(xsample03 ${PROJECT_SOURCE_DIR}/samples/xsample03/xsample03.cpp ${BM_HEADERS})
add_executable(xsample04 ${PROJECT_SOURCE_DIR}/samples/xsample04/xsample04.cpp ${BM_HEADERS})
add_executable(xsample04a ${PROJECT_SOURCE_DIR}/samples/xsample04a/xsample04a.cpp ${BM_HEADERS})
add_executable(xsample05 ${PROJECT_SOURCE_DIR}/samples/xsample05/xsample05.cpp ${BM_HEADERS})
add_executable(xsample06 ${PROJECT_SOURCE_DIR}/samples/xsample06/xsample06.cpp ${BM_HEADERS})
add_executable(xsample07 ${PROJECT_SOURCE_DIR}/samples/xsample07/xsample07.cpp ${BM_HEADERS})
add_executable(xsample07a ${PROJECT_SOURCE_DIR}/samples/xsample07a/xsample07a.cpp ${BM_HEADERS})
add_executable(xsample08 ${PROJECT_SOURCE_DIR}/samples/xsample08/xsample08.cpp ${BM_HEADERS})
add_executable(xsample09 ${PROJECT_SOURCE_DIR}/samples/xsample09/xsample09.cpp ${BM_HEADERS})
add_executable(xsample10 ${PROJECT_SOURCE_DIR}/samples/xsample10/xsample10.cpp ${BM_HEADERS})

add_executable(svsample01 ${PROJECT_SOURCE_DIR}/samples/svsample01/svsample01.cpp ${BM_HEADERS})
add_executable(svsample02 ${PROJECT_SOURCE_DIR}/samples/svsample02/svsample02.cpp ${BM_HEADERS})
add_executable(svsample03 ${PROJECT_SOURCE_DIR}/samples/svsample03/svsample03.cpp ${BM_HEADERS})
add_executable(svsample04 ${PROJECT_SOURCE_DIR}/samples/svsample04/svsample04.cpp ${BM_HEADERS})
add_executable(svsample05 ${PROJECT_SOURCE_DIR}/samples/svsample05/svsample05.cpp ${BM_HEADERS})
add_executable(svsample06 ${PROJECT_SOURCE_DIR}/samples/svsample06/svsample06.cpp ${BM_HEADERS})
add_executable(svsample07 ${PROJECT_SOURCE_DIR}/samples/svsample07/svsample07.cpp ${BM_HEADERS})
add_executable(svsample07a ${PROJECT_SOURCE_DIR}/samples/svsample07a/svsample07a.cpp ${BM_HEADERS})
add_executable(svsample08 ${PROJECT_SOURCE_DIR}/samples/svsample08/svsample08.cpp ${BM_HEADERS})
add_executable(svsample09 ${PROJECT_SOURCE_DIR}/samples/svsample09/svsample09.cpp ${BM_HEADERS})
add_executable(svsample10 ${PROJECT_SOURCE_DIR}/samples/svsample10/svsample10.cpp ${BM_HEADERS})

add_executable(rscsample01 ${PROJECT_SOURCE_DIR}/samples/rscsample01/rscsample01.cpp ${BM_HEADERS})
add_executable(rscsample02 ${PROJECT_SOURCE_DIR}/samples/rscsample02/rscsample02.cpp ${BM_HEADERS})
add_executable(rscsample03 ${PROJECT_SOURCE_DIR}/samples/rscsample03/rscsample03.cpp ${BM_HEADERS})
add_executable(rscsample04 ${PROJECT_SOURCE_DIR}/samples/rscsample04/rscsample04.cpp ${BM_HEADERS})
add_executable(rscsample05 ${PROJECT_SOURCE_DIR}/samples/rscsample05/rscsample05.cpp ${BM_HEADERS})
add_executable(rscsample06 ${PROJECT_SOURCE_DIR}/samples/rscsample06/rscsample06.cpp ${BM_HEADERS})

add_executable(strsvsample01 ${PROJECT_SOURCE_DIR}/samples/strsvsample01/strsvsample01.cpp ${BM_HEADERS})
add_executable(strsvsample02 ${PROJECT_SOURCE_DIR}/samples/strsvsample02/strsvsample02.cpp ${BM_HEADERS})
add_executable(strsvsample02a ${PROJECT_SOURCE_DIR}/samples/strsvsample02a/strsvsample02a.cpp ${BM_HEADERS})
add_executable(strsvsample03 ${PROJECT_SOURCE_DIR}/samples/strsvsample03/strsvsample03.cpp ${BM_HEADERS})
add_executable(strsvsample04 ${PROJECT_SOURCE_DIR}/samples/strsvsample04/strsvsample04.cpp ${BM_HEADERS})
add_executable(strsvsample05 ${PROJECT_SOURCE_DIR}/samples/strsvsample05/strsvsample05.cpp ${BM_HEADERS})
add_executable(strsvsample06 ${PROJECT_SOURCE_DIR}/samples/strsvsample06/strsvsample06.cpp ${BM_HEADERS})
add_executable(strsvsample07 ${PROJECT_SOURCE_DIR}/samples/strsvsample07/strsvsample07.cpp ${BM_HEADERS})
add_executable(strsvsample08 ${PROJECT_SOURCE_DIR}/samples/strsvsample08/strsvsample08.cpp ${BM_HEADERS})


