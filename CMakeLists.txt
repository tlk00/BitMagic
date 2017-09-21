cmake_minimum_required(VERSION 2.8)

project(bm)


set( CMAKE_VERBOSE_MAKEFILE on )

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}/src")


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(flags "-Wall -Wextra -march=core2 -msse4.2 -fPIC")
    set(optflags, -g0 -O3)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(flags, "-tpp7 -march=core2 -restrict -DBM_HASRESTRICT -fno-fnalias -Wall")
    set(optflags, "-g0 -O3 -opt_report_fileopt.txt -opt_report_levelmax")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(flags "/W4 /EHsc")
    set(optflags "-O2")
endif()

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${flags}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flags}")

add_executable(bmtest ${PROJECT_SOURCE_DIR}/tests/stress/t.cpp)
add_executable(bmperf ${PROJECT_SOURCE_DIR}/tests/perf/perf.cpp)

add_executable(bvsample01 ${PROJECT_SOURCE_DIR}/samples/bvsample01/sample1.cpp)
add_executable(bvsample02 ${PROJECT_SOURCE_DIR}/samples/bvsample02/sample2.cpp)
add_executable(bvsample03 ${PROJECT_SOURCE_DIR}/samples/bvsample03/sample3.cpp)
add_executable(bvsample04 ${PROJECT_SOURCE_DIR}/samples/bvsample04/sample4.cpp)
add_executable(bvsample05 ${PROJECT_SOURCE_DIR}/samples/bvsample05/sample5.cpp)
add_executable(bvsample06 ${PROJECT_SOURCE_DIR}/samples/bvsample06/sample6.cpp)
add_executable(bvsample08 ${PROJECT_SOURCE_DIR}/samples/bvsample08/sample8.cpp)
add_executable(bvsample09 ${PROJECT_SOURCE_DIR}/samples/bvsample09/sample9.cpp)
add_executable(bvsample10 ${PROJECT_SOURCE_DIR}/samples/bvsample10/sample10.cpp)

add_executable(svsample01 ${PROJECT_SOURCE_DIR}/samples/svsample01/svsample01.cpp)
add_executable(svsample02 ${PROJECT_SOURCE_DIR}/samples/svsample02/svsample02.cpp)
add_executable(svsample03 ${PROJECT_SOURCE_DIR}/samples/svsample03/svsample03.cpp)


