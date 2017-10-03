cmake_minimum_required(VERSION 2.8)

project(bmjni)


set( CMAKE_VERBOSE_MAKEFILE on )

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}/../../src" "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/../libbm/include" "${PROJECT_SOURCE_DIR}/../libbm/src" "/System/Library/Frameworks/JavaVM.framework/Versions/Current/Headers")


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(flags "-Wall -Wextra -march=core2 -msse4.2 -fPIC -fno-rtti -fno-exceptions ")
    set(optflags, -g0 -O3)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(flags, "-tpp7 -march=core2 -restrict -DBM_HASRESTRICT -fno-fnalias -Wall -fPIC")
    set(optflags, "-g0 -O3 -opt_report_fileopt.txt -opt_report_levelmax")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(flags "/W4 /EHsc")
    set(optflags "-O2")
endif()

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${flags}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flags}")


add_library(bmjni SHARED "src/libbmjni.cpp" "src/bmjni.cpp")




