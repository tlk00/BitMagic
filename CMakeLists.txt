cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(bm VERSION 3.8.0 LANGUAGES CXX)

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}/src")



add_executable(bmtest ${PROJECT_SOURCE_DIR}/tests/stress/t.cpp)
add_executable(bmperf ${PROJECT_SOURCE_DIR}/tests/perf/perf.cpp)
