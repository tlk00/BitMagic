cmake_minimum_required(VERSION 2.8)

project(bmjni)

find_package(JNI)

if (${APPLE})
    execute_process(COMMAND /usr/libexec/java_home OUTPUT_VARIABLE MACOS_JAVA_HOME OUTPUT_STRIP_TRAILING_WHITESPACE)
    set (PLATFORM_JNI_INCLUDES ${MACOS_JAVA_HOME}/include ${MACOS_JAVA_HOME}/include/darwin)
else()
    set (PLATFORM_JNI_INCLUDES ${JNI_INCLUDE_DIRS})
endif()

message(STATUS "JNI include: ${PLATFORM_JNI_INCLUDES}")

set( CMAKE_VERBOSE_MAKEFILE on )

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PLATFORM_JNI_INCLUDES}" "${PROJECT_SOURCE_DIR}/../../src" "${PROJECT_SOURCE_DIR}/include" "${PROJECT_SOURCE_DIR}/../libbm/include" "${PROJECT_SOURCE_DIR}/../libbm/src")


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(flags "-Wall -Wextra -march=core2 -g0 -O2 -msse4.2 -fPIC -fno-rtti -fno-exceptions ")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(flags, "-tpp7 -march=core2 -O2 -restrict -DBM_HASRESTRICT -fno-fnalias -Wall -fPIC")
    set(optflags, "-g0 -O3 -opt_report_fileopt.txt -opt_report_levelmax")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    set(flags "/W4 /EHsc")
    set(optflags "-O2")
endif()

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${flags}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flags} -std=c++11")


add_library(bmcpuidj SHARED "../libbm/src/libbmcpuid.c" "src/bmcpuidj.c")
add_library(bmjni SHARED "src/libbmjni.cpp" "src/bmjni.c")
add_library(bmjni-sse42 SHARED "src/libbmjni.cpp" "src/bmjni.c")
add_library(bmjni-avx2 SHARED "src/libbmjni.cpp" "src/bmjni.c")


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   set_target_properties(bmcpuidj     PROPERTIES COMPILE_FLAGS "-DBMDLLEXPORTS")
   set_target_properties(bmjni        PROPERTIES COMPILE_FLAGS "-DBMDLLEXPORTS")
   set_target_properties(bmjni-sse42  PROPERTIES COMPILE_FLAGS "-DBMDLLEXPORTS -DBMSSE42OPT")
   set_target_properties(bmjni-avx2   PROPERTIES COMPILE_FLAGS "-DBMDLLEXPORTS -DBMAVX2OPT")
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   set_target_properties(bmjni        PROPERTIES COMPILE_FLAGS "-march=core2")
   set_target_properties(bmjni-sse42  PROPERTIES COMPILE_FLAGS "-march=nehalem -msse4.2 -O2 -DBMSSE42OPT")
   set_target_properties(bmjni-avx2   PROPERTIES COMPILE_FLAGS "-march=native -mavx2 -O2 -DBMAVX2OPT")
endif()

target_link_libraries(bmcpuidj ${JNI_LIBRARIES})
target_link_libraries(bmjni ${JNI_LIBRARIES})
target_link_libraries(bmjni-sse42 ${JNI_LIBRARIES})
target_link_libraries(bmjni-avx2 ${JNI_LIBRARIES})





