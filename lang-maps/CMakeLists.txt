cmake_minimum_required(VERSION 2.8)

project(bmlangmap)


set( CMAKE_VERBOSE_MAKEFILE on )

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}/../src" "${PROJECT_SOURCE_DIR}/libbm/include")

MESSAGE( STATUS "CMAKE_CXX_COMPILER_ID:         " ${CMAKE_CXX_COMPILER_ID} )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(LINKER_FLAGS "-nodefaultlibs -lc")
    set(LD_FLAGS ${LINKER_FLAGS})
    set(LDFLAGS ${LINKER_FLAGS})
    set(CMAKE_SHARED_LINKER_FLAGS  ${LINKER_FLAGS})
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
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flags} ")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ldflags}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ldflags}")

add_library(libbm-static STATIC "libbm/src/libbm.cpp")

add_executable(libbmtest ${PROJECT_SOURCE_DIR}/test/libbmtest.c)

#set_property(TARGET target PROPERTY CMAKE_SHARED_LINKER_FLAGS “${CMAKE_SHARED_LINKER} -nodefaultlibs -lc”)
#set_property(TARGET libbmtest PROPERTY CMAKE_SHARED_LINKER_FLAGS “${CMAKE_SHARED_LINKER} -nodefaultlibs -lc”)
#set_property(TARGET libbmtest PROPERTY CMAKE_EXE_LINKER_FLAGS “${CMAKE_EXE_LINKER_FLAGS} ${LINKER_FLAGS}”)
target_link_libraries(libbmtest libbm-static ${LINKER_FLAGS})

MESSAGE( STATUS "LINKER_FLAGS:              " ${LINKER_FLAGS} )
MESSAGE( STATUS "CMAKE_EXE_LINKER_FLAGS:    " ${CMAKE_EXE_LINKER_FLAGS} )


