cmake_minimum_required(VERSION 2.8...3.13)

project(bmlangmap)


set( CMAKE_VERBOSE_MAKEFILE on )
set (CMAKE_CXX_STANDARD 17)


set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/build/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

include_directories("${PROJECT_SOURCE_DIR}/../src" "${PROJECT_SOURCE_DIR}/libbm/include")

MESSAGE( STATUS "CMAKE_CXX_COMPILER_ID:         " ${CMAKE_CXX_COMPILER_ID} )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")

    if ("${BMOPTFLAGS}" STREQUAL "BMSSE42OPT")
        set(bmoptf "-march=nehalem -O2 -msse4.2 -DBMSSE42OPT")
    elseif("${BMOPTFLAGS}" STREQUAL "BMAVX2OPT")
        set(bmoptf "-march=skylake -mavx2 -O2 -DBMAVX2OPT")
    else()
        #set(bmoptf "-march=core2")
    endif()

    set(LINKER_FLAGS "-nodefaultlibs -lc")
    set(LD_FLAGS ${LINKER_FLAGS})
    set(LDFLAGS ${LINKER_FLAGS})
    set(CMAKE_SHARED_LINKER_FLAGS  ${LINKER_FLAGS})
    set(flags "-Wall -Wextra -fPIC -fno-rtti -fno-exceptions ")
    set(optflags "-g0 -O2")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(flags, "-tpp7 -march=core2 -restrict -DBM_HASRESTRICT -fno-fnalias -Wall -fPIC")
    set(optflags, "-g0 -O3 -opt_report_fileopt.txt -opt_report_levelmax")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    if ("${BMOPTFLAGS}" STREQUAL "BMSSE42OPT")
	set(bmoptf "-DBMSSE42OPT")
    elseif("${BMOPTFLAGS}" STREQUAL "BMAVX2OPT")
        set(bmoptf "-DBMAVX2OPT") 
    endif()

    set(flags "/W4 /EHsc /F 5000000 ")
    set(optflags "-O2")
    SET( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} /STACK:\"10000000\"")

endif()

set(CMAKE_C_FLAGS   "${CMAKE_C_FLAGS} ${flags} ${bmoptf}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flags} ${optflags} ${bmoptf}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ldflags}")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ldflags}")

add_library(bm-static STATIC "libbm/src/libbm.cpp")
add_library(bm-dll SHARED "libbm/src/libbm.cpp")
add_library(bmcpuid SHARED "libbm/src/libbmcpuid.c")
add_library(bm-dll-sse42 SHARED "libbm/src/libbm.cpp")
add_library(bm-dll-avx2 SHARED "libbm/src/libbm.cpp")


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
   set_target_properties(bmcpuid      PROPERTIES COMPILE_FLAGS "-DBMDLLEXPORTS")
   set_target_properties(bm-dll       PROPERTIES COMPILE_FLAGS "-DBMDLLEXPORTS")
   set_target_properties(bm-dll-sse42 PROPERTIES COMPILE_FLAGS "-DBMDLLEXPORTS -DBMSSE42OPT")
   set_target_properties(bm-dll-avx2  PROPERTIES COMPILE_FLAGS "-DBMDLLEXPORTS -DBMAVX2OPT")
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
   set_target_properties(bm-dll       PROPERTIES COMPILE_FLAGS "-march=core2")
   set_target_properties(bm-dll-sse42 PROPERTIES COMPILE_FLAGS "-march=nehalem -msse4.2 -DBMSSE42OPT")
   set_target_properties(bm-dll-avx2  PROPERTIES COMPILE_FLAGS "-march=skylake -mavx2 -DBMAVX2OPT")
endif()



add_executable(libbmtest ${PROJECT_SOURCE_DIR}/test/libbmtest.c)

#set_property(TARGET target PROPERTY CMAKE_SHARED_LINKER_FLAGS “${CMAKE_SHARED_LINKER} -nodefaultlibs -lc”)
#set_property(TARGET libbmtest PROPERTY CMAKE_SHARED_LINKER_FLAGS “${CMAKE_SHARED_LINKER} -nodefaultlibs -lc”)
#set_property(TARGET libbmtest PROPERTY CMAKE_EXE_LINKER_FLAGS “${CMAKE_EXE_LINKER_FLAGS} ${LINKER_FLAGS}”)
target_link_libraries(libbmtest bm-static ${LINKER_FLAGS})

MESSAGE( STATUS "LINKER_FLAGS:              " ${LINKER_FLAGS} )
MESSAGE( STATUS "CMAKE_EXE_LINKER_FLAGS:    " ${CMAKE_EXE_LINKER_FLAGS} )


