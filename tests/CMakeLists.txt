# Shared internal test support, also listed explicitly for generated IDE projects.
set(BM_TEST_COMMON_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/common/debug_allocator.h
    ${CMAKE_CURRENT_SOURCE_DIR}/common/gena.h
    ${CMAKE_CURRENT_SOURCE_DIR}/common/neon_test.h
    ${CMAKE_CURRENT_SOURCE_DIR}/common/serialization_allocator_tests.h
    ${CMAKE_CURRENT_SOURCE_DIR}/common/sparse_serial_levels_tests.h
    ${CMAKE_CURRENT_SOURCE_DIR}/common/stacktrace_dbg.h
    ${CMAKE_CURRENT_SOURCE_DIR}/common/test_util.h
)
add_library(bm_test_common INTERFACE)
target_sources(bm_test_common INTERFACE ${BM_TEST_COMMON_HEADERS})
target_include_directories(bm_test_common INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/common)
source_group("Common Test Headers" FILES ${BM_TEST_COMMON_HEADERS})

add_executable(bmtest stress/t.cpp stress/rlebtv.h stress/serial_transform_tests.h ${BM_HEADERS})
add_executable(bmtest64 stress64/t64.cpp ${BM_HEADERS})
add_executable(bmptest test_threads/ptest.cpp ${BM_HEADERS})
add_executable(bmperf perf/perf.cpp perf/simd_bench.h ${BM_HEADERS})
add_executable(bmperf64 perf64/perf64.cpp ${BM_HEADERS})
foreach(test_target bmtest bmtest64 bmperf64)
    target_link_libraries(${test_target} PRIVATE bm_test_common)
endforeach()

add_subdirectory(regression)
