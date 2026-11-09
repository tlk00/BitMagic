# Standalone focused regressions. Keep assertions active in Release builds too.
foreach(regression debug_allocator_test serial_compression_levels
        serial_gap_normalization serial_opcode_registry serial_released_compat)
    add_executable(bm_${regression} ${regression}.cpp ${BM_HEADERS})
    target_link_libraries(bm_${regression} PRIVATE bm_test_common)
    if(MSVC)
        target_compile_options(bm_${regression} PRIVATE /UNDEBUG)
    else()
        target_compile_options(bm_${regression} PRIVATE -UNDEBUG)
    endif()
    # Released compatibility requires an external fixture directory argument.
    if(BUILD_TESTING AND NOT regression STREQUAL "serial_released_compat")
        add_test(NAME ${regression} COMMAND bm_${regression})
    endif()
endforeach()
