# Test layout

The root contains build files, scripts and documentation. Keep C/C++ sources in
the appropriate subdirectory:

- `common/`: shared internal test headers, fixtures, allocators and generators.
- `regression/`: standalone focused regression programs, each with its own main.
- `stress/` and `stress64/`: 32-/64-bit addressing stress drivers and local helpers.
- `perf/` and `perf64/`: performance drivers and local benchmark helpers.
- `test_threads/` and `bench-tpch/`: specialized test programs.

Common headers may depend on other common headers and library headers. Keep
helpers used by just one driver with that driver. Some existing fixture headers
expect types declared by their including driver; they are not standalone APIs.

## CMake

From the repository root, configure an out-of-tree build and select the targets
needed for the check. Native Mac NEON builds can use:

```sh
cmake -S . -B /tmp/bm-tests -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS=-DBMNEONOPT
cmake --build /tmp/bm-tests --target bmtest bmtest64 bm_debug_allocator_test
```

Executables are written to the build directory's `bin/`. Existing stress and
performance target names are preserved. Shared headers appear in the generated
IDE projects under `Common Test Headers`.

The standalone targets are `bm_debug_allocator_test`,
`bm_serial_compression_levels`, `bm_serial_gap_normalization`,
`bm_serial_opcode_registry`, and `bm_serial_released_compat`.
The first four are registered with CTest; build them before running:

```sh
ctest --test-dir /tmp/bm-tests --output-on-failure
```

`bm_serial_released_compat` requires a directory of fixtures generated against a
released library, supplied as its command-line argument. Building that target
does not verify released-format compatibility. The allocator test also accepts
`block-size` and `ptr-size` to deliberately abort and check failure diagnostics;
these negative checks are not part of its ordinary CTest invocation.

## Make

Existing Make entry points remain available. After moving headers in an existing
checkout, regenerate local dependency files with `make dep` in each affected
driver directory (`stress`, `stress64`, and `perf64`) before rebuilding. Their
Makefiles also track the common headers explicitly.
