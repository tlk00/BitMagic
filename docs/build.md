# Building and configuring BitMagic

BitMagic is a **header-only C++17 library**. To use its C++ containers and algorithms, add the `src` directory to your include path and include the appropriate headers. There is no separate BitMagic binary to build or link. You can vendor the headers, use a Git submodule, or manage the source through your project's dependency system.

The repository's Make and CMake builds compile examples, tests, and utilities. They are useful for exploration and validation, but are not prerequisites for integrating the library. Optional language wrappers and tools have their own build requirements.

## Minimal integration

Keep the headers together: they include one another. For example, save this as `example.cpp`:

```cpp
#include "bm.h"

int main()
{
    bm::bvector<> ids;
    ids.set(10);
    ids.set(100);
    return ids.count() == 2 ? 0 : 1;
}
```

GCC or Clang:

```sh
c++ -std=c++17 -O2 -I/path/to/BitMagic/src example.cpp -o example
./example
```

MSVC, from a developer command prompt:

```bat
cl /std:c++17 /EHsc /O2 /I C:\path\to\BitMagic\src example.cpp
example.exe
```

These commands select no explicit BitMagic SIMD backend. The compiler can still optimize for the target selected by its own flags. Application features such as threading or optional third-party integrations may require additional compiler or linker settings.

## Configure in the build or in source code

BitMagic uses preprocessor definitions for compile-time configuration. You can supply them through your compiler, build system, IDE, or a common application configuration header.

For example, enable wider logical indexing from the command line:

```sh
c++ -std=c++17 -O2 -DBM64ADDR -I/path/to/BitMagic/src example.cpp -o example
```

Or define it before any BitMagic header is included:

```cpp
#define BM64ADDR
#include "bm.h"
```

Alternatively, include `bm64.h` first to select the wider addressing configuration.

Use a consistent configuration across translation units that share BitMagic types or inline implementations. In particular, addressing and SIMD definitions should not vary accidentally between source files. A common configuration header or target-wide build definitions help enforce this. Rebuild affected objects after changing configuration.

## Integrating with your own CMake project

You can describe the header dependency with a small interface target, without adding the repository's examples and tests to your build:

```cmake
cmake_minimum_required(VERSION 3.12)
project(my_app LANGUAGES CXX)

add_library(bitmagic_headers INTERFACE)
target_include_directories(bitmagic_headers INTERFACE
    "${CMAKE_CURRENT_SOURCE_DIR}/third_party/BitMagic/src")
target_compile_features(bitmagic_headers INTERFACE cxx_std_17)

# Optional: propagate a consistent addressing mode to consumers.
# target_compile_definitions(bitmagic_headers INTERFACE BM64ADDR)

add_executable(my_app example.cpp)
target_link_libraries(my_app PRIVATE bitmagic_headers)
```

The interface target propagates include paths and configuration; it does not create a library binary.

## Addressing configuration

| Definition | Effect |
|---|---|
| No `BM64ADDR` | Default 32-bit logical indexing domain |
| `BM64ADDR` | Wider index types and the current 48-bit logical indexing domain |

The end sentinel is reserved: `bm::id_max` is `2^32-1` in the default configuration and `2^48-1` with `BM64ADDR`; valid bit positions are below it. These limits concern logical positions, not the amount of physical RAM allocated. Sparse regions do not require a flat allocation covering the whole domain.

## SIMD and CPU configuration

Select one intended BitMagic SIMD backend and supply the compiler flags it requires. A preprocessor definition selects implementation code; it does not automatically enable the corresponding CPU instructions in your compiler.

| Definition | Backend | Configuration notes |
|---|---|---|
| `BMSSE2OPT` | x86 SSE2 | Use an SSE2-capable target |
| `BMSSE42OPT` | x86 SSE4.2 and associated optimized operations | Use matching target flags; the repository provides a CMake preset value |
| `BMAVX2OPT` | x86 AVX2 and associated optimized operations | Match the full target feature set, including supporting bit-manipulation instructions |
| `BMAVX512OPT` | x86 AVX-512 | Specialized configuration; inspect target requirements and validate on the deployment CPU |
| `BMNEONOPT` | Arm NEON through SSE-to-NEON translation | Uses the bundled `src/sse2neon.h` header |
| `BMWASMSIMDOPT` | WebAssembly SIMD through translated intrinsics | Use an Emscripten build with SIMD and the required SSE compatibility flags |

Do not combine backend definitions. The core C++ library does not perform runtime CPU detection and dispatch among separately compiled SIMD variants. For heterogeneous deployments, choose a compatible baseline or implement separately isolated variants and dispatch in the application.

### x86 with GCC or Clang

For a Skylake-class AVX2 target, a direct build can use:

```sh
c++ -std=c++17 -O2 -march=skylake -mavx2 -DBMAVX2OPT \
    -I/path/to/BitMagic/src example.cpp -o example
```

The resulting executable requires a compatible CPU. Avoid `-march=native` for binaries intended to run on an unknown or older machine. With MSVC, set the matching architecture option in addition to the preprocessor definition; GCC/Clang flags are not interchangeable with MSVC options.

### Arm and Apple Silicon

Start with the default implementation, or select `BMNEONOPT`; the required SSE2NEON header is bundled in `src`. Use the compiler's appropriate Arm target settings. Selecting an Arm CPU does not itself select the BitMagic NEON backend.

### WebAssembly

For an Emscripten C++ build using the translated SIMD backend, the relevant settings include:

```sh
em++ -std=c++17 -O2 -msse4.2 -msimd128 -DBMWASMSIMDOPT \
    -sALLOW_MEMORY_GROWTH=1 -sDISABLE_EXCEPTION_CATCHING=0 \
    -I/path/to/BitMagic/src example.cpp -o example.js
```

Match exception handling and runtime options to your application and Emscripten version. File mapping and other operating-system-specific examples may not apply to a WebAssembly runtime.

### Cross-compilation

Select the compiler, sysroot, and target architecture through your toolchain. For targets without a dedicated backend, including RISC-V, start with the portable implementation. Do not inherit host-specific SIMD definitions or `-march=native` flags.

## Other configuration hooks

| Definition or mechanism | Purpose and scope |
|---|---|
| `BM_NO_STL` | Restricted core configuration for integrations that avoid STL; it does not make every algorithm, stream adapter, or example STL-free |
| `BM_HASRESTRICT`, `BMRESTRICT` | Compiler-specific restrict annotation hooks; inspect `bmdef.h` for compiler defaults before overriding |
| Container allocator template parameters | Customize allocation policy; see `samples/bvsample06` |

These hooks are feature-specific. Validate the headers and operations actually used by the application rather than assuming that a core configuration applies to the whole repository. See [`bmdef.h`](../src/bmdef.h) for implementation details.

## Building repository examples with CMake

Run from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target bvsample27 --parallel 4
```

The repository places example executables under `build/bin` (multi-configuration generators can add a configuration directory). Build another named target, or omit `--target` to build the default set. The project fixes its executable output path under the source tree's `build/bin`, even when using another CMake build directory.

For the x86 backends explicitly handled by the current repository CMake configuration:

```sh
cmake -S . -B build-avx2 -DCMAKE_BUILD_TYPE=Release -DBMOPTFLAGS=BMAVX2OPT
cmake --build build-avx2 --target bvsample27 --parallel 4
```

`BMSSE42OPT` is the other explicitly handled `BMOPTFLAGS` value. Do not assume this option is a general dispatcher for every macro in the table above. Other backends need their definitions, compiler settings, and dependencies configured explicitly.

On native x86 GCC/Clang builds, the repository's default configuration currently adds `-march=native`. Thus, an uncustomized repository build is not necessarily suitable for redistribution to older CPUs. For controlled deployment targets, use explicit settings in your application build.

### Xcode and Visual Studio

Generate IDE projects with CMake:

```sh
cmake -S . -B build-xcode -G Xcode
cmake --build build-xcode --config Release --target bvsample27
```

On Windows, run CMake with the installed Visual Studio generator and build with `--config Release`. You can also integrate headers directly into an existing IDE project: configure C++17, the include directory, and consistent preprocessor definitions. There is no requirement to generate the BitMagic repository project.

## Building with GNU Make

From the repository root, build a selected example:

```sh
make -C samples/bvsample27 rebuild
make -C samples/bvsample27 DEBUG=YES rebuild
```

The shared `makefile.in` derives `PROJECT_DIR` from its own location when unset. Sourcing `bmenv.sh` is available for existing workflows but is not required by that root-discovery mechanism.

On a supported x86 platform, the traditional Make interface uses a compiler definition as its option value:

```sh
make -C samples/bvsample27 BMOPTFLAGS=-DBMAVX2OPT rebuild
```

This differs from CMake's `-DBMOPTFLAGS=BMAVX2OPT` spelling. Architecture-flag selection is implemented in individual files under `platforms/`; inspect the selected platform file when adding a target or backend. To build the root collection, use `make rebuild`.

## Validation and further reading

First compile and run a small example in the actual deployment configuration. Then validate the APIs and data distributions your application depends on. Changing a SIMD target, compiler, or addressing mode warrants validation in that configuration.

- [Sample catalogue](../samples/readme.md)
- [Streaming serialization and gather example](../samples/bvsample27/readme.md)
- [Repository CMake configuration](../CMakeLists.txt)
- [Shared Make configuration](../makefile.in)
- [Design and methodology](https://bitmagic.io/design)
- [Technical articles](https://bitmagic.io/articles)

Examples establish integration and demonstrate API usage; they do not replace the full randomized stress suite.
