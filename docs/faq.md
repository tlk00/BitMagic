# BitMagic FAQ

[Project overview](../README.md) · [Build and configuration](build.md) · [Examples](https://github.com/tlk00/BitMagic/tree/master/samples)

This FAQ combines questions raised in public discussions with introductory questions about current BitMagic capabilities. Historical discussions identify user concerns; they do not establish current performance or whether an old defect still exists.

## Contents

1. [Core concepts and container choices](#1-core-concepts-and-container-choices)
   - [What is a succinct data structure?](#what-is-a-succinct-data-structure)
   - [What does “data as an index” mean?](#what-does-data-as-an-index-mean)
   - [Are sparse vectors useful for dense data?](#are-sparse-vectors-useful-for-dense-data)
   - [Does BitMagic support arbitrary sets and unknown logical values?](#does-bitmagic-support-arbitrary-sets-and-unknown-logical-values)
   - [How does BitMagic compare with standard C++ bit containers?](#how-does-bitmagic-compare-with-standard-c-bit-containers)
2. [Search and performance](#2-search-and-performance)
   - [How can I implement Boolean document or record retrieval?](#how-can-i-implement-boolean-document-or-record-retrieval)
   - [How does BitMagic implement rank and select, and are they fast?](#how-does-bitmagic-implement-rank-and-select-and-are-they-fast)
   - [How should I evaluate counting, iteration, and set-operation performance?](#how-should-i-evaluate-counting-iteration-and-set-operation-performance)
3. [Compression and persistence](#3-compression-and-persistence)
   - [What is the difference between compact memory and serialized compression?](#what-is-the-difference-between-compact-memory-and-serialized-compression)
   - [Can I store BitMagic data in a database or a hybrid storage system?](#can-i-store-bitmagic-data-in-a-database-or-a-hybrid-storage-system)
   - [Can I retrieve a few values without restoring the whole container?](#can-i-retrieve-a-few-values-without-restoring-the-whole-container)
   - [Is a deserialization index the same as a search index?](#is-a-deserialization-index-the-same-as-a-search-index)
   - [Does serialization require a complete temporary RAM BLOB, or can it stream to disk?](#does-serialization-require-a-complete-temporary-ram-blob-or-can-it-stream-to-disk)
4. [Integration and configuration](#4-integration-and-configuration)
   - [How do I start using BitMagic as a library?](#how-do-i-start-using-bitmagic-as-a-library)
   - [Why is BitMagic 32-bit by default?](#why-is-bitmagic-32-bit-by-default)
   - [Can I use BitMagic from Python, Rust, or other languages?](#can-i-use-bitmagic-from-python-rust-or-other-languages)
5. [Memory control and detailed API behavior](#5-memory-control-and-detailed-api-behavior)
   - [What is the best runtime allocator for BitMagic?](#what-is-the-best-runtime-allocator-for-bitmagic)
   - [Can I prevent all allocations after initialization?](#can-i-prevent-all-allocations-after-initialization)
   - [Why does flipping a bit-vector produce so many set bits?](#why-does-flipping-a-bit-vector-produce-so-many-set-bits)
6. [Further questions and reporting problems](#6-further-questions-and-reporting-problems)

## 1. Core concepts and container choices

### What is a succinct data structure?

A succinct data structure stores information compactly while retaining efficient operations on that representation. Unlike an opaque compressed file that must first be fully decoded, it provides ways to access or query the represented data directly.

In theoretical computer science, “succinct” has a stronger meaning: space approaches the information-theoretic minimum with lower-order redundancy. BitMagic uses the term to describe its practical compact, operational representations; it is not a universal guarantee that every container and dataset meets that formal space bound.

One example is an RSC vector: a validity bit-vector records which logical positions are assigned, while a compact sequence holds their values. Rank maps an assigned logical position to the compact sequence; select maps back. Bit-transposed columns also support searches on their planes without reconstructing all scalar values. See the [design overview](https://bitmagic.io/design).

Sample references: [svsample01 — bit-transposed integer vectors](https://github.com/tlk00/BitMagic/tree/master/samples/svsample01), [rscsample01 — rank-select compressed values](https://github.com/tlk00/BitMagic/tree/master/samples/rscsample01).

### What does “data as an index” mean?

Succinct data structures retain efficient access and query operations while storing information compactly. BitMagic's bit-transposed columns apply that principle to search: their representation supports logical evaluation of predicates and progressive pruning of candidate positions, without reconstructing every scalar value.

For example, an equality search intersects candidates with the planes for required one-bits and excludes candidates using the planes for required zero-bits. When a candidate block becomes empty, subsequent conditions do not need to be evaluated there. Optimized kernels exploit this pruning to skip eliminated regions and use bitwise and SIMD operations on surviving candidates. Thus, the representation can reduce both memory traffic and the work remaining at later search stages. The gain depends on the data and predicate; retaining most candidates limits pruning.

The column can therefore be both data and a searchable structure, without a separate secondary search index for those predicates. Applications can also build explicit posting-list indexes or combine indexed retrieval with column searches.

This does not imply constant-time search or the absence of auxiliary structures. See [scanner methodology](https://bitmagic.io/sparse-vector-search) and [integer search examples](https://github.com/tlk00/BitMagic/tree/master/samples/svsample10).

### Are sparse vectors useful for dense data?

Yes. Their role is to store searchable columns of integer, string, or floating-point values in compact representations. The container name does not require most positions to be NULL. Dense values can benefit from limited bit width, repeated patterns, small character alphabets, or regularity within bit-planes. RSC additionally compacts unassigned positions. NULL is distinct from a stored numerical zero. High-entropy data may offer limited savings, so measure the intended dataset.

See the [container overview](../README.md#containers-for-sparse-and-dense-data) and [compression methodology](https://bitmagic.io/design).

Sample references: [svsample01 — array import and compact integer storage](https://github.com/tlk00/BitMagic/tree/master/samples/svsample01), [svfsample01 — floating-point vectors](https://github.com/tlk00/BitMagic/tree/master/samples/svfsample01), [bvsample03 — dense bit-vectors and compression strategies](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample03).

### Does BitMagic support arbitrary sets and unknown logical values?

Bit-vectors represent arbitrary sets of integer IDs or positions, not just search results. Boolean set algebra is available, and [`bm3vl.h`](../src/bm3vl.h) supplies Kleene three-valued logic using True, False, and Unknown. See the [three-valued logic example](https://github.com/tlk00/BitMagic/tree/master/samples/bv3vlogic).

### How does BitMagic compare with standard C++ bit containers?

The standard library offers `std::bitset<N>`, whose size is fixed at compile time and which provides bitwise operations and counting, and `std::vector<bool>`, a dynamically sized specialization with a different interface. Neither specifies BitMagic-style adaptive sparse/dense block compression or a compressed binary persistence format. See the standard descriptions of [`bitset`](https://eel.is/c++draft/template.bitset) and [`vector<bool>`](https://eel.is/c++draft/vector.bool).

For short bit sequences and straightforward workloads, a flat representation can use less memory and execute faster than BitMagic's hierarchy of blocks. Allocation, metadata, and indirection have costs. Standard containers can also be appropriate for large dense arrays when their operations match the task; their usefulness is not limited to small inputs.

BitMagic adds compressed set representations, enumeration, aggregation, rank/select, binary metrics, and compressed serialization with selective retrieval. Its advantages become relevant when these capabilities or the dataset's structure justify the additional machinery. Compare both memory footprint and operation throughput on your workload rather than assuming one representation always wins.

Sample references: [bvsample03 — memory statistics and block strategies](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample03), [bvsample08 — STL interoperability](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample08). These illustrate capabilities and integration, rather than a direct performance comparison with standard containers.

---

## 2. Search and performance

### How can I implement Boolean document or record retrieval?

Assign each document or record an integer ID. Represent the IDs satisfying each term or condition as a bit-vector, then combine sets using AND, OR, XOR, or MINUS. Enumerate matching IDs, count them, or use them to select associated column values. Negation needs an explicitly understood universe of valid IDs.

This pattern appears in public questions about [full-text retrieval](https://stackoverflow.com/questions/64286076/fast-search-algorithm) and [complex Boolean conditions over large datasets](https://stackoverflow.com/questions/79644546/bit-manipulation-for-big-data-processing-with-complex-boolean-logic). Start with [set algebra](https://github.com/tlk00/BitMagic/tree/master/samples/bvsetalgebra) and [multi-vector aggregation](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample16).

### How does BitMagic implement rank and select, and are they fast?

`rank(i, index)` counts set bits in the inclusive interval `[0, i]`. `select(k, position, index)` finds the position of the **k-th set bit**, using a one-based rank. For set bits at positions 2, 4, and 5, `rank(4)` is 2 and `select(2)` returns position 4.

BitMagic provides a reusable `rs_index_type`, constructed with `build_rs_index()`. It summarizes population counts so queries can narrow the work to relevant blocks and subdivisions, then use representation-specific counting and selection kernels. Depending on the build and CPU, these paths use hardware bit operations and SIMD rather than rescanning the entire vector.

This is designed for fast repeated queries, but index construction, index memory, cache locality, and the data distribution remain part of the cost. Rebuild an externally maintained rank/select index after modifying its source vector. Use the [rank/select example](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample17) and [implementation article](https://bitmagic.io/rank-select) to explore the method; historical timings in the article are not current performance guarantees.

### How should I evaluate counting, iteration, and set-operation performance?

Measure the operations the application actually needs. If only a cardinality is required, consider count-oriented APIs rather than materializing a result and then counting it. For enumeration, use the bit-vector enumerator or traversal algorithms. For combinations of many sets, examine the aggregator.

Compare memory footprint and throughput on representative densities, runs, ID ranges, and query patterns, using equivalent compiler and CPU settings. Public discussions about [large integer sets](https://stackoverflow.com/questions/14981143/fast-implementation-of-operations-on-large-sets-of-quite-big-integers) and [fast bit operations](https://stackoverflow.com/questions/14432050/high-performance-library-for-bit-wise-operations) illustrate these concerns; their historical measurements are not current benchmarks.

See [counting](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample11), [traversal](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample25), and [SIMD configuration](build.md#simd-and-cpu-configuration).

---

## 3. Compression and persistence

### What is the difference between compact memory and serialized compression?

Live containers support access and supported updates, searches, and logical operations in a compact representation. Serialization adds encoding for storage and transfer. Deserialization reconstructs BitMagic containers; it need not expand the entire dataset into ordinary arrays. See the [two-layer compression explanation](../README.md#compression-for-computation-and-storage).

Sample references: [bvsample03 — in-memory optimization](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample03), [svsample02 — serialization and XOR compression](https://github.com/tlk00/BitMagic/tree/master/samples/svsample02), [rscsample05 — serialization of related columns](https://github.com/tlk00/BitMagic/tree/master/samples/rscsample05).

### Can I store BitMagic data in a database or a hybrid storage system?

Applications can store serialized BLOBs in files or database binary fields and combine database-managed metadata with file-based data. The application supplies database calls, transactions, and storage policy. The stream interfaces and serialization hooks can support custom integrations subject to their documented contracts.

Storing a BLOB in a database does not automatically provide remote selective I/O. The provided C++ stream adapter requires seekable binary streams, including output patching. See [`bmfio.h`](../src/bmfio.h) and the [persistence overview](../README.md#persistence-and-integration-architectures).

Sample references: [bvsample04 — serialized BLOB construction](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample04), [bvsample27 — file I/O](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample27), [strsvsample10 — memory-mapped retrieval](https://github.com/tlk00/BitMagic/tree/master/samples/strsvsample10). These demonstrate storage building blocks; database client integration remains application-specific.

### Can I retrieve a few values without restoring the whole container?

Range and gather deserialization restore selected regions or logical positions. Bookmarks and optional deserialization indexes help skip unrelated encoded blocks. Memory-mapped files and seekable streams provide useful storage access paths. Retrieval remains block-oriented and still performs decoding; benefits depend on selection, layout, and storage behavior.

See [file-based bit-vector retrieval](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample27) and [memory-mapped strings](https://github.com/tlk00/BitMagic/tree/master/samples/strsvsample10).

### Is a deserialization index the same as a search index?

No. It navigates one particular serialized BLOB using marker and bookmark offsets. The application supplies positions to retrieve. A scanner instead evaluates predicates on a live compact container. Navigation indexes can be persisted to avoid rebuilding them, but must be rebuilt when the serialized representation changes.

See [float search followed by indexed retrieval](https://github.com/tlk00/BitMagic/tree/master/samples/svfsample05).

### Does serialization require a complete temporary RAM BLOB, or can it stream to disk?

**BitMagic 9.3.1 supports both approaches.** The RAM-buffer API produces a contiguous serialized BLOB, useful when an application wants to keep compressed data in RAM or pass one binary object to a database or another storage interface. It remains a useful option, but is not a requirement of all serialization paths.

For direct file output, [`bmfio.h`](../src/bmfio.h) provides `bm::streams_encoder`. Serializers write through reusable working buffers without constructing a complete temporary serialized BLOB. Streaming support covers bit-vectors and integer, RSC, string, and floating-point sparse vectors. With equivalent settings, stream and RAM serialization produce the same encoded representation. Stream-based deserializers also avoid loading the entire input BLOB into a temporary RAM buffer.

The supplied adapter requires **seekable binary streams**: serialization may patch headers and bookmarks. Do not use append mode or assume it is a forward-only socket writer. The owning caller completes output with `finish()` and checks errors. Streaming still needs working memory and the source container; it removes the requirement for a second, full serialized copy in RAM.

See the [file I/O example](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample27) and [streaming string-vector example](https://github.com/tlk00/BitMagic/tree/master/samples/strsvsample10). When consulting older examples or comparisons, check the BitMagic version and which API they use.

---

## 4. Integration and configuration

### How do I start using BitMagic as a library?

The C++ library is header-only. Add `src` to the include path, select C++17, and include the required headers. The repository build compiles examples, utilities, and tests; it is not a prerequisite for using the headers. Optional language wrappers have separate requirements. Use a release branch for stable application integration, not development `master`.

See the [build guide](build.md). Public requests include [CMake installation support (#49)](https://github.com/tlk00/BitMagic/issues/49) and [Conan packaging (#54)](https://github.com/tlk00/BitMagic/issues/54). These requests should not be interpreted as confirmation that a particular package-manager integration is available.

For basic use, `#include "bm.h"` provides `bm::bvector<>`; other containers and algorithms have their own headers. Configure macros before any BitMagic header, either through compiler definitions or a common configuration header. For example, `-DBM64ADDR` and `#define BM64ADDR` select the same addressing option.

| Definition | Role |
|---|---|
| `BM64ADDR` | Enable the current 48-bit logical indexing domain |
| `BMSSE2OPT` | Select the x86 SSE2 backend |
| `BMSSE42OPT` | Select the x86 SSE4.2 backend, including hardware POPCNT for population counting, which is critical for performance |
| `BMAVX2OPT` | Select the x86 AVX2 backend with POPCNT and BMI instructions; highly recommended for rank/select acceleration on supported CPUs |
| `BMAVX512OPT` | Select the experimental x86 AVX-512 backend; use at your own risk |
| `BMNEONOPT` | Select Arm NEON through the bundled SSE2NEON translation |
| `BMWASMSIMDOPT` | Select WebAssembly SIMD through translated intrinsics |
| `BM_NO_STL` | Restricted core configuration; not a promise that every auxiliary API is STL-free |
| `BM_HASRESTRICT`, `BMRESTRICT` | Compiler-specific restrict annotation hooks |

Choose one SIMD backend, match its compiler target flags to the deployment CPU, and keep definitions consistent across translation units. The core C++ library does not automatically dispatch among CPU variants. See the [detailed configuration guide](build.md) for commands, dependencies, and limitations.

Sample references: [bvsample01 — bit-vector basics](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample01), [svsample01 — integer-vector basics](https://github.com/tlk00/BitMagic/tree/master/samples/svsample01).

### Why is BitMagic 32-bit by default?

Here “32-bit” describes the logical indexing domain, not the CPU or process pointer width. A default BitMagic build can run in a 64-bit application. Many indexes, partitions, and application datasets fit within the default domain, so they do not require wider logical indices or the larger hierarchy used for wider addressing.

Keeping this configuration as the default can reduce indexing and metadata costs, particularly for edge computing, embedded applications, and WebAssembly. This is not a claim that every allocation is halved: memory consumption depends on the container, occupied blocks, and data distribution.

Define **`BM64ADDR`** before including BitMagic headers, or include **`bm64.h`** first, to use wider index types and the current **48-bit logical domain**. The configuration is not a full 64-bit universe. `bm::id_max` is the reserved end sentinel: `2^32-1` by default and `2^48-1` with `BM64ADDR`; valid bit positions are below it. See the [addressing guide](build.md#addressing-configuration).

Sample reference: [bvsample01_64 — wider addressing and high-position bits](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample01_64).

### Can I use BitMagic from Python, Rust, or other languages?

BitMagic provides a [C interface, `libbm`](../lang-maps/libbm), as a foundation for language bindings. Its [public C API](../lang-maps/libbm/include/libbm.h) currently exposes `bm::bvector<>` functionality. The mapping is incomplete: integer, RSC, string, and floating-point sparse-vector containers are not yet covered. Extending this coverage needs additional contributors. The repository also contains [experimental Java JNI mappings](../lang-maps/readme).

The core supports **`BM_NO_STL`**, which helps make integration through C practical. The C wrapper combines this approach with C-style allocation and an implementation configured without C++ exceptions or RTTI. This allows builds without a C++ runtime dependency, as described in the [language-mapping documentation](../lang-maps/readme). `BM_NO_STL` alone is not a universal runtime-removal switch: the allocator, error handling, compiler, and linker configuration also matter. The wrapper is still compiled from C++ source and presents a C interface to callers.

Community bindings include the historical [Python `bitmagic` package on PyPI](https://pypi.org/project/bitmagic/) and the [Rust `bitmagic` crate](https://docs.rs/bitmagic/latest/bitmagic/). These are independent community efforts, with their own version compatibility, coverage, and maintenance. The Python package describes a Boost.Python-based wrapper; the Rust documentation identifies a release based on BitMagic 7.7.7. Treat them as starting points to evaluate, rather than assuming coverage of current BitMagic APIs. See also the [Python-interface discussion (#55)](https://github.com/tlk00/BitMagic/issues/55).

**Contributors are actively welcome.** Useful work includes extending the C API to other containers, updating Python and Rust bindings, packaging, cross-platform builds, tests, documentation, and examples. Limited contributor capacity, rather than lack of interest in these integrations, has constrained this work. The [C wrapper implementation](../lang-maps/libbm/src/libbm.cpp) and [mapping build configuration](../lang-maps/CMakeLists.txt) provide starting points.

Code reference: [libbmtest.c — C API test calls](https://github.com/tlk00/BitMagic/blob/master/lang-maps/test/libbmtest.c). This is test code for the C wrapper, not a complete Python or Rust binding tutorial.

---

## 5. Memory control and detailed API behavior

### What is the best runtime allocator for BitMagic?

BitMagic is compatible with standard runtime allocators and has also been tested successfully with third-party allocators. Its containers support C++-style allocator customization, allowing an application to integrate its own allocation policy. Custom allocators must preserve the alignment required by the selected SIMD configuration. See [bvsample06 — custom allocator integration](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample06).

Among third-party allocators, we recommend **[jemalloc](https://jemalloc.net/)** based on successful experience with BitMagic workloads. BitMagic commonly allocates fixed-size 8 KiB bitmap blocks and compressed GAP blocks from a small family of sizes. This predictable pattern works well with jemalloc's size classes and block reuse, helping control fragmentation. Its decay and purging mechanisms can also return unused physical memory to the operating system as BitMagic containers release blocks. The timing and extent of reclamation still depend on the platform, allocator configuration, and workload.

### Can I prevent all allocations after initialization?

A logical size limit does not preallocate all storage. BitMagic allocates blocks as needed, and mutations or algorithms can require additional blocks or temporary memory. Custom allocation policies can help control memory management, but a strict no-allocation phase requires auditing the exact operations and providing appropriate storage and scratch-memory policies. Do not treat construction with a maximum bit count as a no-allocation guarantee.

See the [custom allocator example](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample06) and the original [allocation-control question (#82)](https://github.com/tlk00/BitMagic/issues/82).

### Why does flipping a bit-vector produce so many set bits?

Complement is relative to the vector's logical universe, not its last set bit. A default `bvector<>` already has logical size `bm::id_max`: `2^32-1` in the default configuration, with valid positions from `0` through `2^32-2`. The final value is reserved as an end sentinel. With `BM64ADDR`, the corresponding logical size is `2^48-1`.

Physical storage is allocated on demand through a hierarchy of blocks. Blocks can be allocated as needed in arbitrary order, so bits can be set in ascending, descending, or random order. Setting a high-position bit does not require allocating every preceding data block or filling earlier positions. This supports sets whose IDs arrive out of order, including workloads that set bits backwards. Access order can still affect cache locality and throughput.

The default is therefore a full logical universe with dynamically allocated storage, rather than a short vector whose logical size follows its highest set bit. Setting a few bits does not shrink that universe: `flip()` complements every position within the logical size and can turn a large implicit zero tail into ones.

For a complement over a known contiguous domain, set the appropriate logical size before flipping. For an application-specific domain, form `universe MINUS selected` using a bit-vector containing valid IDs. This also handles noncontiguous domains and avoids selecting IDs that do not represent real objects.

The distinction arose in [the `flip()` report (#73)](https://github.com/tlk00/BitMagic/issues/73). The explanation here follows the current constructor and inversion implementation in [`bm.h`](../src/bm.h); it is not a claim about the resolution of that issue.

Sample references: [bvsample02 — set difference for universe-minus-selection](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample02), [bvsample01_64 — full vectors and address boundaries](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample01_64).

---

## 6. Further questions and reporting problems

For a reproducible report, include the release or commit, compiler and platform, addressing/SIMD definitions, a minimal example, and expected versus observed behavior. Performance reports should also describe the data distribution, measured operation, and whether allocation and serialization are included in timing.

Consult [GitHub issues](https://github.com/tlk00/BitMagic/issues), the [website's technical articles](https://bitmagic.io/articles), and [use cases](https://bitmagic.io/use-case). For the project's attribution policy, see [License and contributions](../README.md#license-and-contributions).
