# BitMagic C++ Library

**Compact data. Searchable representations. Hardware-efficient computation.**

BitMagic is a header-only C++17 library for compressed bit-vectors and succinct integer, string, and floating-point vectors. It provides building blocks for information retrieval, data management, data science, and scientific applications that need to search and process large datasets within a practical memory and bandwidth budget.

BitMagic treats data representation as part of the algorithm. Its containers combine memory efficiency with search-friendly layouts, Boolean operations, rank/select support, and data-parallel processing. Applications can use BitMagic to build indexes, search compact columns without separate secondary indexes, or combine both approaches.

The library supports a complete data lifecycle: construct compact containers, search and combine selections, gather values for computation, serialize to storage, and restore selected data when needed.

- **Searchable data:** evaluate supported predicates on bit-transposed columns using scanners and logical operations.
- **Compressed sets and indexes:** represent posting lists, row selections, membership information, and binary features with `bm::bvector<>`.
- **Sparse and dense vectors:** exploit value width, NULL distribution, runs, and other regularities through adaptive representations.
- **Hardware-aware algorithms:** use SIMD, block-oriented processing, and specialized kernels to reduce memory traffic and improve throughput.
- **Flexible persistence:** use RAM BLOBs, seekable file streams, memory-mapped files, or application-managed database storage.

[Website](https://bitmagic.io) · [FAQ](docs/faq.md) · [Examples](https://github.com/tlk00/BitMagic/tree/master/samples) · [API documentation](https://bitmagic.io/doxygen/html/modules.html) · [Technical articles](https://bitmagic.io/articles.html) · [Releases](https://github.com/tlk00/BitMagic/releases)

## Data as an index

Succinct data structures combine compact storage with efficient access and query operations on that representation. BitMagic applies this principle to search: the way values are represented both saves memory and exposes opportunities to eliminate candidates without reconstructing each value. This is the connection between succinct storage and **data as an index**.

Integer and string vectors use a bit-transposed layout. Instead of storing every value as an independent machine word or character sequence, the representation organizes bits into planes backed by compressed bit-vectors. Supported searches can evaluate these planes with logical operations and produce a bit-vector of matching positions, without first expanding the entire column into ordinary scalar values.

This layout enables aggressive **search-space pruning**. For example, an equality search can intersect candidates with planes corresponding to required one-bits and exclude candidates using planes corresponding to required zero-bits. As these logical operations eliminate candidates, a candidate block can become empty. Once no candidates remain in that block, subsequent conditions need no further evaluation there. Optimized search kernels exploit this reduction to skip work on eliminated regions, while processing surviving candidates in groups with bitwise and SIMD operations.

Search efficiency therefore comes from both compactness and the structure of the computation: fewer bytes need to move, and progressive pruning can reduce the amount of data that later stages inspect. The benefit depends on the predicate, value distribution, and how quickly candidates are eliminated; queries that retain most candidates offer less opportunity to skip work.

This makes **data as an index** a useful way to design with BitMagic. A compact column can serve both as the stored data and as the searchable structure. `bm::sparse_vector_scanner<>` supplies search operations for supported container types, including equality and comparison operations, ranges, and string searches. The available predicates depend on the container and value type.

Applications can choose among complementary approaches:

- **Build explicit indexes:** store document IDs or row IDs in bit-vectors and combine posting lists or filter sets.
- **Search the data representation:** use scanners on succinct columns without maintaining a separate secondary search index for those predicates.
- **Combine both:** use an index to identify a population, combine it with column-search results, and gather the selected values.

“Without a separate search index” does not mean constant-time access or the absence of all auxiliary structures. Scanners still perform work over the representation; rank/select indexes and other accelerators can support particular operations. The advantage is that searching can operate directly on the compact data layout.

For the method and experimental background, see [searching bit-transposed integer vectors](https://bitmagic.io/sparse-vector-search) and [searchable scientific dictionaries](https://bitmagic.io/star-search). The latter illustrates combining binary search with logical search over compact string data.

## Performance beyond operation counts

BitMagic's performance comes from the interaction of algorithms, representation, and implementation. Asymptotic complexity matters, but it does not describe cache misses, bytes transferred, branch behavior, allocation overhead, or how much useful work each CPU instruction performs.

A scan over compact bit-planes can be competitive with a theoretically more selective algorithm that performs scattered memory accesses. Which approach wins depends on data size, selectivity, representation, and hardware. BitMagic is designed to exploit the cases where compact storage and regular, data-parallel work reduce the real cost of computation.

Important techniques include:

- **SIMD kernels:** process multiple bits or values per instruction using supported x86, Arm NEON, and WebAssembly SIMD paths.
- **Block-oriented representations:** adapt processing to empty, full, run-compressed, and bitmap blocks.
- **Memory locality:** organize work around reusable blocks and buffers, and offer immutable layouts that reduce fragmentation.
- **Bandwidth efficiency:** avoid moving expanded representations when operations can use compact data.
- **Combined operations:** use aggregators and count-oriented operations to reduce intermediate materialization in supported workflows.
- **Batched extraction and traversal:** decode or gather values in groups and reuse scratch memory, including with `bm::for_each_sparse()`.

These techniques complement algorithmic complexity analysis. Performance and compression should be measured with representative data, query distributions, and target hardware; high-entropy data and dense extraction have different tradeoffs from highly selective queries over structured columns.

The [scanner design and benchmark article](https://bitmagic.io/sparse-vector-search) explains how logical operations, early elimination of candidate blocks, cache blocking, and SIMD work together. The [floating-point study](https://bitmagic.io/bm-svf) adds measurements for financial, linear, and random datasets. These studies describe particular versions, datasets, and machines; their results are evidence for the methods, not universal performance guarantees.

## Containers for sparse and dense data

| Container | Typical role |
|---|---|
| `bm::bvector<>` | Compressed integer sets, posting lists, selection masks, and binary features |
| `bm::sparse_vector<>` | Bit-transposed integer columns with optional NULL support |
| `bm::rsc_sparse_vector<>` | Rank-select compressed columns with many unassigned positions |
| `bm::str_sparse_vector<>` | Compact string collections with search and optional alphabet remapping |
| `bm::sparse_vector_float<>` | Floating-point columns with search, extraction, and serialization |

The name “sparse vector” does not restrict these containers to mostly empty data. A fully populated column can benefit from limited value width, repeated patterns, or regularity within its bit-planes. Actual memory savings depend on the data and include representation overhead.

**NULL and zero are distinct.** Rank-select compressed vectors use a NOT NULL bit-vector to map logical positions to stored assigned values. Related columns can share NULL information, and supported RSC arrangements can also share rank/select indexing. This is useful for groups of columns describing the same observations or entities.

The [compression design overview](https://bitmagic.io/design) explains the relationships among block compression, bit-transposition, rank/select, and string remapping. For floating-point representation and search, see [Using Sparse Vector Floats on Financial Datasets](https://bitmagic.io/bm-svf).

## Information retrieval and set processing

BitMagic's compressed bit-vectors represent arbitrary sets of integer identifiers or positions and support Boolean set algebra. Applications can use them for membership, relationships, selections, and logical inference, including Boolean retrieval:

- AND, OR, XOR, MINUS, and NOT operations.
- Multi-vector aggregation, including combined AND-MINUS operations.
- Cardinality, intersection counts, rank/select, and enumeration.
- Binary similarity and distance calculations.
- Bulk construction, interval traversal, and partitioned processing patterns.

BitMagic also supports **Kleene three-valued logic** through [`bm3vl.h`](https://github.com/tlk00/BitMagic/blob/master/src/bm3vl.h). A compact two-bit-vector representation models True, False, and Unknown, with logical NOT, AND, and OR operations. This allows applications to express logical conditions where information may be missing or unknown. See the [three-valued logic example](https://github.com/tlk00/BitMagic/tree/master/samples/bv3vlogic).

Applications can use these operations for posting lists, document filters, exclusions, cohort selection, and binary feature comparison. Search results remain bit-vectors, making subsequent filtering and combination natural.

Start with [set algebra](https://github.com/tlk00/BitMagic/tree/master/samples/bvsetalgebra), [aggregation](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample16), and [rank/select](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample17).

## Data science, data management, and scientific computing

BitMagic provides components for compact columnar models and selective processing. A typical workflow searches columns, combines row selections, and gathers only the values needed by the next computation.

Examples include:

- **Data management:** nullable columns, shared validity information, Boolean row filters, and selective materialization.
- **Data science:** cohort construction, binary feature comparison, threshold searches, and batch processing of selected observations.
- **Scientific computing:** genomic intervals and variants, categorical sequences, observation catalogues, and compact associations.
- **Memory-constrained applications:** compressed working sets for desktop applications, embedded systems, and WebAssembly.

Extracted values can feed application-specific numerical or statistical routines. BitMagic supplies the compact representation, search, selection, and data movement components around those computations.

Explore [integer comparisons](https://github.com/tlk00/BitMagic/tree/master/samples/svsample10), [shared NULL planes](https://github.com/tlk00/BitMagic/tree/master/samples/rscsample07), [float selection and retrieval](https://github.com/tlk00/BitMagic/tree/master/samples/svfsample05), and [genomic interval representation](https://github.com/tlk00/BitMagic/tree/master/samples/xsample08).

The website's [use-case collection](https://bitmagic.io/use-case) provides application studies covering histograms, genomic intervals, variant search, searchable astronomical dictionaries, and scheduling. The [scientific dictionary study](https://bitmagic.io/star-search) and [floating-point data study](https://bitmagic.io/bm-svf) connect these applications to their representation and search methods.

## Compression for computation and storage

In production information retrieval systems, compression is an architectural decision. It affects how much of an index fits in memory, how much data a query moves, and how efficiently the system uses CPU caches, storage, and network bandwidth.

BitMagic combines complementary compression methods in two layers. The first keeps data compact and operational in memory; the second encodes it more deeply for storage and transfer. Both layers share the same container model, allowing applications to choose a balance among memory footprint, storage size, and processing cost.

| Layer | Techniques | Purpose |
|---|---|---|
| **Operational representation** | Bit-slicing, hierarchical block compression, delta-GAP run-length encoding, character remapping, and rank-select compression | Reduce working memory while retaining access, search, and supported logical operations |
| **Serialized representation** | Adaptive block encoding, Elias-gamma coding, tuned Binary Interpolative Coding, and optional XOR filtering | Reduce persisted size while supporting efficient reconstruction and selective retrieval |

The techniques cooperate rather than act as interchangeable whole-file codecs. Bit-slicing exposes regularities within individual planes. Hierarchical compression represents empty and full regions efficiently. Delta-GAP encoding captures runs, while remapping can reduce the number of active planes. Rank-select compression removes unassigned positions from the stored value sequence while preserving their logical coordinates.

For persistence, serialization selects encoded representations for blocks. Elias-gamma and Binary Interpolative Coding encode integer sequences; XOR filtering can expose similarities between blocks or planes before subsequent encoding. Their effectiveness depends on the distribution and correlations in the data.

The [compression design overview](https://bitmagic.io/design) explains how these methods fit together. The [XOR compression article](https://bitmagic.io/bm-xor) develops the method for correlated bit-transposed vectors, using aligned biological sequences as an example.

### Compression that preserves selective access

Bookmarks and optional deserialization indexes provide navigation through the serialized representation. They let range and gather deserialization locate relevant blocks and skip unrelated payloads, so applications can retrieve portions of a compressed dataset without restoring the complete container.

Selective access is a property of the serialization architecture, rather than of each codec in isolation. Retrieval remains block-oriented, and encoding dependencies can require additional decoding. A deserialization index identifies navigation points; it does not map every individual value directly to independent compressed bytes.

### Choosing the compression balance

The API exposes choices such as in-memory optimization, serialization compression levels, XOR filtering, bookmark spacing, and construction or persistence of deserialization indexes. These let applications tune several distinct costs:

- Operational memory footprint.
- Serialized size.
- Serialization and deserialization time.
- Memory and I/O traffic during selective retrieval.
- Space and preparation costs of navigation metadata.

More compact output can improve retrieval by reducing memory traffic or I/O, but may also require more encoding or decoding work. Denser bookmarks improve positioning precision at a storage cost. The appropriate configuration depends on whether the workload favors frequent updates, repeated searches, bulk transfer, sequential restoration, or sparse retrieval.

BitMagic offers a coordinated set of controls within a common API. Compression depth is not a single setting with a fixed size-versus-speed tradeoff: the effects depend on the data and workload.

## Persistence and integration architectures

BitMagic separates container representation and serialization from the application's choice of storage system. Serialized containers are binary objects that an application can place in ordinary files, RDBMS BLOB columns, key-value stores, document databases with binary-value support, or object storage.

This permits several integration architectures:

- **File-based persistence:** serialize containers to files and restore them through binary streams or memory mapping.
- **Database-managed persistence:** store serialized values in a relational or post-relational database, using that system for keys, transactions, metadata, and lifecycle management.
- **Hybrid storage:** combine database-managed metadata and identifiers with compressed data files, mapped working sets, or application-managed caches.
- **Custom adapters:** build storage integration around the serialization interfaces and channel hooks, using the open-source implementation as a reference.

Database integration is at the serialized-data boundary: applications supply the database client calls and storage policy. Selective physical I/O depends on the backend's access capabilities; storing a BLOB in a database does not by itself provide block-level remote access.

### Live containers and serialized data

Live containers support compact in-memory access, updates, search, and supported logical operations. Serialization applies additional encoding for storage or transfer. Range and gather deserialization reconstruct selected data into BitMagic containers, avoiding a mandatory expansion of the whole dataset into ordinary arrays.

Search scanners and deserialization indexes have different roles. A scanner evaluates predicates on live compact columns. A deserialization index records navigation information for a particular BLOB and accelerates retrieval of positions already selected by the application.

### Streaming and selective retrieval

The new [`bmfio.h`](https://github.com/tlk00/BitMagic/blob/master/src/bmfio.h) interfaces support buffered serialization and deserialization through seekable binary C++ streams. Direct serialization avoids constructing a complete temporary RAM BLOB; serializers retain working memory that can be reused across calls.

Bookmarks and optional persistent deserialization indexes help locate relevant encoded blocks. Memory-mapped BLOBs can serve as sources for repeated gathers, allowing applications to keep large datasets on SSD and retrieve subsets on demand. Access is block-oriented, and selected values still require decoding.

The provided stream adapter requires seeking, including header and bookmark patching on output. It is not a forward-only network writer. Custom integrations must honor the channel's positioning, buffering, completion, and error-handling contracts; these are documented alongside the implementation.

See [bit-vector file I/O and indexed gather](https://github.com/tlk00/BitMagic/tree/master/samples/bvsample27), [memory-mapped string retrieval](https://github.com/tlk00/BitMagic/tree/master/samples/strsvsample10), and [persistent float retrieval indexes](https://github.com/tlk00/BitMagic/tree/master/samples/svfsample05).

### Format compatibility

New readers retain support for older serialized data. Older readers are not guaranteed to understand newer format versions. BitMagic 9.3.1 adds explicit sparse-vector logical-size information, preserving trailing NULL positions, using schema versions 3 and 4 for the respective addressing modes.

Persisted deserialization indexes belong to a specific serialized representation and must be rebuilt when that representation changes.

## Getting started

BitMagic is a header-only C++17 library. Add the [`src`](https://github.com/tlk00/BitMagic/tree/master/src) directory to your compiler's include path and include the headers for the containers and algorithms you use.

| Header | Purpose |
|---|---|
| `bm.h` | Core bit-vector container |
| `bmaggregator.h` | Operations over groups of bit-vectors |
| `bmsparsevec.h` | Integer sparse vectors |
| `bmsparsevec_compr.h` | Rank-select compressed vectors |
| `bmstrsparsevec.h` | String vectors |
| `bmsparsevec_float.h` | Floating-point vectors |
| `bmsparsevec_algo.h` | Scanners and sparse-vector algorithms |
| `bmserial.h` | Bit-vector serialization and deserialization indexes |
| `bmsparsevec_serial.h` | Sparse-vector serialization |
| `bmsparsevec_float_serial.h` | Floating-point vector serialization |
| `bmfio.h` | C++ stream I/O |
| `bmintervals.h` | Interval operations and traversal |
| `bm3vl.h` | Three-valued logic |

To build an example using CMake:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target bvsample27
```

The [sample catalogue](https://github.com/tlk00/BitMagic/tree/master/samples) is the recommended starting point for learning usage patterns. Stress tests deliberately exercise unusual and inefficient cases and should not be treated as application templates.

## Platforms and configuration

BitMagic is a **header-only C++17 library**: add its headers to your project's include path and use them directly. No separate BitMagic library binary is required. It can be integrated into an existing application, IDE, or build environment without adopting the repository's build system.

Configure addressing, SIMD, and other supported options through preprocessor definitions in your build environment or in a common configuration header included before BitMagic headers. Keep these settings consistent across translation units.

See the [build and configuration guide](docs/build.md) for direct compiler commands, application CMake integration, repository Make/CMake builds, platform settings, and the configuration macro reference.

## Documentation, validation, and releases

The [BitMagic website](https://bitmagic.io) complements the source repository with design explanations, methodology, benchmark studies, and application notes. Choose a starting point according to what you want to understand:

| Topic | Reading |
|---|---|
| How the representations fit together | [Compression design overview](https://bitmagic.io/design) |
| How data can be searched without a separate secondary index | [Search with sparse vectors](https://bitmagic.io/sparse-vector-search) |
| Searchable string collections and scientific identifiers | [Dictionary compression and search](https://bitmagic.io/star-search) |
| Float representation, range search, and measured workloads | [Sparse vector floats on financial datasets](https://bitmagic.io/bm-svf) |
| Correlated columns and alignment compression | [XOR compression of bit-transposed vectors](https://bitmagic.io/bm-xor) |
| Application architectures and worked studies | [Use cases and design patterns](https://bitmagic.io/use-case) |
| Further algorithms and optimization material | [Technical articles](https://bitmagic.io/articles) |
| Boolean retrieval foundations | [Algebra of sets tutorial](https://bitmagic.io/set-algebra.html) |
| API details and runnable code | [API documentation](https://bitmagic.io/doxygen/html/modules.html) and [examples](https://github.com/tlk00/BitMagic/tree/master/samples) |
| Version-specific changes and distributions | [GitHub releases](https://github.com/tlk00/BitMagic/releases) |

Some technical articles document earlier releases. Use them for design rationale and workload-specific benchmark evidence; consult the current headers, examples, and release notes for current API behavior.

The project includes randomized stress tests and performance tests covering different representations, addressing modes, and build configurations. Use a tagged release for reproducible deployments, and evaluate memory consumption and performance with representative workloads.

## License and contributions

BitMagic is distributed under the [Apache License 2.0](https://github.com/tlk00/BitMagic/blob/master/LICENSE).

**Project attribution requirement:** Use of BitMagic requires explicit mention of the BitMagic project in derived materials, including the product or project documentation and published materials that incorporate or describe work based on the library. Include a reference and link to [BitMagic](https://bitmagic.io) on your product or project page.

Contributions to the library, examples, documentation, tests, and language bindings are welcome through the [GitHub project](https://github.com/tlk00/BitMagic).

Follow [BitMagic on Twitter/X (@bitmagicio)](https://twitter.com/bitmagicio) for project news and updates.
