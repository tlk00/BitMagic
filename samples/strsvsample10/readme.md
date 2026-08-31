# strsvsample10: mmap-backed gather deserialization

This example demonstrates how to keep a large `bm::str_sparse_vector<>` as a
serialized file and later retrieve a small fraction of its elements with gather
deserialization.

The intended use case is a large string vector that is queried only
occasionally. Instead of keeping the succinct vector resident in RAM, an
application can prepare and serialize it once, evict the live vector, and later
use the serialized BLOB on disk as the retrieval source.

## Workflow

The sample is written as two stages.

Stage 1 is preparation. The program constructs a large string sparse vector,
remaps and optimizes it, serializes it with bookmarks enabled, and saves the
serialized BLOB to `strsvsample10_data.bm`.

Stage 2 is retrieval. The program opens the saved file read-only, maps it with
`mmap()`, constructs a deserialization index from the mapped BLOB, and performs
several gather deserializations. The original source vector is not used during
this stage.

In a real application these stages can be separated. Stage 1 can run offline on
a machine with enough memory to build the vector. Stage 2 can run later in a
lower-RAM process that only maps the serialized file and retrieves selected
elements.

## Platform requirements

The mmap part of the sample uses POSIX APIs: `open()`, `fstat()`, `mmap()`,
`munmap()`, and `close()`. It is compiled on Unix-like platforms, including
Linux and macOS.

On platforms where these APIs are not available, the program still builds. It
prints a short message and exits successfully without trying to serialize or map
the file.

## Serialization

The source vector is filled through `str_sparse_vector<>::back_insert_iterator`,
then `remap()` and `optimize()` are called before serialization. Remapping
reduces the character alphabet in each string plane, and optimization compresses
the bit-vector planes.

The sample uses `bm::sparse_vector_serializer<>` with bookmarks enabled:

```cpp
serializer.set_bookmarks(true, 16);
serializer.serialize(str_sv, layout);
```

Bookmarks add skip points to the serialized BLOB. They increase serialized size
slightly, but make later range and gather deserialization more efficient.

## mmap access

After serialization, the sample writes the BLOB to disk and releases the source
vector and temporary serialization layout. The file is then opened and mapped
read-only:

```cpp
mapped = mmap(0, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
```

The mapped address is passed directly to BitMagic deserialization APIs as a
serialized buffer pointer.

## Gather ranges

Gather deserialization is controlled by a `bm::bvector<>` mask. Every set bit in
the mask identifies one logical string-vector element to retrieve.

Each request in the sample creates three small islands. Every island contains
roughly 128-256 adjacent elements. This models sparse access to a few localized
areas of a much larger vector.

## Deserialization index

The deserialization index is built from the serialized BLOB:

```cpp
deserializer.construct_deserialization_index(dindex, mapped.data());
dindex.optimize();
```

The index records serialized marker offsets and bookmark landing points. It acts
as a map of the serialized stream, allowing index-assisted gather
deserialization to jump near requested blocks instead of scanning and decoding
large unrelated regions.

This reduces two costs:

- disk paging, because fewer cold serialized pages need to be touched;
- CPU work, because unrelated compressed blocks are skipped instead of decoded.

The index itself uses RAM, so this technique trades a small auxiliary structure
for lower memory residency of the large vector and more selective access to the
serialized representation.

## Gather deserialization

For every request, the sample creates a fresh output vector, attaches the
deserialization index, enables index use, and deserializes only the masked
positions from the mapped BLOB:

```cpp
deserializer.set_deserialization_index(&dindex);
deserializer.set_deserialization_index_use(true);
deserializer.deserialize(str_sv_out, mapped.data(), mask_bv);
```

The output vector contains the selected logical positions. The sample prints the
requested island ranges but intentionally avoids correctness checks or timing
comparisons so the API flow remains easy to follow.

## Understanding the output

A typical run prints:

```text
Stage 1: constructing source string sparse vector
  source vector size = 4194304 elements
  plain std::vector<std::string> estimate = 160.70 MB
  optimized str_sparse_vector<> memory = 30.05 MB
Stage 1: serializing vector with bookmarks
  serialized file = strsvsample10_data.bm
  serialized BLOB size = 17.88 MB
  serialized BLOB is 59.5% of optimized str_sparse_vector<> memory
Stage 2: opening serialized BLOB with mmap
  mapped BLOB size = 17.88 MB
Stage 2: building deserialization index
  index resident memory = 0.08 MB (0.444% of mapped BLOB)
Stage 2: gather deserialization from mapped BLOB
  gather 1: requested 576 elements in [23757..23948] [1080252..1080443] [2136747..2136938]
  ...
  gather 10: requested 576 elements in [237570..237761] [1365336..1365527] [2493102..2493293]
Stage 2: resident memory used by selected results
  average gathered str_sparse_vector<> memory = 0.64 MB
Summary: storage and resident-memory progression
  plain std::vector<std::string> estimate = 160.70 MB
  optimized str_sparse_vector<> memory = 30.05 MB
  serialized BLOB on disk = 17.88 MB
  resident deserialization index = 0.08 MB
  average gathered output vector = 0.64 MB
Removed generated file: strsvsample10_data.bm
```

The first stage confirms that the sample built a logical string vector with
4,194,304 elements. This is the high-RAM preparation step: the vector exists in
memory while it is being populated, remapped, optimized and serialized.

The plain `std::vector<std::string>` number is an estimate, not an allocator
measurement. It adds the vector object, one `std::string` object per element,
and the generated character data. It provides a simple baseline for how much
RAM a direct STL representation would need.

The optimized `str_sparse_vector<>` memory is measured with BitMagic
`calc_stat()` after remapping and optimization. It shows the Stage 1 in-memory
succinct representation: still randomly accessible, but much smaller than the
plain string-vector model.

The serialized BLOB size is the compressed representation written to disk. It
is usually smaller than the optimized succinct vector, but it is no longer the
normal mutable in-memory container. The BLOB is intended for persistence,
transfer, or selective deserialization. The percentage line compares the BLOB
size with the optimized source vector memory to show the additional compression
achieved by serialization.

The second stage opens the same file with `mmap()`. The mapped size should match
the serialized size because the mapped region covers the whole BLOB. At this
point the sample is no longer using the original `str_sparse_vector<>`; it is
using the serialized file as the data source.

The index resident memory line describes the small in-memory navigation
structure built over the mapped BLOB. The sample prints this as megabytes and
as a percentage of the mapped BLOB. This is the main tradeoff demonstrated by
the sample: keep a small map of the serialized stream in RAM instead of keeping
the full source vector resident.

Each gather line shows one retrieval request. The sample asks for three islands
per request, and each island contains 192 adjacent logical elements. That is why
every request reports 576 elements total. The bracketed ranges are positions in
the original logical vector, not byte offsets in the serialized file.

For example:

```text
gather 1: requested 576 elements in [23757..23948] [1080252..1080443] [2136747..2136938]
```

This means the mask requested three small ranges from widely separated parts of
the vector. The deserializer uses the mask, the deserialization index and the
serialized bookmarks to reconstruct only those selected positions into the
output `str_sparse_vector<>`.

The average gathered output memory reports the measured memory of the small
`str_sparse_vector<>` objects produced by the gather requests. In a service
that retrieves only a small fraction of a large vector, this output vector plus
the deserialization index can be the main resident BitMagic memory during
retrieval.

The summary repeats the high-level progression: plain strings, optimized
succinct vector, serialized BLOB on disk, the resident deserialization index,
and the average gathered output vector. This is the intended way to read the
sample output. It is not a performance benchmark; it is a memory-residency and
API workflow demonstration.

The final line reports cleanup of the generated BLOB. The sample removes the
file to avoid leaving large temporary data behind. If you want to inspect the
file manually, comment out the `std::remove()` call near the end of the sample.

## Design tradeoffs

This technique trades CPU and storage access for lower memory residency. A full
`str_sparse_vector<>` keeps the data ready for random access in RAM. The mmap
gather approach keeps only the mapped compressed BLOB, a small deserialization
index and the current gathered result vector resident in BitMagic structures.

That tradeoff is attractive when the data set is large and requests touch only
a small fraction of elements. It is less attractive when the application needs
frequent dense access, because repeated decompression and disk paging can cost
more than keeping the optimized succinct vector in memory.

The deserialization index and bookmarks reduce this penalty but do not remove
it completely. The index gives the deserializer a map of useful stream
positions, and bookmarks provide serialized skip points, so the code can avoid
walking most unrelated blocks. Some CPU work is still needed to decode selected
regions into an output vector.

Applications with predictable request patterns can add their own policy around
this sample. For example, a service can keep hot vectors resident, use the
mmap-backed path for cold vectors, or use worker threads to prefetch mapped
pages or prepare deserialization indexes before user requests arrive. Those
policies are intentionally outside the sample so the core API workflow remains
clear.

## Limitations

This is an API clarity sample, not a benchmark. It does not compare mmap access
with in-memory BLOB access, does not validate every restored string against the
source vector, and does not tune bookmark intervals for a particular storage
device.

For production use, choose the vector size, bookmark interval, and index
lifetime based on request selectivity, available RAM, storage latency, and the
expected reuse of the mapped BLOB.
