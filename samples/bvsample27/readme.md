# bvsample27: file I/O and gather deserialization

This example fills a `bm::bvector<>`, uses `bm::random_subset` to select one
quarter of its set bits into a second vector, and writes the vectors to two
separate files through `bm::streams_encoder` from `bmfio.h`.
It then restores each vector through `bm::streams_deserializer`, and gathers
selected positions from both its file and an equivalent RAM BLOB.

**One serializer is constructed and reused for both files.** Each file gets its
own stream and encoder. The first channel is finished and its file closed before
the second file is opened. The serializer retains its working buffer between
these calls, avoiding a new serializer and repeated scratch-buffer allocation.

## Interface demonstrated

```cpp
bm::random_subset<bm::bvector<>> sampler;
sampler.sample(second, first, first.count() / 4);

bm::serializer<bm::bvector<>> serializer;
serializer.set_bookmarks(true, 16);
// Do not pin a compression level: retain the library default.

const bm::bvector<>* vectors[] = {&first, &second};
for (unsigned i = 0; i < 2; ++i)
{
    std::ofstream file(paths[i], std::ios::binary | std::ios::trunc);
    bm::streams_encoder output(file);
    if (!serializer.serialize(*vectors[i], output) || !output.finish())
        throw std::runtime_error("File serialization failed");
}
```

The complete source creates the paths in a private temporary directory, checks
file opening and closing, verifies that the second vector is a subset with the
requested cardinality, and cleans up the files.

- The serializer owns a reusable BM allocator-typed working buffer. The channel's
  RAM encoder performs existing block encoding, and the serializer flushes only
  at safe boundaries. The write path allocates no complete serialized RAM BLOB.
- Successful `serialize()` means the object and its patches were accepted. It
  does not call `finish()`.
- Successful `flush()` makes working memory immediately reusable. A future
  asynchronous adapter must own queued bytes before returning.
- The caller invokes `finish()` for each channel before closing that channel's
  file. It completes submitted work but does not close the file or guarantee
  disk durability.
- Header and bookmark patches require binary, seekable output without
  `std::ios::app`. Forward-only sockets and compressed streams are not supported
  by this file adapter.
- Errors return `false` and leave the channel's `is_good()` false. Stream
  exceptions propagate if enabled. The sample reports failures, does not retry
  failed I/O, and attempts to remove partial files during exception unwinding.

Each file contains one ordinary BitMagic BLOB, identical to RAM serialization
of the same vector with the same settings. There is no additional framing or
collection header. Bookmarks use interval 16; compression stays at the default
for forward compatibility.

## Stream deserialization

`deserialize_file()` shows a complete restore using a seekable binary input:

```cpp
std::ifstream file(path, std::ios::binary);
bm::streams_decoder input(file);
bm::streams_deserializer<bm::bvector<>> reader;
bm::bvector<> restored;
if (!reader.deserialize(restored, input))
    throw std::runtime_error("File deserialization failed");
```

The reader owns reusable temporary memory instead of loading the complete BLOB
into RAM. The sample reuses one reader across both files. A successful call
leaves the underlying stream immediately after the BLOB, including when input
buffering has read ahead. Deserialization into a nonempty target combines bits
with OR; these helper functions take a mutable output parameter and clear it
before deserialization. Gather output must be distinct from the selection mask.

## Gather deserialization: restore only requested positions

A selection mask describes logical bit positions wanted by the application.
Gather deserialization restores the blocks containing those positions, then
an AND with the mask removes unwanted bits within those blocks. Its exact result
is therefore **source AND mask**. The mask need not be a subset of the source:
requested positions where the source has a zero remain zero.

The sample creates disjoint intervals over nearby and distant blocks.
`build_file_index()` scans each file once to construct a
`bm::deserialization_index`, recording marker/bookmark offsets without restoring
the complete vector. `save_dindex()` serializes the small index into the
library's RAM buffer and writes its standalone record to a `.didx` file.
`load_dindex()` reads that file into RAM and restores it with
`bm::deserialization_index_deserializer`. The sample checks the restored index
with `equal()` and deliberately uses it for both gather paths.

`gather_file()` uses the restored index and a block digest of the mask:

```cpp
bm::bvector<> digest, gathered;
mask.build_block_digest(digest);
reader.set_deserialization_index_use(&index);
reader.set_block_digest_vector_use(&digest);
const bool ok = reader.deserialize(gathered, input);
reader.unset_block_digest_vector();
reader.unset_deserialization_index();
if (!ok)
    throw std::runtime_error("File gather deserialization failed");
gathered &= mask;
```

The digest selects blocks of 65,536 logical bits; it does not encode the exact
within-block selection. The final AND is essential. The index describes positions
in one serialized representation, not values or an independently searchable
database. Reuse it for repeated queries against the same BLOB; changing the
serialized bytes requires rebuilding it.

Index persistence currently uses RAM serialization rather than the streaming
adapter. This is a reasonable tradeoff because the index is much smaller than
the bvector BLOB it describes. The standalone index format has its own header;
the load helper checks that deserialization consumes the complete file.

`gather_ram()` demonstrates the same index/digest setup with the existing
`bm::deserializer<bvector_type, bm::decoder>` and
`reader.deserialize(gathered, blob.data())`. The sample serializes the source
into a RAM buffer with identical settings, so its bytes and index offsets match
the file. For both vectors it explicitly checks:

```cpp
bvector_type expected(source, bm::finalization::READWRITE);
expected &= mask;
const bool file_equal = expected.equal(from_file);
const bool ram_equal = expected.equal(from_ram);
```

Both comparisons must print `equal()=true`; a mismatch fails the sample.

### Why use gather for files?

- **I/O:** index/bookmark offsets allow seeking over unrequested serialized
  records instead of loading the entire BLOB into a RAM buffer.
- **CPU:** skipped blocks avoid reconstruction and subsequent filtering of a
  complete source vector.
- **Memory:** only selected blocks are materialized, alongside the mask, digest,
  index and reusable input buffer. The source vector and full restore kept here
  for demonstration and comparison are not prerequisites for file gathering.

These benefits depend on the selection and serialization layout. Building the
index costs an initial scan, best amortized across multiple queries. Buffer
read-ahead can fetch unrequested bytes, and the current reader scans the remaining
tail to locate the exact BLOB end after the last selected region. Consequently,
gather does not promise to read only requested bytes or to outperform a full
restore for every query. RAM gather can save decoding and result allocations,
but its complete serialized BLOB is already resident; it saves no file-loading
I/O by itself.

## Build and run

From this directory:

```sh
make
./sample27 --verify
```

Or compile directly:

```sh
c++ -std=c++17 -O2 -I../../src sample27.cpp -o sample27
./sample27 --verify
```

From the repository root:

```sh
cmake -S . -B build
cmake --build build --target bvsample27
./build/bin/bvsample27 --verify
```

The executable accepts `[--verify] [directory]`. The optional directory is an
existing parent directory; it defaults to the working directory. The program
creates a unique `sample27-*` subdirectory containing `first.bv` and `subset.bv`,
plus `first.didx` and `subset.didx`. It prints their sizes, then **deletes all
four files and the private directory**.
Existing files in the parent are not overwritten or removed. Cleanup is also
attempted on errors; deletion failures are reported.

Use `--help` for usage. Success returns zero; argument, I/O, or verification
errors return a nonzero exit status.

## Optional compatibility verification

`--verify` reads each completed file into RAM, compares its length and bytes
with the existing RAM serializer, and restores its vector through the RAM
deserializer. Both files are verified before deletion. It prints:

```text
Verification OK: exact RAM bytes and both vectors restored
Deleted all demo files and their temporary directory
```

The RAM gather demonstration always creates a complete RAM BLOB. `--verify`
additionally loads each file for byte-for-byte compatibility checks. The file
serialization and gather functions themselves do not need either full-file
buffer. The sampled subset can vary with the random generator; verification
compares against the actual generated vectors rather than fixed expected bytes.

## 64-bit addresses

```sh
c++ -std=c++17 -O2 -DBM64ADDR -I../../src sample27.cpp -o sample27_64
./sample27_64 --verify
```

The source pattern starts at `2^40` in this build, and the random subset retains
those high bit positions. Verification uses the matching 64-bit RAM interfaces.

## Related serialization examples

- [bvsample04](../bvsample04/readme.md): RAM serialization, buffer ownership, and ordinary deserialization.
- [bvsample14](../bvsample14/readme.md): set algebra and count operations on serialized RAM BLOBs.
- [bvsample22](../bvsample22/readme.md): bookmarks and selective range deserialization.
- [bvsample01_64](../bvsample01_64/readme.md): large-address mode basics; build this sample with `-DBM64ADDR` to serialize patterns above `2^32`.

The RAM, operation, and range deserializers consume an in-memory BLOB.
This example preserves that format and adds buffered file reading through
`bm::streams_decoder` and `bm::streams_deserializer`.
