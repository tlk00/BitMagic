# svsample08: Selective sparse-vector deserialization

This example serializes a NULL-able `bm::sparse_vector<>` and restores selected
logical positions from the compressed BLOB.

A buffered back inserter creates the sample vector, including explicit NULLs.
The serializer enables block bookmarks, which help the deserializer skip
unneeded regions. The program then compares mask-based `deserialize()` with
`deserialize_range()` for one closed interval.

The destination retains the source's logical size and NULL positions, while
unselected assigned values read as zero. Range deserialization is the clearer
and faster choice when the selection is contiguous.

The sample also shows index-assisted gather deserialization. A
`bm::sparse_vector_deserializer<>::deserialization_index_type` object is built
once from the serialized BLOB, optimized, and attached to a deserializer before
calling mask-based `deserialize()`.

The deserialization index is a compact navigation map for the serialized BLOB.
It records stream marker positions and bookmark landing points so repeated
selective deserialization calls can jump closer to requested blocks and avoid
decoding unrelated payloads.

This example is intentionally small and is not a benchmark. Its purpose is to
show the API sequence and object lifetime:

1. serialize with bookmarks enabled;
2. build and optimize the deserialization index for the BLOB;
3. keep the index alive while indexed deserializers use it;
4. pass a gather mask to restore only selected logical positions.

Larger sparse vectors with repeated sparse gather requests are the better use
case for this technique.

See `strsvsample10` for a more advanced example of the same technique using a
large `bm::str_sparse_vector<>`, a serialized file on disk, `mmap()` access and
index-assisted gather deserialization from the mapped BLOB.

The tradeoff is memory versus CPU. Keeping data in serialized compressed form
can substantially reduce resident memory, especially when only a small subset
is needed at a time. The cost is that requested values must be decoded from the
compressed BLOB on demand, so the application spends more CPU during retrieval.
Bookmarks and the deserialization index reduce unnecessary scanning and
decompression, but they do not make selective deserialization free.
