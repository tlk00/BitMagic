# strsvsample03: Bulk loading, remapping and serialization

This example bulk-loads a large sorted string collection into
`bm::str_sparse_vector<>` with its buffered back insert iterator.

It measures memory use, remaps character codes to improve compression,
optimizes the bit-planes and demonstrates an alternative inserter that builds
the remap tables while loading. It also freezes a vector into immutable mode
and creates another vector with copied remap tables so related projections use
the same encoding.

The final section serializes the remapped vector, writes the BLOB to `test.sv`
in the current directory, deserializes it and verifies equality. Explicit
`flush()` calls make each buffered load's commit point visible.
