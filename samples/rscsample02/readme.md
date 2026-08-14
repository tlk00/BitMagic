# rscsample02: Selective RSC deserialization

This example builds a NULL-able `bm::rsc_sparse_vector<>` with its buffered
back insert iterator, calls `sync()` to prepare rank-select access and
serializes it with block bookmarks.

It then restores only selected data from the compressed BLOB in two ways:

- a bit-vector mask passed to `deserialize()`; and
- a closed interval passed to `deserialize_range()`.

The destination retains the original logical size and NULL shape, while
unselected assigned positions contain zero. Bookmarks trade a small increase
in serialized size for faster skipping during repeated range reads.
