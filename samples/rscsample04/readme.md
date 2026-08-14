# rscsample04: Predefined non-NULL positions

This example constructs `bm::rsc_sparse_vector<>` from a bit-vector that
defines every logical position allowed to contain a value.

After `sync()` builds the rank-select index, the sample uses `set()` and
`inc()` only at those predefined non-NULL positions. Because these writes
change values but not the NULL shape, the rank-select structure remains valid,
which is verified with `in_sync()`.

This pattern is useful when the sparse shape is known in advance and values
must be updated repeatedly without rebuilding rank-select metadata. Adding or
removing non-NULL positions would require synchronization again.
