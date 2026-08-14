# rscsample01: Rank-select compressed sparse vectors

This example introduces `bm::rsc_sparse_vector<>`, a compact mostly read-only
representation for NULL-able scalar vectors. It removes unassigned values from
the bit-transposed value planes and maps logical positions through a
rank-select index.

The program builds a NULL-able `bm::sparse_vector<>`, loads it into an RSC
vector with `load_from()`, reads values with `get()` and `try_get()`, optimizes
and serializes the result, and validates deserialization.

Finally, `load_to()` expands the RSC vector back into an editable sparse
vector, demonstrating a complete mutable-to-succinct-to-mutable round trip.
