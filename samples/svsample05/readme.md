# svsample05: Set-to-set remapping

This example uses a `bm::sparse_vector<>` as a translation table and
`bm::set2set_11_transform<>` to map every set bit in one `bm::bvector<>` to an
index in another.

With `bm::use_null`, unassigned table positions have no mapping and their input
bits are omitted from the result. Without NULL support, unassigned sparse
positions read as zero, so the same input is considered mapped to index 0.

The sample demonstrates both complete set transformation with `run()` and a
single-value lookup with `remap()`, making the effect of NULL semantics on the
image of a set explicit.
