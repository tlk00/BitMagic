# strsvsample04: NULL-able string vectors

This example constructs `bm::str_sparse_vector<>` with `bm::use_null` and
demonstrates unassigned string values.

It adds NULLs through a back insert iterator, a null pointer and `set_null()`;
checks values through random access and `const_iterator`; and uses
`clear_range(..., true)` to turn a closed range into NULLs. It also applies
`set_null()` with a bit-vector of logical indexes for efficient bulk updates.

The final traversal uses `get_string_view()`: a view whose data pointer is null
represents an unassigned element. `optimize()` recompresses the changed
bit-planes and can release blocks made empty by bulk clearing.
