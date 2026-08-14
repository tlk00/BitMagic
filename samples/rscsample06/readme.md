# rscsample06: Gathering RSC values

This example uses `rsc_sparse_vector<>::gather()` to extract several logical
positions in one call, including positions requested in random order.

The vector is constructed from a predefined non-NULL mask, synchronized and
populated without changing that mask. The gather call receives an index array,
an equally sized scratch array for rank-select coordinate translation, an
output values array and the input sort-order hint.

For a requested NULL position, the corresponding scratch coordinate is set to
`bm::id_max`. If indexes are already ascending, passing `bm::BM_SORTED` instead
of `bm::BM_UNKNOWN` can select a faster path.
