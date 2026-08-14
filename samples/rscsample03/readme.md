# rscsample03: RSC vector traversal

This example traverses a NULL-able `bm::rsc_sparse_vector<>` with its
`const_iterator`.

It demonstrates iteration from `begin()` to `end()`, NULL checks, constructing
an iterator at a specific logical position, `valid()`, `advance()` and
repositioning with `go_to()`. It also shows `bm::for_each_sparse()` as a
visitor-based full-scan interface that uses internal decode buffers.

The call to `sync()` after buffered construction is essential: it rebuilds the
rank-select index required by the iterator and bulk decode methods.
