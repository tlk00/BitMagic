# svsample07: Lower-bound insertion sort

This example builds a sorted `bm::sparse_vector<>` from shuffled unsigned
integers.

For each value, a reusable `bm::sparse_vector_scanner<>` performs binary
lower-bound search with `bfind()` and returns the insertion position.
`sparse_vector<>::insert()` then keeps the succinct vector ordered. After all
edits, `optimize()` recompresses its bit-planes.

The finished vector is compared element by element with `std::sort()` output.
The example illustrates ordered search and insertion in compressed memory;
bulk loading pre-sorted data remains preferable when that option is available.
