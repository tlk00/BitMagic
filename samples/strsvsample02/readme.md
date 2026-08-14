# strsvsample02: Compressed insertion sort

This example maintains a sorted `bm::str_sparse_vector<>` while loading an
unsorted collection of numeric strings.

For each input string, a reusable `bm::sparse_vector_scanner<>` calls
`lower_bound_str()` to find the insertion position, and `insert()` places the
value there. The finished succinct vector is optimized and compared element by
element with the result of `std::sort()`.

The sample demonstrates how ordered string data can be built and searched in
compressed memory, although repeated middle insertion is inherently more
expensive than bulk loading already sorted input.
