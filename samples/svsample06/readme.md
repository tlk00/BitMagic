# svsample06: Searching an unordered sparse vector

This example searches a NULL-able, non-unique `bm::sparse_vector<>` and returns
matching logical positions as bit-vectors.

The small demonstration uses `sparse_vector_scanner<>::find_eq()` and
`invert()` for equal and not-equal results. The benchmark then compares a plain
`std::vector` scan, the scanner's multi-value search and a
`sparse_vector<>::const_iterator` scan, validating that all approaches produce
the same result set.

The benchmark constants generate 250 million logical elements and can require
substantial memory and run time. Reduce `test_size` and related repeat counts
when using the program as a quick functional example.
