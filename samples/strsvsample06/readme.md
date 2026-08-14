# strsvsample06: Substrings, mismatch and search

This example works with a NULL-able string sparse vector containing formatted
telephone-like values.

It demonstrates full and substring-mode iterators, direct access to the
NOT-NULL bit-vector, and `sparse_vector_find_first_mismatch()` for locating the
first difference between two vectors. Substring iterators transpose only the
requested character range and can avoid extracting complete strings.

The search section reuses `bm::sparse_vector_scanner<>` for exact string
matching and prefix matching. Matches are returned as a bit-vector of logical
positions and traversed with a set-bit enumerator.
