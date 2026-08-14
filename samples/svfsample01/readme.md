# svfsample01: Floating-point sparse-vector basics

This example introduces `bm::sparse_vector_float<>`, which represents
floating-point values as separate succinct sign, exponent and mantissa
components.

It demonstrates `empty()`, `size()`, `push_back()`, random `set()` and `get()`,
bulk `import()`, `optimize()` and `swap()`. Assigning beyond the current end
automatically grows the vector, with intervening positions reading as zero.

The two demo functions cover both individual edits and array-oriented loading,
then print the values before and after exchanging two float vectors.
