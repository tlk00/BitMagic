# svsample01: Integer sparse-vector basics

This example introduces `bm::sparse_vector<>`, a bit-transposed succinct
container for scalar integer values.

It bulk-imports C arrays, appends another array at an offset, optimizes the
bit-planes and reads values with `at()` and `const_iterator`. Separate demos
show unsigned and signed integer vectors; signed values are transformed from
their normal complementary representation into a layout better suited to
bit-plane compression.

For large batches, `import()` avoids the overhead of repeated random element
updates and is the preferred construction method shown by the sample.
