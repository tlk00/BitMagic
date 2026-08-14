# bvsample12: Setting and clearing bits

This example surveys `bm::bvector<>` mutation APIs and benchmarks several bulk
construction strategies.

The functional part covers initializer lists, `set()`, `set_bit()`,
`set_bit_conditional()`, `set_range()`, `clear_bit()`, `reset()`, `flip()`,
vector `swap()`, `extract_next()` and `bm::combine_or()`. It also demonstrates
`set_bit_no_check()`; this faster low-level method requires an explicit
`init()` before it is used.

The benchmark compares repeated `set_bit()`, `set_bit_no_check()`,
`bm::combine_or()` over STL ranges and the array form of `bvector<>::set()`.
The source enables `BM64ADDR`, so its index type follows BitMagic's 64-bit
address configuration.
