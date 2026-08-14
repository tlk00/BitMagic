# bvsample11: Population-count techniques

This example benchmarks several ways to count set bits in a large mixed
BIT/GAP bit-vector.

It compares:

- full-vector `count()`;
- `count_range()` with and without a prebuilt rank-select index;
- accelerated `count_to()` and a range count formed from two prefix counts;
- `bm::count_and()` with a range mask; and
- `counted_enumerator`, which tracks a running count while traversing bits.

The program builds the acceleration structure with `build_rs_index()` and
prints all timings through `bm::chrono_taker`. A rank-select index is most
useful when a vector changes infrequently and serves many count queries; it
must be rebuilt after relevant mutations.
