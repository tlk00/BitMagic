# bvsample07: Operations with arrays of indexes

This example applies set operations directly between a `bm::bvector<>` and a
plain integer array, including an unsorted array with duplicate values.

It demonstrates:

- `bm::combine_and()` for an unsorted input range;
- `bm::combine_or()`, which naturally produces a sorted unique set;
- `bm::combine_and_sorted()` after sorting the input for a faster
  intersection; and
- `bm::combine_sub()` for set difference.

These algorithms avoid first constructing a temporary bit-vector for the
array. Sorted-input variants should be preferred when the caller can provide
ordered indexes.
