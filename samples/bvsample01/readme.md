# bvsample01: Basic bit-vector operations

This example introduces `bm::bvector<>` by building a small set of integer
indexes and reading the set bits back.

It demonstrates:

- initializer-list construction and `set()`;
- population counting with `count()`;
- traversal with `get_first()` and `get_next()`;
- clearing a vector and using the `operator[]` bit reference;
- checking bits with `test()`; and
- exchanging two bit positions with `swap(index1, index2)`.

The assertions at the end verify that swapping positions changes the bit
values without changing the vector's population count.
