# bvsample16: Group operations with `bm::aggregator<>`

This example uses `bm::aggregator<>` to execute logical operations over groups
of bit-vectors. The aggregator applies cache blocking and memory-bandwidth
optimizations that are useful for multi-vector operations.

It demonstrates:

- adding non-owned vector pointers to an argument group;
- `combine_or()` and `combine_and()` over the same group;
- `reset()` before configuring a new operation; and
- fused `combine_and_sub()`, which computes the intersection of group 0 and
  subtracts the union of group 1.

The aggregator does not take ownership of vectors passed to `add()`, so every
source vector must remain alive until the operation or reset completes.
