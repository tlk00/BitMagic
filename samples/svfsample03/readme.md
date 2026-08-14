# svfsample03: Float-vector comparison and composition

This example demonstrates operations that clear, compare and combine
`bm::sparse_vector_float<>` containers.

It covers full `clear()`, zeroing a closed interval with `clear_range()`,
vector equality and element-to-float `compare()` with a tolerance. It then
shows `join()`, destructive `merge()`, `copy_range()` and conversion to
immutable storage with `freeze()` and `is_ro()`.

`join()` combines the underlying bit-planes rather than performing numerical
addition; overlapping nonzero IEEE-754 components can therefore produce an
invalid or unintended float. `merge()` may also modify the source because it
borrows its storage.
