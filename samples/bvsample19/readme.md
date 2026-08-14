# bvsample19: Destructive merge

This example demonstrates `bvector<>::merge()`. Its set result is equivalent
to logical OR, but the operation may transfer memory blocks from the source
instead of allocating and copying them.

That optimization makes `merge()` destructive: after the call, the source
vector's content is undefined and it should be cleared, destroyed or rebuilt
before reuse. The destination contains the union of both original sets.

A typical use case is combining partial bit-vectors produced by independent
or partitioned workers when those temporary vectors are no longer needed.
