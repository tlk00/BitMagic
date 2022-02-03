# sample26

This example shows how to turn `bm::bvector<>` into a `READ-ONLY` immutable vector.


## Why?

There are a few reasons to use immutable structures.

Default, mutable `bm::bvector<>` is always making some internal memory space allocation reservations to be able to chnage without calling allocations too often (can be a performance killer). 

Many systems and applications (especially in the area of information retrieval or large data analysis) have a point when data structures are complete and can be made read-only (if it has any intrinsic benefits). 

Succinct data structures for example can be farther optimized to improve memory performance once it is a read-only. 

Another advantage of READ-ONLY mode is heap defragmentation. BitMagic vectors on construction allocate in fixed block sizes, this strategy minimizes heap fragmentation but only to an extent. Turning read-only we can not just avoid edit space reservations, but we can also place blocks of memory togther in some arena space, so that logical operations on bit-vectors go faster, because adjacent memory would be better pre-fetched and preserve CPU cache coherency.

### Recap

Immutable bit-vectors:

- reduce memory consumption. 
This example shows how to measure it using synthetic yet illustrative bit distributions
- improve performance of your application by contiguous memory allocation versus heap fragmentation (this example does not benchmark it).

### What if I need bit-vector to be writable again, can it do it?

Yes, bit-vector can be turned mutable or you can construct a mutable copy, this example shows how to do it.


## Summary 

Exmaple shows ways how to memory compress bit-vectors, measure its memory allocation statistics, turn bit-vector to read-only state and convert it back to writable if necessary.