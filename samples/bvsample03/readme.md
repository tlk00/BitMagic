# bvsample03: Compression strategies and memory statistics

This example fills two bit-vectors with the same dense distribution while
using different allocation strategies. One vector uses the default strategy;
the other requests GAP/RLE blocks with `set_new_blocks_strat(bm::BM_GAP)`.

`calc_stat()` reports the bit count, BIT blocks, GAP blocks, memory use and
maximum serialization buffer size for both representations. The default
vector is then passed through `optimize()` so its post-compression footprint
can be compared with the vector that used GAP blocks from the start.

The sample is intended to illustrate the memory trade-off of block strategy;
the exact numbers depend on the generated bit distribution and build.
