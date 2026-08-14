# bvsample24: Splitting by rank

This example uses `bm::rank_range_split()` to partition a bit-vector into
closed index ranges containing approximately equal numbers of set bits.

The motivating case is parallel post-processing of a query result represented
as a bit-vector: if each result has similar cost, equal-population ranges make
more balanced jobs than equal-width index ranges.

After computing the target rank from the total population count, the sample
prints each returned range and uses `get_enumerator(range_start)` to traverse
only the set bits that belong to it.
