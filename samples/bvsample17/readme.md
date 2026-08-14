# bvsample17: Rank and select

This example builds a `bm::bvector<>::rs_index_type` with `build_rs_index()`
and uses it to answer rank and select queries.

`rank(position, index)` returns the number of set bits up to a position, while
`rank_corrected()` adjusts the result when the queried position itself is set.
`select(rank, position, index)` performs the inverse lookup and finds the bit
position associated with a one-based rank.

The rank-select index accelerates repeated queries on a stable vector. If the
vector is modified, rebuild the index before relying on it again.
