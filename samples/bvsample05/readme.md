# bvsample05: Enumerating set bits

This example uses `bm::bvector<>::enumerator` to decode the indexes of set
bits efficiently.

It covers traversal from `first()` to `end()`, use with `std::for_each()` and
random positioning with `go_to()`. When positioned on an arbitrary index, the
enumerator advances to that bit if it is set or to the next available set bit.

The direct prefix-increment loop is the preferred traversal form in the
sample. Generic STL algorithms can perform additional iterator-distance work,
which may be more expensive for a bit-vector enumerator.
