# bvsample08: STL interoperability

This example shows how `bm::bvector<>` enumerators and inserters participate
in standard C++ iterator algorithms.

The sample copies set-bit indexes into `std::vector` and `std::list`, combines
a bit-vector with a sorted STL vector through `std::set_union()`,
`std::set_intersection()` and `std::set_difference()`, and copies an unsorted
sequence into `bvector<>::inserter()`.

Because a bit-vector represents a set, inserting the source sequence also
orders the indexes and removes duplicates in the resulting vector.
