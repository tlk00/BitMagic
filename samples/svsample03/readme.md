# svsample03: Joining and decoding sparse vectors

This example combines two `bm::sparse_vector<>` containers and extracts a
contiguous block of scalar values.

It demonstrates array import, memory optimization, partial value access with
`get_unsigned_bits()`, dynamic growth through `set()` and bit-plane union with
`join()`. The two source vectors occupy separated logical regions, so joining
them produces one vector containing both regions.

Finally, `decode()` transposes 16 consecutive elements into an ordinary
`std::vector<unsigned>`, which is more efficient than repeated random reads for
a range.
