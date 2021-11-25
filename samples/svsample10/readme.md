# Scanner searches 

BitMagic implements fast index-free searches on its succinct memory structures.

This example illustrates GT, GE, LT, LE and closed range searches using `bm::sparse_vector_scanner<>` - fast search solver which uses optimized logical operations between bit-slices to come up with a result bit-vector.

This example is not intended to benchmark index-free search, just show cases the usage basics.

