## sample25

This example illustrates various bit-vector traversal techniques when we need 
to get access to each and every set bit (set index) in bit-vector.

Testing each bit is obviously very naive and inefficient method, 
so this example does not even try it.

The default and often preferred method is to use `bm::bvector<>::enumerator`
which is an iterator class to travel decode bit-vector.
While convenient to use if cannot be very fast, because it is a complex object 
has to safe/restore state on each step so all of this produces a significant 
challange to even modern compilers. 

This example shows how to use C-style callback and C++ visitors to traverse bit-vectors.
Algorithms `bm::visit_each_bit` and `bm::for_each_bit` make C++ compiler brew different 
code more based on multiple nested loops, rather than iterator state machine, whih can be
quite more efficient. 

`bm::for_each_bit` exposes some internal details (not too comlicated) and should 
be considered as the fastest method for bit-vector traversal.

`bm::visit_each_bit` is an adaptation for C-style code for cases where compiled code size
is important.

### Ranges traversals

`bm::bvector<>::enumerator` can be positioned for range traversal using
`bm::bvector<>::get_enumerator` or with `bm::bvector<>::enumerator::go_to(idx)`.

BitMagic offers efficicnt algorithms for [from..to] closed range traversals.
See `bm::for_each_bit_range` and `bm::visit_each_bit_range` usage in this example.

