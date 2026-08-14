# bvsample01_64: 64-bit address mode

This example shows how to use `bm::bvector<>` with indexes beyond the default
32-bit address space. It enables the mode by including `bm64.h`, creates a
vector containing `bm::id_max - 1`, finds its active range and traverses it
with the 64-bit enumerator type.

It also demonstrates a full bit-vector, clearing its last addressable bit and
combining 64-bit vectors with logical AND.

Important constraints shown in the source:

- `BM64ADDR` or `bm64.h` must be selected consistently for the compilation
  unit;
- 32-bit and 64-bit BitMagic modes cannot be mixed in one compilation unit;
- the current implementation uses a 48-bit internal address space; and
- a 32-bit serialized vector can be read in 64-bit mode, but the reverse is
  not supported.
