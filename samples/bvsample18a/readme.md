# bvsample18a: Importing an external bit stream

This example imports a packed array of 32-bit words into `bm::bvector<>` with
`bm::bit_import_u32()` from `bmbvimport.h`.

Three bits are placed on both sides of a 64K boundary to show how word and
block offsets are translated into bit-vector indexes. The import call enables
on-the-fly optimization, allowing each completed block to be compressed while
its data is still cache-resident.

Assertions verify the imported population count and the exact destination
indexes.
