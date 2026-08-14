# strsvsample01: String sparse-vector basics

This example introduces `bm::str_sparse_vector<>`, a bit-transposed succinct
container for zero-terminated strings.

It demonstrates setting, assigning, appending, inserting and erasing values;
automatic resizing; memory optimization; random access; and bulk traversal
with `const_iterator`. The configured string length is an initial capacity,
while storage is allocated sparsely.

A second vector enables `bm::use_null` and shows `try_get()`, explicit NULL
insertion, iterator positioning with `go_to()` and two ways to recognize NULL
values during traversal.
