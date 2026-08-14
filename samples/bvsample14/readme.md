# bvsample14: Operations on serialized bit-vectors

This example performs set algebra directly between a live `bm::bvector<>` and
a serialized, compressed bit-vector BLOB.

It first optimizes and serializes two vectors, then demonstrates ordinary
deserialization followed by `bm::operation_deserializer<>` operations. The
examples include AND and SUB, plus `set_COUNT_SUB_AB` to compute the result
population count without materializing the result vector.

The program checks that BLOB-based AND reproduces the expected source vector
and that count-only subtraction matches the count of a fully materialized
subtraction. This approach is useful when compressed vectors are stored in a
file, database or network buffer and do not need to be expanded first.
