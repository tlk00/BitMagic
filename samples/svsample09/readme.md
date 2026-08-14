# svsample09: Finding the first mismatch

This example uses `bm::sparse_vector_find_first_mismatch()` to locate the first
logical difference between two sparse vectors.

The first case compares ordinary vectors, first while equal and then after one
value is changed. The second case uses `bm::use_null` to show that assignment
state is part of equality: an explicit numeric zero and a NULL element are a
mismatch even though random value access may produce zero for both.

The algorithm returns a Boolean indicating whether a difference exists and
writes its logical position to the supplied index variable.
