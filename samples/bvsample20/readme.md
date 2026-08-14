# bvsample20: Right shifts and insertion

This example shows sequence-style insertion operations on a bit-vector.

`shift_right()` moves every set bit one position toward a larger index.
`insert(position, value)` creates a new bit at an arbitrary position and moves
all existing positions at or above it to the right. Inserting `false` at
position zero is therefore equivalent to a one-bit right shift.

The sample prints the set indexes after each transformation so the positional
effect of inserting both zero and one bits is visible.
