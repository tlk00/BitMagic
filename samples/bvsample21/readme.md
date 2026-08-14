# bvsample21: Left shifts and erasure

This example shows sequence-style deletion operations on an optimized
bit-vector.

`shift_left()` moves set bits one position toward zero and reports a carry for
the bit shifted out. `erase(position)` removes a position and shifts every
higher position left; erasing position zero is therefore equivalent to a
one-bit left shift.

Repeated shift, insert or erase operations can make the internal block layout
less compact. The sample calls `optimize()` after editing to restore memory
compression.
