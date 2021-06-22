# Three-valued logic

BitMagic implements 3-VL set algebra (Kleene) functions ("bm3vl.h") using compact (and fast) bit-vectors.

[https://en.wikipedia.org/wiki/Three-valued_logic
](https://en.wikipedia.org/wiki/Three-valued_logic)

This example shows how to initialize 3VL vectors and perform logical operations AND, OR, NOT.


# Implementation details

BitMagic implements three-valued logic using two bit-vectors:

1. Bit-vector of values (known 1s|0s) 1s here represent known TRUE, 0s - known false OR UNKNOWNs 
2. Bit-vector of knowns (assigned values) or NULL vector contains 1s in positions where value was explicitly assigned (known true or known false), 0s represent NULL (NIL, N/A) values.

BitMagic allows direct access to both bit-vectors, with a requirement not to assign 1 (true) value for UNKNOWN (NOT NULL) elements, because UNKNOWN true does not make sense and underlying logical algorithms does not take this into account and will produce incorrect results if it happens.

For example:
------------
	1001000000000 - values bit-vector
	1101000000000 - NULL vector

	idx=0: true  (1)
	idx=1: false (-1)
	idx=2: NULL  (0)
	idx=3: true  (1)
	idx=4+: NULL


# Logical operations

3VL algebra follows Kleene (and Priest) truth tables.

More on that:
[https://en.wikipedia.org/wiki/Three-valued_logic
](https://en.wikipedia.org/wiki/Three-valued_logic)
