# bvsample02: Set algebra and comparison

This example treats `bm::bvector<>` objects as sets of integer indexes. It
demonstrates both expression and in-place forms of:

- intersection (`&` and `&=`);
- union (`|` and `|=`);
- set difference (`-` and `-=`); and
- symmetric difference (`^` and `^=`).

The XOR result is used as a mismatch vector. The sample then shows how
`find_first_mismatch()` can locate the first difference without materializing
the complete XOR product. It finishes with equality and lexicographical
comparison using `operator==` and `compare()`.
