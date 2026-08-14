# svsample08: Selective sparse-vector deserialization

This example serializes a NULL-able `bm::sparse_vector<>` and restores selected
logical positions from the compressed BLOB.

A buffered back inserter creates the sample vector, including explicit NULLs.
The serializer enables block bookmarks, which help the deserializer skip
unneeded regions. The program then compares mask-based `deserialize()` with
`deserialize_range()` for one closed interval.

The destination retains the source's logical size and NULL positions, while
unselected assigned values read as zero. Range deserialization is the clearer
and faster choice when the selection is contiguous.
