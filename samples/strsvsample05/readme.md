# strsvsample05: Selective string-vector deserialization

This example serializes a large remapped `bm::str_sparse_vector<>` and restores
only a selected part of the compressed BLOB.

The serializer enables block bookmarks, which add a small amount of metadata
so the deserializer can skip unrelated regions more efficiently. The same
closed interval is restored first with a bit-vector mask and then with
`deserialize_range()`, the clearer and faster form for one contiguous range.

The program compares every restored string with the source and reports a
successful gather check. This pattern supports paging or projection without
fully materializing a large succinct vector.
