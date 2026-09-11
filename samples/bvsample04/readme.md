# bvsample04: Bit-vector serialization

This example serializes compressed `bm::bvector<>` objects into memory and
restores them with `bm::deserialize()`.

It demonstrates:

- optimizing a vector before serialization and using `calc_stat()` to size a
  raw output buffer;
- reusing `bm::serializer<>` and configuring compact output;
- using the serializer's RAII-managed `buffer` type;
- accumulating multiple serialized vectors into one destination, where
  repeated deserialization performs a logical OR; and
- `optimize_serialize_destroy()` for destructive serialization when the
  source vector is no longer needed.

The program validates the restored result against an explicit OR of the two
source vectors.

## Related serialization examples

- [bvsample14](../bvsample14/readme.md): set algebra and count operations on serialized RAM BLOBs.
- [bvsample22](../bvsample22/readme.md): bookmarks and selective range deserialization.
- [bvsample27](../bvsample27/readme.md): buffered file output, caller-owned finish(), and RAM compatibility checks.

The RAM, operation, and range deserializers consume an in-memory BLOB.
The file-output example preserves that format; it does not introduce a file
deserializer or make those readers accept a C++ stream.
