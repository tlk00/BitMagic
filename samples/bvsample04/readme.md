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
