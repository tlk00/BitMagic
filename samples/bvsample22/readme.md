# bvsample22: Ranges and intervals

This example uses an RLE-friendly `bm::bvector<>` as a compact collection of
closed integer intervals.

It demonstrates `set_range()`, `is_all_one_range()`, `any_range()`,
`bm::is_interval()`, reverse search and interval-boundary discovery with
`bm::find_interval_start()` and `bm::find_interval_end()`. It also edits and
extracts ranges with `clear_range()`, `keep_range()` and `copy_range()`.

The serialization section enables block bookmarks and restores only a closed
range with `bm::deserialize_range()`. Bookmarks slightly increase BLOB size
but can make repeated range extraction faster.

## Related serialization examples

- [bvsample04](../bvsample04/readme.md): RAM serialization, buffer ownership, and ordinary deserialization.
- [bvsample14](../bvsample14/readme.md): set algebra and count operations on serialized RAM BLOBs.
- [bvsample27](../bvsample27/readme.md): buffered file output, caller-owned finish(), and RAM compatibility checks.

The RAM, operation, and range deserializers consume an in-memory BLOB.
The file-output example preserves that format; it does not introduce a file
deserializer or make those readers accept a C++ stream.
