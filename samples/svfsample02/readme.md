# svfsample02: Float serialization and iteration

This example covers serialization, selective restoration and traversal of
`bm::sparse_vector_float<>`.

The serialization section uses `sparse_vector_float_serial_layout` and shows
its size, capacity, buffer and memory-release interfaces. It performs full
deserialization, closed-range deserialization and masked deserialization, then
checks the restored values.

The iterator section demonstrates construction through `begin()`, value and
position access, comparisons, `advance()`, increment, `go_to()` and
`bm::for_each_sparse()`. A serialized float vector consists of a header plus
serialized sign, exponent and mantissa components.
