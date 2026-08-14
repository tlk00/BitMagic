# svfsample04: Float extraction and back insertion

This example extracts values from `bm::sparse_vector_float<>` into ordinary
float arrays and demonstrates its back insert iterator.

It validates three extraction patterns:

- `decode()`/`extract()` for a contiguous run;
- `extract_range()` for a smaller slice copied to the start of an output
  buffer; and
- `gather()` for arbitrary logical indexes with a sort-order hint.

The second demo appends values through both `add()` and assignment. It also
shows copy and move construction of back inserters; moving disconnects the
original iterator from the destination.
