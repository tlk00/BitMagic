# bvsample18: Bulk insert iterator

This example contrasts `bvector<>::insert_iterator` with
`bvector<>::bulk_insert_iterator`.

The regular inserter changes the destination immediately. The bulk inserter
buffers indexes so it can submit them more efficiently, especially for sorted
data; consequently, the destination count may not change while assignments
are still being buffered.

Buffered values are committed when the iterator is destroyed or when
`flush()` is called explicitly. The explicit form makes the commit point clear
and ensures any exception is reported from normal code rather than a
destructor.
