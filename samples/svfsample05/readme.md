# svfsample05: Search result sets and indexed gather deserialization

This example demonstrates a two-stage workflow for `bm::sparse_vector_float<>`
data. It searches a floating-point sparse vector for simulated anomalies, saves
the search hits as compact bit-vector BLOBs, serializes the float vector itself,
serializes a deserialization index for that vector, and later restores only the
selected values using gather deserialization assisted by the restored index.

The sample is intentionally not a benchmark. It is meant to make the API pattern
clear for applications where the full vector is expensive to keep resident, but
only a small subset of values is needed later.

## Workflow

Stage 1 builds a synthetic signal. Most values are small baseline values, while
periodic positive and negative spikes simulate events or anomalies.

The sample prints a plain `std::vector<float>` payload estimate first. This is
the simple uncompressed baseline: one 32-bit float per logical element, before
considering allocator overhead or any additional application metadata.

The scanner searches the vector and produces `bm::bvector<>` result sets:

- positive spikes, found with `bm::sparse_vector_scanner<>::find_gt_float()`;
- negative spikes, found with `bm::sparse_vector_scanner<>::find_lt_float()`;
- all anomalies, created by combining the positive and negative result sets.

The result sets are serialized independently from the float vector. In a real
system these BLOBs can be saved as query results, filters, segments, alert sets,
or other compact address lists.

Stage 1 also serializes the float sparse vector with bookmarks enabled.
Bookmarks add a small amount of metadata to the serialized sparse-vector planes.
The extra metadata gives the later deserializer places where it can skip ahead
instead of walking the whole compressed stream.

After the vector BLOB is available, Stage 1 builds a deserialization index over
it and serializes that index into a RAM buffer. In a real application the index
BLOB can be stored next to the vector BLOB, so startup does not need to rebuild
the index from the full serialized vector.

Stage 2 represents a later retrieval phase. The source sparse vector is gone, so
the program keeps only:

- the serialized float-vector BLOB;
- the serialized result-set BLOBs;
- the serialized deserialization-index BLOB.

The sample restores the deserialization index from its BLOB, restores each
result set, and uses the result set as the gather mask. The deserializer then
reconstructs only the values addressed by the mask.

## Deserialization index

The deserialization index is a resident helper object built from the serialized
float-vector BLOB. It is not a copy of the sparse vector. It is closer to a map
of the compressed stream: it records where useful regions and bookmark-assisted
jump points are located.

The sample serializes the index with
`bm::sparse_vector_float_deserialization_index_serializer<>` and restores it
with `bm::sparse_vector_float_deserialization_index_deserializer<>`. The sample
also checks the RAM round trip with `equal()` before using the restored index.

For repeated sparse retrievals, this reduces both disk or memory traffic and CPU
spent decoding unrelated compressed blocks. The tradeoff is that the index uses
some RAM and the selected values still need to be decompressed when gathered.

This pattern is useful when an application chooses to trade CPU for memory:
store the large vector in compressed serialized form, release the in-memory
vector, and periodically retrieve only a small fraction of elements. More
advanced systems can add prefetching or multi-threaded scheduling around the
same idea, but this sample keeps the flow single-threaded and direct.

## Result processing

After gather deserialization, the sample computes simple aggregate statistics
over the restored values:

- number of selected values;
- sum;
- average;
- minimum and maximum.

The aggregation is deliberately simple. Its purpose is to show that once the
selected values are restored into a sparse vector, ordinary application logic can
work with them without loading the full source vector back into memory.

The code uses the filtered overload of `bm::for_each_sparse()` for this step.
The restored result-set bit-vector is passed as the filter, so the visitor sees
only the values selected by the original scanner search.

## Expected output

A typical run prints output similar to this:

```text
Stage 1: build float sparse vector
  vector size = 524288 values
  plain std::vector<float> payload estimate = 2.000 MB
Stage 1: search for synthetic anomalies
  positive spikes: 127 positions, serialized result-set BLOB = 16 bytes
  negative spikes: 85 positions, serialized result-set BLOB = 16 bytes
  all anomalies: 212 positions, serialized result-set BLOB = 352 bytes
Stage 1: serialize float sparse vector with bookmarks
  serialized float-vector BLOB = 1583768 bytes
Stage 1: build and serialize deserialization index
  deserialization index memory = 12512 bytes
  serialized deserialization index BLOB = 178 bytes

Stage 2: restore serialized deserialization index
Stage 2: gather values from serialized BLOB
  positive spikes: count=127, sum=2191.400, avg=17.255, min=15.710, max=18.965
  negative spikes: count=85, sum=-1617.110, avg=-19.025, min=-20.745, max=-17.370
  all anomalies: count=212, sum=574.290, avg=2.709, min=-20.745, max=18.965
```

The first stage shows the size of the logical signal and the compact serialized
forms produced from it. The positive and negative result-set BLOBs are very
small because the synthetic anomalies are periodic and compress well as
bit-vectors. The combined anomaly set is larger because it contains both
patterns.

The second stage shows the later retrieval flow. The source vector is no longer
used. The deserialization index is restored from its serialized BLOB, result-set
BLOBs are restored into bit-vectors, and gather deserialization restores only
the selected positions. The aggregate values are computed from those restored
positions.
