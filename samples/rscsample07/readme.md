# rscsample07: shared NULL plane for RSC data frames

This example demonstrates a data-frame style construction where several sparse
containers have identical NULL/not-NULL shape and therefore share one NOT NULL
bit-vector.

The leading column is an `rsc_sparse_vector<>`. It owns the NOT NULL plane and
is responsible for its lifetime. Secondary columns attach to that plane after
construction:

- another `rsc_sparse_vector<>`, which can also reuse the leader's rank-select
  index after the leader is synchronized;
- a plain `sparse_vector<>`, which shares the NULL plane but keeps its own value
  bit-planes.

## Construction pattern

1. Construct all participating vectors as NULL-able containers.
2. Load or build all vectors with the same logical NULL shape.
3. Call `sync()` on the RSC owner when rank-select access is needed.
4. Attach secondary vectors to the owner's NOT NULL plane.
5. Keep the owner alive longer than all attached vectors.

Secondary vectors do not own the shared NOT NULL plane. Destruction, clear,
assignment and optimization of a secondary vector must not delete the master
plane.

## Mutations

A write through any attached vector that changes NULL/not-NULL state changes the
shared master plane. After such a shape-changing write, RSC vectors using that
plane need their rank-select state invalidated and rebuilt before synced access.
Value-only writes to already NOT NULL positions do not require rebuilding NULL
state.

## Serialization

The master vector is serialized as a self-contained object because it owns the
NOT NULL plane. Attached follower vectors skip externally owned NULL planes by
default, reducing serialized size.

Deserialization must follow the same build pattern:

1. Deserialize the master first.
2. Attach all follower vectors to the restored master's NOT NULL plane.
3. Deserialize follower vectors whose external NULL plane was skipped.

Loading a skipped-NULL follower into a standalone, non-attached vector is a
configuration error and is expected to throw an exception.
