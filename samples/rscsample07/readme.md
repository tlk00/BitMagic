# rscsample07: RSC shared NULL plane

This sample demonstrates a data-frame style group of BitMagic sparse vectors where one `rsc_sparse_vector<>` owns the NOT NULL plane and secondary vectors attach to it.

The example covers:

- building multiple columns with the same NULL/not-NULL shape;
- attaching RSC and plain sparse-vector followers to the master's NULL plane;
- serializing followers without the external NULL plane;
- restoring a connected group in the required order: master first, attach followers, then load followers;
- using `set_bookmarks(true, 64)` for faster range deserialization;
- range restore where the same `[from, to]` interval is applied to the master and every attached follower;
- verifying that attached vectors keep the shared NULL plane and that RSC followers can reattach to the master's rank-select index after sync.

For connected groups, construct or attach all columns before deserializing followers whose external NULL plane was skipped. During range restore, use the same range for all columns so the restored value planes match the master-owned NULL plane.
