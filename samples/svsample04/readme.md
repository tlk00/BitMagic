#svsample04

This example illustrates succinct vector manipulation methods to set/get and clear elements.

Conventinal random access methods work ok, but succinct container also provids methods for bulk extraction and bulk settting to NULL or clearing elements.

Thise operations are implented as logical ops on bit-slices, SIMD accelarated and ayutomatically detect some possibilities to release sparse vector blocks.



