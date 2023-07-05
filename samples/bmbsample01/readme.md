# basic_bmatrix<> serailization example

This example discusses serialization techniques for succinct vectors.
basic_bmatrix<> used to be an internal class used for compressive structures 
construction. It turns out it can useful by its own. 

Example shows how to:

1. Simple serialization
2. Deserialization as immutable (read-only) structure
3. Deserialization by range
4. Deserialization using AND mask

All the techniques are fully applicable for other sparse vectors in BM library.
