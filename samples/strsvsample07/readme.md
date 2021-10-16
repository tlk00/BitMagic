# strsvsample07

Example for bit-transposed succinct vector search using bm::scanner<>::pipeline.

BitMagic Succinct string container bm::str_sparse_vector<> uses bit-slicing as a method 
to save memory, but it also provides an array of efficient search algorithms based on
bit-vector logical operations.


Example illustrates use of a bm::scanner<>::pipeline to do a quick bulk search in the unordered 
succinct sparse vector.

Bulk search can be used to construct search bit-vectors (inverted lists) or compute 
population counts (can be used for histogram construction) a lot faster than it can be done
using looped single item search.

BitMagic applies a range of optimization techniques: SIMD, cache L1/L2 reuse, 
multi-way memory reads (bandwidth optimizations), algorithmic optimizations to reuse uncompressed GAP blocks.

This example showcases several search methods and benchmarks various search methods.


## Build

Build this example as a part of BitMagic library or use ./build_all.sh to prepare benchmarks for SSE4.2 and AVX2.


## Sample run

Use ./run_all.sh script to run the benchmarks:

 
> 
> Regular build
> Generating the test data... OK
> 
> Remapping the data to create compressed vector OK
> 
> Running benchmark tests..
> 
> PASS = 0 -- remap/optimized
> 
> scanner<>::find_eq_str(); 12.9 sec
> 
> scanner::pipeline find_eq_str(); 4.755 sec
> 
> scanner::pipeline find_eq_str()-count(); 4.474 sec
> 
>...


 
