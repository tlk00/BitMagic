# strsvsample07

Example for bit-transposed succinct vector search using `bm::sparse_vector_scanner<>::pipeline`.

BitMagic Succinct string container `bm::str_sparse_vector<>` uses bit-slicing as a method 
to save memory, but it also provides an array of efficient search algorithms based on
bit-vector logical operations.


Example illustrates use of a bm::sparse_vector_scanner<>::pipeline to do a quick bulk search in the unordered succinct sparse vector.

Given a vector of values (column 1) we search for two search values: "str1" and "str4".
For each search value scanner produces a result, conside it a filter for the vector (column 0)
where 1s correspond to positions in the vector matching the search request.


| vector      | "str1"  | "str4" |  UNION ALL  |
| ----------- | -------:|------:|-------------:|
| "str1"      | 1       | 0     | 1 
| "str2"      | 0       | 0     | 0
| "str1"      | 1       | 0     | 1
| "str4"      | 0       | 1     | 1
| "str2"      | 0       | 0     | 0
| "str3"      | 0       | 0     | 0
|             |         |       |
| COUNT:      | 2       | 1     | 3


Bulk search can be used to construct search bit-vectors (inverted lists) or compute 
population counts (can be used for histogram construction) a lot faster than it can be done
using looped single item search. 

Scanner offers options not to produce the materialize the search results as bit-vectors, 
but rather compute population count (2 and 1 in the example above). 

Another option is to merge all the result together as OR (SET UNION).


## Build

Build this example as a part of BitMagic library or use `./build_all.sh` to prepare benchmarks for regular, SSE4.2 and AVX2. For Arm consider `build_all_arm.sh` for regular and NEON accelerated enabled builds (tested with R-Pi).

## Benchmarks

BitMagic applies a range of optimization techniques: SIMD, cache L1/L2 reuse, 
multi-way memory reads (bandwidth optimizations), algorithmic optimizations to reuse uncompressed GAP blocks.

This example showcases several search methods to benchmark different search methods.


## Sample run

Use ./run_all.sh script to run the benchmarks with expected output similar to:

 
> Generating the test data... OK
> 
> Remapping the data to create compressed vector OK
> 
> Running benchmark tests..
> 
> PASS = 0 -- remap/optimized
> 
> scanner<>::find_eq_str(); 11.04 sec
> 
> scanner::pipeline find_eq_str(); 6.028 sec
> 
> scanner::pipeline+MASK find_eq_str(); 1.817 sec
> 
> scanner::pipeline find_eq_str()-count(); 7.017 sec
> 
> scanner::pipeline find_eq_str()-OR(); 7.904 sec
> 
> PASS = 1 -- NOT remapped
> 
> scanner<>::find_eq_str(); 10.02 sec
> 
> scanner::pipeline find_eq_str(); 5.114 sec
> 
> scanner::pipeline+MASK find_eq_str(); 1.942 sec
> 
> scanner::pipeline find_eq_str()-count(); 7.563 sec
> 
> scanner::pipeline find_eq_str()-OR(); 7.3 sec
>
>...


 
