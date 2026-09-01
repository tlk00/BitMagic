# BitMagic Library Samples

This directory contains small applications demonstrating BitMagic containers,
serialization APIs, search algorithms and larger end-to-end use cases.

## Bit-Vector Operations

- [`bvsample01`](./bvsample01) - basic operations to set/get bits with `bvector<>`, find cardinality (bit count).
- [`bvsample02`](./bvsample02) - set algebra operations, unions, intersections, equivalence of sets and lexicographical comparison with `find_first_mismatch()`.
- [`bvsample03`](./bvsample03) - use of different in-memory bitset compression options, calculate memory footprint.
- [`bvsample04`](./bvsample04) - serialization of `bvector<>` to save compressed BLOB to a file or a database.
- [`bvsample05`](./bvsample05) - use of `bvector<>::enumerator`, a fast iterator to decode a bit-vector into indexes of 1 bits.
- [`bvsample06`](./bvsample06) - allocator example.
- [`bvsample07`](./bvsample07) - logical operations between arrays and bit-vectors.
- [`bvsample08`](./bvsample08) - STL interoperability and set operations with iterators.
- [`bvsample09`](./bvsample09) - simple binary distance functions and pipeline for complex binary distance formulas.
- [`bvsample10`](./bvsample10) - extraction of a random subset for Monte Carlo simulations.
- [`bvsample11`](./bvsample11) - population counts for ranges in a bit-vector: `count_range()`, `count_to()` and `count_and()`.
- [`bvsample12`](./bvsample12) - review and comparison of different methods to set and clear bits.
- [`bvsample14`](./bvsample14) - serialization of `bvector<>`, logical operations on compressed BLOBs.
- [`bvsample15`](./bvsample15) - `bvector<>::find`, search for first and last set bit, dynamic range detection.
- [`bvsample16`](./bvsample16) - `bm::aggregator<>`, a utility class for fast logical operations on `bvector<>` groups. Tech note: <http://bitmagic.io/aggregator.html>
- [`bvsample17`](./bvsample17) - rank-select operations on `bvector<>` using `rs_index`.
- [`bvsample18`](./bvsample18) - `bvector<>::bulk_insert_iterator` for efficient bit-vector construction.
- [`bvsample18a`](./bvsample18a) - import of `bvector<>` from an external bit-stream with `bmbvimport.h`.
- [`bvsample19`](./bvsample19) - `bvector<>::merge`, merge of bit-vectors for partitioned processing.
- [`bvsample20`](./bvsample20) - `bvector<>::shift_right()` and `bvector<>::insert()`, bit-shifting and insertion.
- [`bvsample21`](./bvsample21) - `bvector<>::shift_left()` and `bvector<>::erase()`, bit-shifting and bit deletion.
- [`bvsample22`](./bvsample22) - algorithms for ranges, intervals (`bmintervals.h`) and `bvector<>::find_reverse()`.
- [`bvsample23`](./bvsample23) - `bm::interval_enumerator<>`, iterator for interval traversal of a bit-vector (`bmintervals.h`).
- [`bvsample24`](./bvsample24) - `bm::rank_range_split()` algorithm to split a bit-vector into equal rank/weight ranges.
- [`bvsample25`](./bvsample25) - `bm::visit_each_bit()`, `bm::visit_each_bit_range()`, `bm::for_each_bit()` and `bm::for_each_bit_range()` traversal algorithms.
- [`bvsample26`](./bvsample26) - immutable bit-vectors, construction, measuring memory savings, `bvector<>::freeze()`.
- [`bvsetalgebra`](./bvsetalgebra) - tutorial for Algebra of Sets operations. Tutorial: <http://bitmagic.io/set-algebra.html>
- [`bv3vlogic`](./bv3vlogic) - three-valued logic (Kleene) implemented on bit-vectors.
- [`bvsample01_64`](./bvsample01_64) - basic operations with 64-bit bit-vectors and initialization for 48-bit address space.

## Sparse Vectors

Sparse vector samples demonstrate bit-transposed vectors for scalar integer
types.

- [`svsample01`](./svsample01) - basic example of `sparse_vector<>`.
- [`svsample02`](./svsample02) - `sparse_vector<>` serialization, XOR compression and deserialization into read-only immutable succinct vector. Details: <http://bitmagic.io/bm-xor.html>
- [`svsample03`](./svsample03) - `sparse_vector<>` import, join and extract methods.
- [`svsample04`](./svsample04) - `sparse_vector<>` operations with nullable values.
- [`svsample05`](./svsample05) - set transformation algorithm to translate one set to another using a memory-efficient bit-transposed sparse vector. Tech note: <http://bitmagic.io/set2set-assoc-remap-opt.html>
- [`svsample06`](./svsample06) - load data and search for values in `sparse_vector<>`; demonstrates back-insert iterator, `const_iterator`, `bm::sparse_vector_scanner<>` and search/scan comparisons. Tech note: <http://bitmagic.io/sparse-vector-search.html>
- [`svsample07`](./svsample07) - insertion sort using `bm::sparse_vector_scanner<>::lower_bound()`.
- [`svsample08`](./svsample08) - `sparse_vector<>` deserialization to extract specific elements and ranges; shows bookmarks for faster range deserialization.
- [`svsample09`](./svsample09) - `bm::sparse_vector_find_first_mismatch`.
- [`svsample10`](./svsample10) - `bm::sparse_vector_scanner<>::find_gt()` searches: GT, GE, LT, LE and RANGE.

## Float Sparse Vectors

Float sparse-vector samples demonstrate bit-transposed vectors for floating
point values.

- [`svfsample01`](./svfsample01) - basic example of `sparse_vector_float<>`.
- [`svfsample02`](./svfsample02) - `sparse_vector_float<>` serialization, deserialization and const iterator.
- [`svfsample03`](./svfsample03) - `sparse_vector_float<>` comparison, join, merge and freeze.
- [`svfsample04`](./svfsample04) - `sparse_vector_float<>` extract, range extract, gather and back-insert iterator.
- [`svfsample05`](./svfsample05) - search `sparse_vector_float<>` with scanner, serialize result sets, and later retrieve selected values using deserialization-index assisted gather.

## RSC Sparse Vectors

RSC sparse-vector samples demonstrate rank-select compressed vectors for scalar
integer types.

- [`rscsample01`](./rscsample01) - basics of `rsc_sparse_vector<>`: load, unload and serialize.
- [`rscsample02`](./rscsample02) - back-insert iterator for `rsc_sparse_vector<>`; gather and range deserialization of a rank-select compressed container.
- [`rscsample03`](./rscsample03) - `rsc_sparse_vector<>::const_iterator`.
- [`rscsample04`](./rscsample04) - construction of `rsc_sparse_vector<>` with known NOT NULL elements and random `set()` / `inc()` in fast rank-select index mode.
- [`rscsample05`](./rscsample05) - serialization of a group of sparse vectors (data frame) using XOR compression. See also: <http://bitmagic.io/bm-xor.html>
- [`rscsample06`](./rscsample06) - `rsc_sparse_vector<>::gather()`, random or sorted-order extraction.
- [`rscsample07`](./rscsample07) - shared NULL plane for RSC/sparse-vector data-frame groups; serialization with skipped external NULL planes.

## String Sparse Vectors

String sparse-vector samples demonstrate bit-transposed vectors for ASCII(z)
strings.

- [`strsvsample01`](./strsvsample01) - basic `str_sparse_vector<>` operations: add values, optimize and iterate.
- [`strsvsample02`](./strsvsample02) - insertion sort using `str_sparse_vector<>` and `sparse_vector_scanner<>::lower_bound_str()`.
- [`strsvsample02a`](./strsvsample02a) - fast comparator for sorting succinct vector in compressed mode with `std::sort`.
- [`strsvsample03`](./strsvsample03) - `str_sparse_vector<>::back_insert_iterator`, remap compression and saving to disk using serialization API.
- [`strsvsample04`](./strsvsample04) - `str_sparse_vector<>` with NULL/unassigned values.
- [`strsvsample05`](./strsvsample05) - selective gather and range deserialization for `str_sparse_vector<>`; shows bookmarks for faster range deserialization.
- [`strsvsample06`](./strsvsample06) - `str_sparse_vector<>::const_iterator` in substring mode, `bm::sparse_vector_scanner<>` string search and result-set iteration.
- [`strsvsample07`](./strsvsample07) - `bm::sparse_vector_scanner<>::pipeline` for thousands of string searches, statistics and succinct-vector remapping.
- [`strsvsample08`](./strsvsample08) - `bm::sparse_vector_scanner<>::bfind_eq_str()` binary search with binding, `remap()`, `optimize()` and `freeze()`.
- [`strsvsample10`](./strsvsample10) - mmap-backed gather deserialization from a serialized `str_sparse_vector<>` BLOB using bookmarks and deserialization index-assisted retrieval.

## Complex Examples

These examples combine multiple BitMagic techniques and demonstrate larger use
cases.

- [`xsample01`](./xsample01) - advanced methods for handling super sparse sets and benchmarks of set operations. Details: <http://bitmagic.io/case-ER-join.html>
- [`xsample02`](./xsample02) - sparse vector based counting sort and histogram construction techniques compared to `std::sort()`. Details: <http://bitmagic.io/hist-sort.html>
- [`xsample03`](./xsample03) - search in human genome data using `bm::sparse_vector<>` and `bm::rsc_sparse_vector<>`. Details: <http://bitmagic.io/succinct-snp-search.html>
- [`xsample04`](./xsample04) - fingerprint pattern matching with SHIFT-AND Bitap algorithm; DNA substring search. Details: <http://bitmagic.io/dna-search.html>
- [`xsample04a`](./xsample04a) - fast construction of bit-indexes using bulk insert iterator, multi-threaded partitioning and merge. Details: <http://bitmagic.io/dna-search-idx.html>
- [`xsample05`](./xsample05) - memory compressed dictionary search using a part of the NED catalog of celestial objects. Details: <http://bitmagic.io/star-search.html>
- [`xsample06`](./xsample06) - DNA 2-bit per bp compression using `sparse_vector<>` and `bm::sparse_vector_find_first_mismatch()`. Details: <http://bitmagic.io/dna-compare.html>
- [`xsample07`](./xsample07) - DNA k-mer counting for short k-mers; demonstrates `bvector<>`, succinct `rsc_sparse_vector<>`, search scanner and multi-threaded map-reduce style techniques.
- [`xsample08`](./xsample08) - data model for Genomics Viewer using bit-intervals and succinct vectors, with ASCII rendering. Details: <http://bitmagic.io/gen-layout.html>
- [`xsample09`](./xsample09) - histogram construction using memory compact vectors. Details: <http://bitmagic.io/bm-hist-compr.html>
- [`xsample10`](./xsample10) - compressive scrolling for Model View Controller style data sets; demonstrates `bm::sparse_vector_serializer<>::set_bookmarks()` and range deserialization. Details: <http://bitmagic.io/bm-mvc.html>
