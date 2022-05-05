# strsvsample08

Example for string succinct vector binary search using `bm::sparse_vector_scanner<>`

For the best sorted search performance scanner class can be 
pre-created with the "sampling factor" parameter and binded to a particular 
succinct string vector. Binded search would create an approximate uncompressed 
index (it comes with a certain memory penalty, but not much).

Samplig factors can be: 4, 8, 16 (default), 32 and 64. 
Lowes sampling is 4 (offers minimal memory penalty) and 64 - higher penalty 
but faster speed. 


strsvsample08 provides some benchmarking. It is based on synthetic data, 
so it is a good idea to test on your real data.

Another important use case here is use of remap()-optimize()-freeze()
to turn succinct vector into the most compact memory mode.
 
