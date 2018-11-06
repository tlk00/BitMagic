Example demonstrates fast fingerprint/index construction based on DNA sequence.

1. Download a FASTA file. 
We recommend to use human chromosome 1, since it is big enough (250 million base pairs, or approximately 8% of human genome).

2. Build
For a regular in-place incidental Linux or MacOS build you can use 
./build_all.sh 
This will compile regular, SSE4.2 and AVX2 versions

3. Run benchmarking cases:

4-threaded version
./xsample04a -fa NC_000001.11.fa -j 4 -t

8-threaded AVX2 build
./xsample04a_avx2 -fa NC_000001.11.fa -j 8 -t


