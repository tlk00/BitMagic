Example on how to build (short) k-mers for DNA sequences,
translate into a k-mer presence-absence bit-vector and run k-mer counting algorithms.
K-mer counts (Term Frequency vector) can be saved to disk as succinct/compressed BLOB.

Based on TF vector example can compute top N percent of highly represented k-mers
which are most likely parts of repeats



How to get test data:
---------------------

>wget bitmagic.io/data/NC_000913.3.fa.gz
>wget bitmagic.io/data/NC_000001.11.fa.gz



How to build:
--------------

0. use project make file

OR 

1. apply environment variables at the BitMagic project root:

>. ./bmenv.sh
or
>source ./bmenv.sh


2. Build regular version:
>make rebuild 

or AVX2

make BMOPTFLAGS=-DBMAVX2OPT rebuild


OR just use (it will use GCC to create all build variants)

>./build_all.sh



How to run:
------------

Help:
>./xsample07 -h

Generate k-mer fingerprint (4 threads):
>./xsample07_avx2 -kd test10.kd -fa NC_000001.11.fa -k 10  -t -j 4

Generate k-mer fingerprint with diagnostics checks(slower)
./xsample07 -kd test10.kd -fa NC_000001.11.fa -k 10 -diag

Generate k-mer fingerprint and count all k-mers (8threads):
>./xsample07_avx2 -kd test10.kd -kdc test10.kdc -fa NC_000001.11.fa -k 10  -t -j 8

Generate k-mer fingerprint and count all k-mers and compute frequent k-mer vector for top 10% of all k-mers.
k-mer frequency histogram is reported to a file (hmap.tsv) 
build and save the k-mer fingerprint cleaned from over-represented k-mers (test.kdc)
:
>./xsample07_avx2 -kd test.kd -kdf test.kdf -kdc test.kdc  -fa NC_000001.11.fa -k 16  -t -j 4 -kh hmap.tsv -fpc 10