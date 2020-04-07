Example on how to build (short) k-mers for DNA sequences,
translate into a k-mer presence-absence bit-vector
and run k-mer counting algorithm using DNA binary fingerprints (Bitap algorithm).


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

Generate k-mer fingerprint and count all k-mers (8threads):
>./xsample07_avx2 -kd test10.kd -kdc test10.kdc -fa NC_000001.11.fa -k 10  -t -j 8