Example on how to build (short) k-mers for DNA sequences,
translate into a k-mer presence-absence bit-vector
and run k-mer counting algorithm using DNA binary fingerprints (Bitap algorithm).


How to get test data:
---------------------

>wget bitmagic.io/data/NC_000913.3.fa.gz
>wget bitmagic.io/data/NC_000001.11.fa.gz



How to build:
--------------

1. apply environment variables at the BitMagic project root:

>. ./bmenv.sh
or
>source ./bmenv.sh


2. Build regular version:
>make rebuild 

or AVX2

make BMOPTFLAGS=-DBMAVX2OPT rebuild


OR just use

>./build_all.sh



How to run:
------------

>./xsample07_avx2 -kd test.kd -fa NC_000001.11.fa -k 10  -t -j 4