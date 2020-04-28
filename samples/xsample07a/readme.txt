Example on how to build (short) k-mers for DNA sequences,


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


Process FASTA file with a set of DNA sequences to build k-mers:
./xsample07a -fa all.fasta -kd pa_vect_coll.kd  -k 16 -t -j 4 


Load k-mer collection to memory:
./xsample07a  -kd pa_vect_coll.kd  -k 16 -t -j 4 

