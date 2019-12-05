Example on how to construct sparse vector lexicographical compare functions

Sparse vector in this example is used to keep string of DNA sequences using
compression close to 2-bits per string character (base pair).

This example show how to do mismatch search and construct comparison function
on succinct sparse vector. 

Technical notes:
http://bitmagic.io/dna-compare.html



How to build:
--------------

1. apply environment variables at the BitMagic project root:

>. ./bmenv.sh
or
>source ./bmenv.sh


2. Build regular version:
make rebuild 

or AVX2

make BMOPTFLAGS=-DBMAVX2OPT rebuild


How to run:
------------

>./xsample06

size() = 20 : CTTGGAANNNNNNGCCCTAA
Generate test data.
generated mismatches=15271
SV mismatch search test
STL interator mismatch test
::strncmp test
SV compare test

Performance:
1. SV mismatch; 33.34 sec
2. STL iterator; 753.7 sec
3. strcmp() test ; 183.6 sec
4. sv compare; 34.52 sec
