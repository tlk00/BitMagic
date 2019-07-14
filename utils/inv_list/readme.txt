Test utility for benchmarking compression of test collections.


Input format as:

Document identifier data set
Packaged by D. Lemire on April 3rd 2014
Based on data sets prepared by L. Boytsov

The
main data files are in a binary flat file format
described simply as follows. They contain a sequence 
of 32-bit unsigned integers (written in little endian
mode). These integers make up lists of integers. Each list 
starts with an integer indicating its length followed by 
the corresponding number of integers in sorted order. 


HOW TO BUILD (in-place make):

1. Export root path for the build (BitMagic root dir)
run:
> . ./bmenv.sh  ( from the root directory  dot is important!)
OR just
>PROJECT_DIR=`pwd`; export PROJECT_DIR

2. build using:

make BMOPTFLAGS=-DBMAVX2OPT
OR
make BMOPTFLAGS=-DBMSSE42OPT
OR 
make

2. HOW TO RUN


To run this test you need to download or prepare benchmark files (integer sets)


2.1 Read the input set, create the output (bit-vector compressed format):

./bminv -c -u32in /gov2.sorted -bvout /gov2.sorted.bv 

2.2 Run validation comparison to make sure there are no errors
./bminv -verify -u32in /gov2.sorted -bvin /gov2.sorted.bv 


Help:
./bminv -h









