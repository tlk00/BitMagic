Example to illustrate dataset scrolling using different de-serialization
techniques with succinct vectors.

Example to illustrate how to use a sparse_vector_float in a real-world
setting

https://bitmagic.io/bm-9.1.0


How to build and run (x86 centric):

./build_all.sh

./run_all.sh


How to run:
-------------
./xsample11


Application notes:
------------------

Example reads the data in the EURUSD_H1.csv and the USDJPY_H1.csv 
The data in the .csv's is the exchange rate information between EUR to USD and 
USD to JPY for each hour of every day from 
June 24, 2010 5:00am to July 7, 2026 3:00pm, (15:00)
This information includes each hours open, close, high, low and volume.

This program then calculates each hour's percent change based on the close, and
stores all of the information in a str_sparse_vector for times, sparse_vector
for volume, and sparse_vector_float's for the exchange rates and percent changes,
Optimizes the sparse_vectors, and freezes them as no new information will be added.

Information on the memory used of the sparse_vector_float



Example generates a sample data frame, serializes it into memory using compressive
techniques and setting bookmarks to support range deserialization.

The use case illustrated is the the main body of data remains compressed in RAM
we create a range specific data frame.

Benchmarks simulate forward scrolling, every new range overlaps with a previous 
range (assumed that application allows relative scrolling within ranges).
Under the assumption of overlap between ranges it demonstrates various 
range deserialization techniques on how to preserve overlap and only deserialize
the new range.

Sample uses only forward scrolling, more complex application should handle random
ranges with or without overlap. 


