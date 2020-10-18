Example to illustrate dataset scrolling using different de-serialization
techniques with succinct vectors.

How to build and run (x86 centric):

./build_all.sh

./run_all.sh



Application notes:

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


