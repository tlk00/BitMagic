Example to illustrate dataset scrolling using different de-serialization
techniques with succinct vectors.

http://bitmagic.io/bm-mvc.html



How to build and run (x86 centric):

./build_all.sh

./run_all.sh


How to run:
-------------
./xsample10_sse42 -bookm 16 -check

  Serialization bookmark at:16
  ID BLOB size = 125490402
  SEQ BLOB size = 20354073
  POS BLOB size = 9662368
Total = 155506843

01. Scrolling test range; 19.66 sec
02. Scrolling test clear/range/merge; 13.05 sec
03. Scrolling merge/keep_range; 13.25 sec




Application notes:
------------------

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


