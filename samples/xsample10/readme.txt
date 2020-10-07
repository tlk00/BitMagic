Density histogram construction using memory compact vectors.

http://bitmagic.io/bm-hist-compr.html


Example uses simulated data imitation human variations on chr1 (size wise).
Rank-Select compressed vector used to keep the test data set.
Then we run a few algorithms to construct density histograms as bit-transposed 
memory compact vectors or two types (SV and RSC).

After construction example serializes both containers to disk to understand 
compressed footprint.

Example runs benchmarks to illustrate random access performance to compact 
histogram data.


How to build:

Use ./build_all.sh script (it is x86 centric as builds SIMD versions).

How to run:

./run_all.sh



Expected output sample:
-------------------------


Number of elements in the data set: 29411768
Histogram 1 (SV)
  size: 100001
  RAM size: 119104
  file size: 71727
Histogram 2 (RSC)
  size: 250000001
  NOT NULL count: 100001
  RAM size: 1157552
  file size: 235580
Histogram 3 adaptive (RSC)
  size: 249999385
  NOT NULL count: 305
  RAM size: 168400
  file size: 1052

Access sample size = 1961594

01. Histogram 1 construction (SV)) ; 16.2047 ms
02. Histogram 2 construction (RSC) ; 410.607 ms
03. Histogram 3 (adaptive) construction (RSC) ; 24.0972 ms
04. Verification; 18.4433 ms
05. Serialize and save the histogram 1 (SV); 4.19348 ms
06. Serialize and save histogram 2(RSC); 5.69941 ms
07. Serialize and save the adaptive histogram 3 (RSC); 0.855228 ms
08. Random access test H1 (SV); 82.6797 ms
09. Random access test H2 (RSC); 110.062 ms
10. Random access test H3 adaptive (RSC); 43.9542 ms

