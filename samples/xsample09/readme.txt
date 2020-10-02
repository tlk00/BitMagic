Density histogram construction using sparse vectors.

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



