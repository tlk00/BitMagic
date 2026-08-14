# bvsample10: Random subsets

This example builds a bit-vector with values in two widely separated ranges
and uses `bm::random_subset<>` to select a requested number of its set bits.

A single sampler instance is reused to create two independent 20-element
subsets. Each result is another `bm::bvector<>` whose set indexes come only
from the source vector.

This pattern is useful for random sampling and Monte Carlo workflows where a
compact set must be sampled without expanding it into a separate array first.
