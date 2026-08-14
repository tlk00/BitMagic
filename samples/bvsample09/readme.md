# bvsample09: Binary distance metrics

This example computes similarity and distance values between two bit-vectors
without materializing intermediate result sets.

It uses `bm::count_xor()` for Hamming distance and combines
`bm::count_and()` with the source population counts to compute the Dice
coefficient. The same Dice calculation is then repeated with a
`bm::distance_metric_descriptor` pipeline and `bm::distance_operation()` so
the intersection and both source counts are collected in one pass.

The pipeline form is useful when a formula needs several related counts for
the same pair of vectors.
