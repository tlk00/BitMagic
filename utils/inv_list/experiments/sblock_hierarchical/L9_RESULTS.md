# Full Gov2 level-9 wave selector comparison

Run: `/Volumes/USB_1TB/BitMagic-experiments/gov2-l9-wave-policy-20261114`

All six streams passed full raw-ID verification and independent decode of
1,107,205 vectors and 5,979,715,441 IDs per dataset. The same frozen `bminv`
binary and ten compression threads were used throughout. Each mode was run
once; timings are observations, not controlled performance estimates. The
exact mode also collected per-leaf features, adding instrumentation cost.

| Dataset | Mode | Stream bytes | Bits/ID | Encode s | Verify s | Decode s | Bytes above exact |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| sorted | 108 preset | 2,642,908,243 | 3.535831 | 1,073.558 | 43.413 | 23.701 | 276,888 |
| sorted | exact 108/109 | 2,642,631,355 | 3.535461 | 1,170.764 | 46.826 | 23.970 | 0 |
| sorted | tree 108/109 | 2,642,794,003 | 3.535679 | 1,043.161 | 46.762 | 25.697 | 162,648 |
| unsorted | 108 preset | 5,331,827,332 | 7.133219 | 1,394.869 | 81.017 | 63.755 | 4,648 |
| unsorted | exact 108/109 | 5,331,822,684 | 7.133213 | 1,612.919 | 83.261 | 65.717 | 0 |
| unsorted | tree 108/109 | 5,331,824,468 | 7.133215 | 1,436.840 | 80.300 | 63.659 | 1,784 |

The L9 tree's emitted size regret matches the leaf-stage estimate exactly:
162,648 bytes sorted and 1,784 bytes unsorted. Sorted leaf 6 accounts for
125,312 bytes of regret, including one 2,768-byte miss; the other leaves
together account for 37,336 bytes. Unsorted leaf 6 accounts for 1,312 of
1,784 bytes, with a maximum miss of four bytes in every leaf.

The tree saves 114,240 bytes against the sorted 108 preset and 2,864 bytes
against the unsorted preset. Exact selection saves another 162,648 and 1,784
bytes respectively. These gains are small compared with the full streams.
The single-pass encode result favors the tree over exact, but the exact run's
feature-export work prevents a clean measurement of selector cost. The
unsorted tree took 42 seconds longer than the preset in this run. Independent
decode timings are close and do not establish a material policy difference.

The evidence supports using the tree as the L9 low-cost selector. A targeted
exact fallback for the sorted leaf-6 uncertainty region is the next policy
candidate to measure; the present run does not justify full exact selection
for all groups at L9.

Following this comparison, the level-9 preset uses the tree and level 10 uses
exact shared 108/109 scoring. Explicit selector overrides remain available for
research comparisons.
