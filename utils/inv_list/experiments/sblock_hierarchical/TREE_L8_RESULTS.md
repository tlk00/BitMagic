# Experimental L8 tree selector, 2026-11-14

The source-metadata tree chooses one of opcodes 108 and 109. It uses the
number of level-zero GAP blocks, total GAP length, residual set-bit count,
occupied-block run count, and occupied-block span. The chosen format still
competes against the ordinary and adaptive superblock records by complete
padded bytes. Opcode 98 and other L8 options remain enabled.

Both full Gov2 sets passed raw-ID verification, independent decode, and
opcode accounting. Each stream holds 1,107,205 vectors and 5,979,715,441 IDs.
All four streams used the same frozen native NEON build with ten workers.

| Dataset | Policy | Bytes | Bits/ID | Encode s | Verify s | Decode s |
|---|---|---:|---:|---:|---:|---:|
| sorted | Exact | 2,664,247,433 | 3.564380 | 537.787 | 43.747 | 22.754 |
| sorted | Tree | 2,664,410,081 | 3.564598 | 427.255 | 43.426 | 22.698 |
| unsorted | Exact | 5,331,828,266 | 7.133220 | 1046.737 | 81.843 | 66.342 |
| unsorted | Tree | 5,331,830,050 | 7.133222 | 916.059 | 88.403 | 63.808 |

Sorted regret was 162,648 bytes (0.0061%) and encoding took 110.532 fewer
seconds (20.6%). Unsorted regret was 1,784 bytes (0.000033%) and encoding
took 130.678 fewer seconds (12.5%). Decode observations differed by -0.056 s
and -2.534 s. These are single-pass whole-phase timings.

| Dataset | Opcode | Exact records | Tree records |
|---|---:|---:|---:|
| sorted | 68 | 712,387 | 715,313 |
| sorted | 92 | 6,203 | 6,508 |
| sorted | 98 | 6,405 | 6,405 |
| sorted | 108 | 1,004,047 | 961,772 |
| sorted | 109 | 34,819 | 73,863 |
| unsorted | 68 | 2,163,014 | 2,163,310 |
| unsorted | 92 | 0 | 0 |
| unsorted | 98 | 4,016 | 4,016 |
| unsorted | 108 | 5,488 | 4,528 |
| unsorted | 109 | 1,141 | 1,805 |

Exact scoring also applied the frozen tree to every candidate and accumulated
complete-byte regret after the incumbent fallback. The leaf regrets sum exactly
to both emitted-stream size differences.

| Leaf | Choice | Sorted groups | Sorted losses | Sorted regret B | Worst B | Unsorted regret B |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 108 | 54,491 | 1,214 | 4,900 | 16 | 56 |
| 1 | 109 | 2,787 | 306 | 1,224 | 4 | 0 |
| 2 | 108 | 7,666 | 679 | 3,016 | 28 | 116 |
| 3 | 109 | 93,332 | 7,015 | 28,060 | 4 | 260 |
| 4 | 109 | 310,835 | 28 | 112 | 4 | 32 |
| 5 | 109 | 188 | 6 | 24 | 4 | 8 |
| 6 | 108 | 1,288,157 | 12,923 | 125,312 | 2,768 | 1,312 |

Leaf 6 accounts for 77% of sorted regret, including the worst miss. A future
L9 exact fallback should be evaluated on that leaf first. Its 1.29 million
sorted groups also make CPU cost important. Leaf 3 is the next largest source
of sorted regret at 28,060 bytes.

The collection-wide sorted regret was much larger than the earlier sampled
cross-validation prediction of 2,848 bytes. The sample missed uncommon costly
cases in leaf 6. The tree was integrated into the L8 preset after this paired trial.
-wave-tree remains an explicit override for other levels.

Complete opcode byte counts and all leaf rows are in
/Volumes/USB_1TB/BitMagic-experiments/gov2-tree-l8-20261114-final/comparison.txt.
Reproduce with run_tree.sh FRESH_USB_DIR /path/to/bminv 10 after make neon.
