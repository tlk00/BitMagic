# Full L8 gate16 results, 2026-11-14

All four full streams passed raw-ID verification, independent decoding and opcode/stream accounting.
Each dataset contains 1,107,205 vectors and 5,979,715,441 IDs. Native ARM NEON, ten workers.

The gate selects only 109 when at most 16 source blocks contain set bits.
Larger groups use exact shared competition between 108 and 109. Every group
retains its ordinary/adaptive fallback. Each enabled format searches six widths.
The rule uses logical occupancy independently of GAP/BIT representation.

| Dataset | Selection | BV bytes | Bits/ID | Encode s | Verify s | Decode s |
|---|---|---:|---:|---:|---:|---:|
| sorted | Exact | 2,664,247,433 | 3.564380 | 537.990 | 44.127 | 22.332 |
| sorted | Gate16 | 2,664,532,485 | 3.564762 | 427.825 | 43.290 | 23.209 |
| unsorted | Exact | 5,331,828,266 | 7.133220 | 1036.048 | 81.218 | 63.989 |
| unsorted | Gate16 | 5,331,830,458 | 7.133223 | 1035.958 | 81.474 | 63.784 |

sorted: size loss 285,052 bytes (0.010699%). Observed encoding change -20.48%, decoding change +3.93%.
unsorted: size loss 2,192 bytes (0.000041%). Observed encoding change -0.01%, decoding change -0.32%.

The sorted loss exactly matches the 285,052-byte survey prediction.
Encoding improves on sorted and is effectively unchanged on unsorted.
Sorted decode time increased in this pass. Unsorted decode changed by only
0.205 s. These are single-pass whole-phase timings. Compression includes
input conversion, serialization, profiling and output. Decode includes file
reading. These observations do not establish intrinsic per-record decode costs.

## Superblock opcode counts

| Dataset | Opcode | Exact | Gate16 |
|---|---:|---:|---:|
| sorted | 68 | 712,387 | 739,762 |
| sorted | 92 | 6,203 | 6,296 |
| sorted | 98 | 6,405 | 6,405 |
| sorted | 108 | 1,004,047 | 788,985 |
| sorted | 109 | 34,819 | 222,413 |
| unsorted | 68 | 2,163,014 | 2,163,407 |
| unsorted | 92 | 0 | 0 |
| unsorted | 98 | 4,016 | 4,016 |
| unsorted | 108 | 5,488 | 3,304 |
| unsorted | 109 | 1,141 | 2,932 |

Some former 108 records fall back to 68/92 when 109 loses to that incumbent.
Many equal-size 108/109 choices also switch to 109 under the gate.

## Level-policy discussion

L7/L8 decoder-biased coefficients still need matched-record decode timing.
The gate provides an encoder benefit on sorted; its decode result does not
yet justify using it as a decoder-speed preference.

L7 and L8 can share a wave policy while other techniques distinguish them.
L9 can use cheap preferences in clear regions and exact dry scoring in an
uncertainty zone. L10 can score every permitted candidate exactly and write
the winner once. Full trial output gives no extra size-selection accuracy
over the exact scorer for the current BIC/gamma/dictionary formats.

This gate is a partial-selection policy: direct 109 for small occupied-block
counts, exact competition elsewhere. Learned coefficients and level defaults
remain a separate design decision.

[Complete opcode byte/count report](/Volumes/USB_1TB/BitMagic-experiments/gov2-l8-gate16-20261114/comparison.txt)

The focused scorer regression passed 552 decisions under ASan/UBSan.
Build separately with make neon in utils/inv_list, then use
run_gate16.sh --bin ./bminv --threads 10 --out FRESH_USB_DIR.
