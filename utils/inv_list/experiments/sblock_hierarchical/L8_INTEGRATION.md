# L8 tree preset and L9 policy measurements, 2026-11-14

L8 now selects one of opcode 108 and 109 using source block metadata. The
selected wave format still competes with the ordinary/adaptive superblock
record on complete padded bytes. Explicit serializer policy 1 retains exact
108/109 scoring for measurement and overrides. The L9 and L10 presets retain
their prior wave selector while their policies are studied separately.

The full normal-preset Gov2 run used native NEON and ten workers. Both full
collections passed raw-ID verification, independent decoding and opcode/stream
accounting. Preset streams were byte identical to the previous explicitly
selected tree streams, verified with SHA-256.

| Dataset | Vectors | IDs | Bytes | Bits/ID | Encode s | Verify s | Decode s |
|---|---:|---:|---:|---:|---:|---:|---:|
| sorted | 1,107,205 | 5,979,715,441 | 2,664,410,081 | 3.564598 | 430.345 | 43.179 | 22.455 |
| unsorted | 1,107,205 | 5,979,715,441 | 5,331,830,050 | 7.133222 | 907.739 | 82.078 | 65.550 |

The earlier paired exact L8 controls measured 2,664,247,433 sorted bytes
and 5,331,828,266 unsorted bytes. The L8 tree therefore costs 162,648 and
1,784 bytes respectively. Paired encode times were 537.787 versus 427.255
seconds sorted and 1046.737 versus 916.059 seconds unsorted. These timings
are single-pass whole-phase observations. The present normal-preset run has
matching stream sizes and additional timing observations above.

The exact L8 controls also evaluated the tree decision at every eligible
group. Regret below is extra complete bytes after ordinary/adaptive fallback.

| Leaf | Tree choice | Sorted groups | Sorted regret B | B/group | Unsorted groups | Unsorted regret B | B/group |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 108 | 54,491 | 4,900 | 0.090 | 15,825 | 56 | 0.004 |
| 1 | 109 | 2,787 | 1,224 | 0.439 | 0 | 0 | 0 |
| 2 | 108 | 7,666 | 3,016 | 0.393 | 1,521 | 116 | 0.076 |
| 3 | 109 | 93,332 | 28,060 | 0.301 | 4,777 | 260 | 0.054 |
| 4 | 109 | 310,835 | 112 | 0.0004 | 54,072 | 32 | 0.0006 |
| 5 | 109 | 188 | 24 | 0.128 | 3,402 | 8 | 0.002 |
| 6 | 108 | 1,288,157 | 125,312 | 0.097 | 2,090,046 | 1,312 | 0.0006 |

Trying both formats on leaves 1, 2 and 3 would recover 32,300 sorted bytes
across 103,785 groups, and 376 unsorted bytes across 6,298 groups, if the
L8 distribution and fallback are retained. Leaf 6 accounts for 125,312
sorted regret bytes across 1,288,157 groups, with a 2,768-byte worst miss.
Its large group count makes an unconditional exact trial expensive. Leaf 4
has negligible regret. A narrower gate inside leaf 6 requires per-record
measurements and a paired timing trial. These are L8 observations; L9 also
engages BWT and other selectors, so its policy needs L9-specific measurement.

Full run logs and opcode counts are under
/Volumes/USB_1TB/BitMagic-experiments/gov2-l8-integrated-tree-20261114.
The exact controls and per-leaf TSV files are under
/Volumes/USB_1TB/BitMagic-experiments/gov2-tree-l8-20261114-final.
