# Full opcode-109 comparison, 2026-11-14

All six full L10 runs passed raw-ID verification, independent decoding and
opcode/stream accounting. Each dataset contains 1,107,205 vectors and
5,979,715,441 IDs. Native ARM NEON build, ten workers.

109-only excludes the flat-wave opcode 108. Other L10 competitors remain enabled.

| Dataset | Wave formats | BV bytes | Bits/ID | Encode s | Verify s | Decode s |
|---|---|---:|---:|---:|---:|---:|
| sorted | 108 only | 2,642,702,642 | 3.535556 | 1422.821 | 58.529 | 43.356 |
| sorted | 108 + 109 | 2,642,425,754 | 3.535186 | 1251.707 | 43.566 | 24.675 |
| sorted | 109 only | 2,646,099,354 | 3.540101 | 911.206 | 37.777 | 20.035 |
| unsorted | 108 only | 5,331,826,516 | 7.133218 | 1744.977 | 84.963 | 66.298 |
| unsorted | 108 + 109 | 5,331,821,868 | 7.133211 | 1792.964 | 82.719 | 64.925 |
| unsorted | 109 only | 5,331,830,720 | 7.133223 | 1408.448 | 82.477 | 64.166 |

Adding 109 saved 276,888 bytes on sorted and 4,648 bytes on unsorted
against the 108-only control. Retiring 108 increased size by 3,673,600 bytes
on sorted and 8,852 bytes on unsorted relative to allowing both formats.

For sorted, 109-only was 0.139024% larger than both, with observed encode time 27.20% lower.

| sorted opcode | 108-only records | Both records | 109-only records |
|---|---:|---:|---:|
| 68 | 717,865 | 712,387 | 948,255 |
| 92 | 7,049 | 6,203 | 6,799 |
| 98 | 6,105 | 6,105 | 6,105 |
| 108 | 1,032,542 | 1,004,047 | 0 |
| 109 | 0 | 34,819 | 802,402 |

For unsorted, 109-only was 0.000166% larger than both, with observed encode time 21.45% lower.

| unsorted opcode | 108-only records | Both records | 109-only records |
|---|---:|---:|---:|
| 68 | 2,163,762 | 2,163,014 | 2,164,312 |
| 92 | 1 | 0 | 0 |
| 98 | 4,016 | 4,016 | 4,016 |
| 108 | 5,880 | 5,488 | 0 |
| 109 | 0 | 1,141 | 5,331 |

Timings are single-pass observations. The control overlapped focused checks
during part of its compression. Both candidate configurations ran sequentially.
Control report assembly was recovered from completed validation logs and
reconciled measurement metadata after the running shell script was edited.
The frozen harness completed both candidate configurations without that issue.

Focused checks passed ASan/UBSan for width/payload/header variants, dictionary
hits and literal escapes, OR, range, indexed skip, final superblock and 32/64-bit
address builds. The portable codec check and released-opcode registry also passed.

[Full opcode bytes and counts, including all block formats](/Volumes/USB_1TB/BitMagic-experiments/gov2-op109-comparison-20261114.txt)
