# Horizontal XOR survey of Gov2 posting lists

From `utils/inv_list`, build with `make gov2-posting-xor` and run
`python3 experiments/horizontal/run_posting_xor.py`. The default is 10,000
vectors per dataset on USB_1TB. Use `--help` for paths, sample size and seed.
The runner indexes raw record lengths and chooses one vector per equal-count
stratum with the same seed on each collection. It saves the zero-based ordinals,
per-dataset metrics and a paired report, with no serialized output files.

Each occupied group of 256 blocks is rebased to zero and optimized. Search
compares every native BIT target with every earlier native BIT source in that
group. GAPs are inventoried by actual run count and allocation level. Their
standalone record costs are hypothetical when the group codec wins.

The population scorer chooses the mask giving the smallest whole-block minority
population. The run-aware scorer uses BM's wave complexity descriptor. Each
keeps four positive-scoring sources (lowest slot first on ties). Trial the
proposed mask and full mask for each, sharing identical trials. Separate exact
wave counts exclude constants and require an earlier block. Different-offset-only
counts target waves with no exact match at the same offset.

**Reading report.tsv.** Baseline bytes sum actual in-memory level-10 streams
for isolated groups, including their individual headers. Block savings compare
native records with actual residual records plus projected metadata: two bytes
for opcode and source, plus eight for a partial mask. An exact copy omits the
residual. Sources can themselves be reconstructed earlier targets. Group
savings subtract winning block savings from the measured ungrouped stream size
and compare with the cheaper of ordinary level 10 and the ungrouped control.
The union credits each target once. This measures the incremental opportunity
under the assumed framing, not a whole-collection wire size. The detailed TSV
also gives search-only timings (excluding shared preparation, ranking and encoding),
formula bytes, winning distances and mask sizes. Missing histogram bins are zero.

**Gather and range decoding remain open.** A requested target can depend on an
earlier source outside the requested range, transitively. A future format must
define dependency discovery or bounded group replay, index/bookmark behavior,
and its access cost before adoption. This survey adds no opcode or decoder.

## Whole-vector follow-up on Gov2.sorted, November 11, 2026

Keep the existing level-10 selection for opcode 92 (`bm::set_sblock_bienc_adaptive`).
The experimental selector tried individual blocks after a mixed BIT/GAP group
selected opcode 92. It won 73 of 164 trials across 50,000 posting lists
(269,477,079 IDs). Both streams decoded to the original vectors.

| Whole-vector measurement | Level 10 | With constituent trial |
|---|---:|---:|
| Bytes | 117,173,972 | 117,139,265 |
| Single-thread encode seconds | 325.570 | 327.737 |

The trial saved 34,707 bytes (0.0296%) at an additional 2.167 seconds
(+0.67%). The isolated-group projection was 34,774 bytes. The measured
whole-vector gain did not justify adding the selector to level 10.

[Isolated-group profile](/Volumes/USB_1TB/BitMagic-experiments/gov2-sorted-group-profile-20261111-102932/report.tsv) · [Whole-vector report](/Volumes/USB_1TB/BitMagic-experiments/gov2-sorted-whole-group92-20261111-104718/report.txt)

## Displaced-wave copy trial, November 11, 2026

On the same 10,000 posting lists per dataset, the survey measured copying
nonconstant 1,024-bit waves from an earlier BIT block at a different offset.
The formula charged two bytes per target block and 20 bits per copied wave,
rounded up to whole bytes per formula.

| Dataset | Exact matches | Copies in winning blocks | Net block saving | Net saving after superblock selection |
|---|---:|---:|---:|---:|
| Gov2.sorted | 859 | 73 | 245 B | 174 B |
| Gov2.unsorted | 17 | 1 | 3 B | 3 B |

These projected gains were too small to pursue a copy opcode. The trial was
retired. The horizontal survey still counts exact repeats at different offsets.
[Paired profiling report](/Volumes/USB_1TB/BitMagic-experiments/gov2-displaced-waves-20261111-111257/report.tsv)

## Arithmetic payloads after superblock wave compaction (retired, 2026-11-13)

Two distributed 20K-vector samples compared actual arithmetic chunk sizes with
the level-10 superblock selector across six wave widths (128 through 4,096 bits).
Trials used raw or minimum-gap-corrected compacted bitmaps. Savings below are
projections from complete record sizes with the existing production winner retained.

| Sample | Groups surveyed | Byte arithmetic saving | Extra nibble saving | Groups beating production |
|---|---:|---:|---:|---:|
| Gov2.sorted | 31,652 | 620 B | Untested | 2 |
| Gov2.unsorted | 39,146 | 0 B | 416 B | 3 |

The sorted run also scored a replacement selector that disabled opcode 90 and
tagged the BIC payload. It increased the 46,817,900-byte baseline by 6,328 bytes.
The unsorted trial retained opcode 90 and projected 416 bytes saved on a
96,676,902-byte baseline. The gains did not justify adding payload trials to
production. Retain fixed-wave opcode 90 and adaptive-wave opcode 97 with BIC.

[Sorted report](/Volumes/USB_1TB/BitMagic-experiments/gov2-wave-payload-20261113-082148/report.tsv) ·
[Unsorted report](/Volumes/USB_1TB/BitMagic-experiments/gov2-wave-payload-20261113-083212/report.tsv)
