# Full Gov2 L8 per-vector zstd replay (tree preset)

Input: the validated, independently decoded L8 tree streams at
`/Volumes/USB_1TB/BitMagic-experiments/gov2-l8-integrated-tree-20261114`.
Every one of 1,107,205 serialized vectors in each collection was passed through
an independent zstd 1.5.7 level-3 frame and decompressed byte-for-byte. The
original BM stream is kept whenever its independent frame is no smaller.
Reports, all-vector TSVs, and the top 100 absolute and percentage examples are
under `/Volumes/USB_1TB/BitMagic-experiments/gov2-zstd-full-tree-l8-20261114`.
The same directory retains `retained_wave_diagnostics.tsv` for selected
examples and `all_sorted_winner_opcodes.tsv` /
`all_unsorted_winner_opcodes.tsv` for every winning vector.

| Dataset | BM bytes | Winning vectors | Selective saving | Saving / BM bytes | All zstd frames minus BM |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gov2.sorted | 2,664,410,081 | 738 | 113,408 B | 0.004256% | +10,753,757 B |
| Gov2.unsorted | 5,331,830,050 | 4 | 168 B | 0.00000315% | +11,116,469 B |

The earlier v4 L8 108-only stream had 784 sorted winners and 123,688 B of
selective saving. The tree stream has 46 fewer winners and 10,280 B less
remaining zstd saving. The same four unsorted ordinals still win, by 36–50 B
each on 132–397 KiB vectors. The total sorted selective gain remains only
about 0.0043% of the BM collection.

## Largest current winners

| Sorted ordinal | BM bytes | zstd saving | Dominant records |
| ---: | ---: | ---: | --- |
| 14831 | 14,136 | 3,398 B | 108: 12,325 B; 68: 1,809 B |
| 13047 | 9,132 | 3,273 B | two 108 records: 8,589 and 541 B |
| 17076 | 15,552 | 2,212 B | 108: 12,213 B; 68: 3,337 B |
| 14205 | 15,160 | 2,118 B | 108: 12,949 B; 68: 2,209 B |
| 3672 | 89,240 | 1,740 B | two 98 records: 46,449 and 42,789 B |
| 5164 | 104,740 | 1,660 B | two 98 records: 55,949 and 48,789 B |
| 22487 | 5,551 | 1,445 B | 109: 5,545 B |
| 7950 | 41,372 | 1,427 B | 98: 22,625 B; 108: 18,745 B |

Independent zstd of the opcode 108 record alone saves 3,316 B on ordinal
14831 and 3,310 B on ordinal 13047; their opcode 68 or second short 108
records do not compress independently. The two leading records both use
128-bit waves, an empty-wave BIC map, and the repeated run/gap-pair position
payload (108 layout mode 10). Their map costs are only 5,354 and 2,126 bits,
while their position payloads cost 93,155 and 66,500 bits. Thus the large
zstd saving cannot come from map bytes alone.

In the corrected position lists, the eight most frequent `(ones, gap)` pairs
account for 16,223 of 17,624 pair occurrences on ordinal 14831 and 18,003 of
19,034 on ordinal 13047. The most frequent pairs are `(1, 7)` occurring 3,596
times and `(2, 6)` occurring 10,760 times, respectively. The pair streams still
contain 11,556 and 7,731 distinct adjacent groups: the existing writer
compacts consecutive repetitions, while many of these pairs recur later.
Repeated 8-byte windows cover 11.2% and 10.5% of the raw opcode 108 records.
Five similarly sized, nonwinning opcode 108 controls had no repeated 8-byte
windows and gained nothing when compressed as separate records. This supports
repeated nonadjacent payload sequences as the main zstd opportunity; it is
not an attribution of zstd's internal match decisions.

Exact L8 108/109 scoring would select the 109 cap-8 pair-dictionary format on
these two groups, shrinking their wave records by 2,076 and 1,872 B. The
whole vectors then still offer 1,608 and 1,250 B of zstd saving. Ordinal
17073 illustrates the stronger case: the old 14,213-byte opcode 108 record
became a 7,433-byte mode-3 opcode 109 record under the tree preset. Its whole
vector's zstd opportunity fell from 4,898 to 300 B. The remaining winners
therefore include both tree misses and repetition beyond the cap-8 dictionary.

Other winners are not all wave records. Independently compressed records in
the winning sorted vectors have positive savings in opcode 92 (adaptive
superblock), opcode 98 (large superblock), and opcode 63 (block position list)
as well. Their per-record zstd gains are diagnostic only: separate frame
overhead and matches between records prevent summing them into the collection
gain. The full per-vector numbers above are the valid end-to-end comparison.

This replay suggests no general zstd wrapper for L8 vectors. Further codec
work, if warranted on other datasets, should measure higher-order pair/tag
sequences in the wave position payload and compare complete record sizes with
the existing 109 dictionary and ordinary fallback.
