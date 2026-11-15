# Individual-vector zstd survey

Run from the BitMagic root with
`python3 utils/inv_list/experiments/vector_zstd/run_survey.py`.
Uses independent zstd level-3 frames on frozen v3 L8 BLOBs for both datasets,
replaying the saved distributed 20K ordinals. Every frame is decompressed and
compared byte-for-byte. Output stays on USB_1TB.

Retains up to 100 winners by bytes saved and 100 by percentage saved
(percentage ranking requires a 1 KiB source by default). The union is saved
as individual .bv and .bv.zst files with ordinal playlists and SHA-256 hashes.
Selective savings assume existing BLOB fallback and exclude a future wrapper
selector/index. Compression times cover native codec calls only.

## 2026-11-13 result

| Dataset | Sample bytes | Winning vectors | Selective saving | Forced frame bytes |
|---|---:|---:|---:|---:|
| sorted | 47,170,499 | 24 | 2,685 B (0.00569%) | 47,363,931 |
| unsorted | 96,676,984 | 0 | 0 B | 96,877,803 |

Whole-file zstd gains do not reproduce in independent sample-vector frames.
The largest sorted winner is ordinal 8872 (620 bytes, 4.36%); ordinal 14321
saves 189 bytes from a 955-byte BLOB (19.79%), dominated by fixed-wave opcode
90. Other top winners contain adaptive correction or fixed/adaptive waves.
Opcode composition alone does not locate the bytes zstd exploits.

All 24 saved sorted BLOBs passed the BM opcode profiler. Ignore AppleDouble
sidecars when enumerating the USB examples. Raw opcode ID counters are zero
because this profiler uses its default byte-accounting mode.

[Reports and retained examples](/Volumes/USB_1TB/BitMagic-experiments/gov2-vector-zstd-20261113-114602)

## Wave-map and position-run follow-up

Full frozen L8 replay covered 1,107,205 vectors per dataset. Gamma maps plus
repeated run-length pairs, with BIC fallback and three-bit hypothetical
selection, projected 951,132 bytes saved on sorted and 5,664 on unsorted.
Arithmetic maps, hierarchical maps, arithmetic-coded hierarchical literals,
and the standalone ordinary-run option were retired from the survey sources.

Experimental opcode 108 now combines six adaptive widths with BIC/gamma maps
and BIC/repeated-pair positions. The initial trial compared its complete dry-scored
record against 90/97 and other superblock candidates. Opcodes 90/97 were
subsequently retired after the full consolidation comparison. Map polarity,
uniform correction and BIC endpoint alternatives remain available. All fields
share one accumulator and final padding. The subsequent full L10 run saved
1,164,008 bytes on sorted and 16,448 on unsorted against v3. See the
shared-selector benchmark notes in ../sblock_wave_adaptive/RESULTS.md.

[Full measurement reports](/Volumes/USB_1TB/BitMagic-experiments/gov2-full-wave-runs-20261113)

## Full L10 local pair-dictionary measurement

The 2026-11-14 measurement covered both full v4 L10 collections. Sorted had
18,852 paired-run opcode-108 records in 18,035 vectors; unsorted had none.
Dictionary caps of 4, 8 and 16 entries saved 69,828, 78,624 and 75,128 bytes
respectively with current-108 fallback. Best-per-record selection saved
87,844 bytes (0.003324% of sorted), winning on 1,499 records. All entry counts,
pair values, tags, gamma ranks, repeat counts, two extension-subformat bits
and final word padding were charged. No candidate streams were emitted.

Full replay took 126.388 seconds sorted and 111.461 seconds unsorted; paired
record parsing/scoring accounted for 0.361 seconds. Vector 17073's large
record projected a 6,872-byte saving (14,213 to 7,341 bytes). This remains a
promising localized pattern, with modest measured whole-collection benefit.
Decision: defer the opcode-108 dictionary extension. For future datasets use
a cap of 8 entries as the single-cap benchmark. On sorted Gov2 it captured
78,624 of the 87,844 bytes available from choosing among all three caps
(89.5%), or 0.002975% of the collection. The top 100 cap-8 records contributed
50,832 bytes. Unsorted Gov2 had no paired-run records to test. The dry scorer
and its runner are retired; the reports retain the per-record comparisons.

[Dictionary measurement report](/Volumes/USB_1TB/BitMagic-experiments/gov2-pair-dictionary-v4-l10-20261114/report.txt)

## Full L10 block/wave occupancy-map measurement

The 2026-11-14 measurement covered both full v4 L10 collections. Sorted
superblock records contained 217 full bit-blocks, all in 214 positive-polarity
opcode-98 records. Unsorted superblock records contained no full bit-blocks.
Neither collection emitted inverse-polarity opcode 98.

For opcode 108 at its emitted width, the first two-level map scorer saved
17,136 bytes on 4,139 sorted records and 44 bytes on 11 unsorted records with
current-format fallback. Forced replacement increased complete size by
21,412,032 and 68,100 bytes respectively. Replay took 540.634 seconds sorted
and 132.821 seconds unsorted. These are measurement timings, not codec timings.

The root block map and each partial local wave map independently select BIC
or gamma runs and occupied/empty polarity. BIC uses the full coordinate bounds.
Each map pays a gamma-coded list count even for run coding. There is a root
all-occupied flag and a per-occupied-block all-waves-occupied flag; local partial
maps pay two selector bits. Two extension-subformat bits and complete-record
padding are included. All fields are assumed to share one accumulator.

The result remains under review before retiring H-map. Further framing analysis
should consider deriving counts for gamma runs, sharing or deriving selector
bits, explicit raw-mask alternatives, corrected BIC coordinates and bounds,
and grouping local status fields. The survey held wave width and the position
payload fixed. No candidate streams were emitted.

[Occupancy survey report](/Volumes/USB_1TB/BitMagic-experiments/gov2-occupancy-v4-l10-20261114/report.txt)

### FULL-block extraction on Gov2.sorted

The v4 L10 dry replay covered the same 214 positive opcode-98 records
containing 217 FULL blocks (4,816,518 emitted bytes). With original-98 fallback,
extracting FULL blocks and BIC-coding the residual saved 19,032 bytes on 135
records. The global minimum-gap correction was zero on every residual.

Extracting FULL blocks and trying six wave widths saved 33,116 bytes on 188
records. Its complete record was smaller than the residual-BIC candidate on
all 214. The residual-BIC trial is retired. The six-width result amounts to
0.69% of the selected records' bytes and 0.00125% of the sorted collection.
The unsorted superblock census found no FULL blocks. A FULL-only selector is
a research option pending evidence from workloads with more FULL blocks.

[Paired record scores](/Volumes/USB_1TB/BitMagic-experiments/gov2-full-residual-v4-l10-20261114/records.tsv)
and [dry-trial report](/Volumes/USB_1TB/BitMagic-experiments/gov2-full-residual-v4-l10-20261114/report.txt).

### Two-level wave-map baseline

The map-only trial used six wave widths. Store the
root block map, concatenate the surviving partial blocks, and encode one wave
occupancy map over that shortened domain using BIC or gamma runs. With fallback
to the emitted record, complete-record dry savings were 116,172 bytes against
sorted opcode 108, 1,632 bytes against unsorted opcode 108, and 33,116 bytes
against opcode 98 on the separate 214-record sorted FULL subset.

Retire the nested BM bit-block selector trial. Adding that choice increased
complete size relative to the two-level BIC/gamma layout by 1,998, 153 and
17 bytes respectively. On sorted opcode 108 records, embedded block encodings
recovered 7,570 bytes while signaling the choice cost 9,568 bytes after
padding and fallback.

[Map-only BIC/gamma trial](/Volumes/USB_1TB/BitMagic-experiments/gov2-concat-map-v4-l10-20261114/report.txt)
and [nested-selector comparison](/Volumes/USB_1TB/BitMagic-experiments/gov2-inner-map-v4-l10-20261114/report.txt).

### Candidate for future datasets

Retain one design for further measurement: the two-level BIC/gamma wave map
with an optional local dictionary of up to eight run/gap pairs in the position
payload. Score the dictionary at each of the six wave widths. Defer the wire
opcode decision until the post-cleanup discussion.

On the full v4 level-10 replay, complete-record dry savings against emitted
opcode 108 were 207,148 bytes on sorted Gov2 and 1,632 bytes on unsorted.
The dictionary added 90,976 bytes of savings over the map alone on sorted
and zero on unsorted. Keeping a separate dictionary path on the original
wave map would add just 708 bytes on sorted and zero on unsorted.
Replay took 111.22 s and 19.36 s respectively with eight workers.

[Combined trial report](/Volumes/USB_1TB/BitMagic-experiments/gov2-map-dictionary8-v4-l10-20261114/report.txt)

The subsequent emitted-stream experiment uses opcode 109. All six full L10
runs validated. See [the three-way comparison](../sblock_hierarchical/RESULTS.md)
for the measured gain from adding 109 and the separate cost of retiring 108.

[Shared exact scoring and L8 predictor survey](../sblock_hierarchical/SHARED_SCORING.md) records the byte-identical refactor and strict-win/tie diagnostics.
