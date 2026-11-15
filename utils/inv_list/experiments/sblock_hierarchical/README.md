# Experimental hierarchical-wave superblocks

Opcode 109 (`set_sblock_bienc_wave_hierarchical`) stores an occupied block
map followed by an occupied wave map in the concatenated occupied blocks.
Each map selects BIC with polarity or gamma-coded runs. Six wave widths
compete. Corrected positions select BIC, explicit-endpoint BIC, repeated
run/gap pairs, or pairs with a local dictionary capped at eight entries.

The format has one reserved-zero extension bit after the width. A reserved
FULL-list flag is zero in this sparse experiment; fewer than 65,536 source
positions cannot contain a FULL block. Opcode 98 remains a separate selector.

Build with `make neon` in `utils/inv_list`. Run three fresh output directories
through `experiments/current/run_gov2_highest.sh` with `--datasets both
--level 10 --threads 10`, selecting `--wave-mode 108`, `--wave-mode both`
or `--wave-mode 109`. These modes affect wave formats only; other preset
competitors remain enabled. Level 8 and 9 presets use the source-metadata
tree; level 10 exactly scores both 108 and 109. Explicit `--wave-mode`
overrides the preset for comparisons.

For a fresh sequential matrix, use `run_compare.sh --bin ./bminv --out FRESH_USB_DIR`
from `utils/inv_list`. It runs all three configurations and saves `comparison.txt`.

The harness checks raw-ID equality, independent decoding and opcode accounting.
Run `report_compare.py CONTROL BOTH HIERARCHICAL` after all three complete.
It reports absolute sizes, savings against 108-only, retirement cost relative
to both, encode/verify/decode times and every emitted opcode.

Focused regression sources are `tests/regression/serial_sblock_hierarchical.cpp`
and `serial_sblock_hierarchical_routes.cpp`. They exercise all six widths and
four payloads, profile traversal, OR, range and indexed decoding.

[Full sorted/unsorted results](RESULTS.md). Opcode 108 remains available pending the retirement decision.

[Shared exact scoring and L8 predictor survey](SHARED_SCORING.md).


## Occupied-block gate experiment

Mode `gate16` skips opcode 108 scoring when at most 16 source blocks contain
set bits. It still compares the selected 109 record against ordinary/adaptive
fallback. Larger groups use exact shared selection. This is an explicit
experiment; compression-level presets are decided separately.

Build with `make neon` in `utils/inv_list`, then run
`bash experiments/sblock_hierarchical/run_gate16.sh --bin ./bminv --threads 10 --out FRESH_USB_DIR`.
The runner validates exact and gated L8 streams on both full Gov2 collections
and writes `comparison.txt` with size changes, timings and all opcode counts.

## Selection rules and metrics

```text
N = set-bit count
B = occupied source block count
incumbent = best ordinary/adaptive record

formats = {108, 109}
if gate16 and B <= 16:
    formats = {109}

for width in {128, 256, 512, 1024, 2048, 4096}:
    compact occupied waves and apply minimum-gap correction
    score each ordinary position payload once
    for each enabled format:
        score its maps, metadata and allowed payloads
        record_bytes = 1 + 4 * ceil(total_field_bits / 32)
        retain this format's first strict complete-byte improvement

winner = incumbent
if best108.bytes < winner.bytes: winner = best108
if best109.bytes < winner.bytes: winner = best109
write winner once
```

At the same width, both formats have identical wave-compacted corrected
positions. Opcode 108 describes occupied/empty waves in the original superblock.
Opcode 109 stores a root block map and a wave map in the concatenated occupied
blocks. Ordinary position payloads are BIC, BIC with explicit endpoints and
repeated run/gap pairs. Opcode 109 additionally permits a cap-8 pair dictionary.

The score includes header fields, map counts/selectors, the position payload
and final word padding. Gamma lengths and recursive BIC lengths are exact
dry counts. A small bit difference can become a four-byte difference after
padding. Incumbents win equal byte counts; 108 wins byte ties against 109.

The gate uses B alone. N, occupied-wave count, map cost, correction and payload
cost govern the subsequent exact width/layout choice. Gate regret is the
difference between gated and exact complete sizes with incumbent fallback.
Search statistics and emitted-opcode counts have different scopes.

[Full L8 gate16 results and level-policy discussion](GATE16_RESULTS.md).

### FULL-block extension trial

The opcode-109 extension flag stores a FULL-block map in original
block coordinates. The existing root map and wave map cover only the remaining
partial blocks. Levels 8–10 compare complete records against opcode 98 and the
ordinary stream, retaining strict fallback. The first implementation accepts 2–131,071 residual
IDs; zero and one residual ID remain for a later format decision. FULL groups
remain outside the 108/109 tree training.

The full Gov2 Level 8 trial passed raw-ID verification and independent decode
on both datasets. Sorted fell from 2,664,410,081 to 2,664,376,468 bytes
(33,613 bytes saved); opcode-109 count rose by 193. Unsorted stayed at
5,331,830,050 bytes and its stream was byte-identical to the control.
Encode times were 427.898 vs 430.345 seconds sorted and 897.297 vs 907.739
seconds unsorted, with normal run-to-run noise. See the retained reports in
`/Volumes/USB_1TB/BitMagic-experiments/gov2-l8-hier-full-20261114`.

The standalone predictor_export.cpp samples the same every-55th vector in
both raw Gov2 files, using the level 8 exact selector. Its arguments are the
raw input and FULL-group output file. Redirect stdout to sorted.features or
unsorted.features. train_predictor.py takes their common directory. Costs in
the last two columns include the incumbent fallback. Regret is measured at
wave selection, before any later large-superblock replacement.

[Full L8 tree selector results](TREE_L8_RESULTS.md).

[L8 preset integration and L9 leaf measurements](L8_INTEGRATION.md).

[Full L9 selector comparison](L9_RESULTS.md).
