> Historical experiment notes. The exploratory measurement code and old
> runners described below were retired on 2026-11-01. Use
> [current highest-compression verification](../current/README.md).
> Retained codec correctness tests remain where applicable.

# Superblock wave-size and occupancy-map diagnostics

Optional `-sblock-map-trial` observes final records for **68, 90 and 92**.
It decodes their original positions and measures hypothetical wave-compacted
records without changing the serializer's decisions or emitted bytes.
This includes unsorted records where opcode 90 currently loses.

Widths: **128, 256, 512, 1024, 2048, 4096 bits**. Every width recomputes
occupied/empty maps, stable compaction, global minimum-delta correction,
residual bounds and BIC costs. All indexes use unsigned 32-bit integers:
128-bit waves need 131,072 possible wave indexes.

Map codecs: ordinary BIC, min-delta BIC, hierarchical groups of 8/32/128
(two-bit zero/one/literal tag), and gamma run lengths over the occupancy bitmap.
Both map polarities and both residual boundary layouts are tried. Map count
remains in the header even for codecs that could derive it.

Three independent arms:

| Arm | Choices | Additional selector bits |
|---|---|---:|
| Wave size | Six widths, ordinary BIC map | 3 |
| Map codec | Six codecs, fixed 1024-bit width | 3 |
| Joint | Six widths times six codecs | 6 |

Each compares **complete records** against the actual winning 68/90/92 record,
including one opcode byte, metadata and final 32-bit accumulator padding.
Losses and ties retain the existing record. Separate fixed-width BIC rows omit
selector bits to estimate the best constant-width alternative; rows/arms are
not additive. The 1024-bit BIC control must reproduce emitted opcode-90 sizes.
No-empty-wave cases simply retain baseline; no candidate is invented.

Reports include source opcode, winning width/codec counts, savings, collection
weighted bits per ID, header/map/residual/padding costs and summed probe-worker
seconds. Retained baseline records are opaque: their bytes appear in the header
accounting bucket, not as a fabricated map or residual breakdown. Detailed raw
rows remain in compress.log; the final tables are printed once per dataset.

## Build and full paired run

From `utils/inv_list`:

```bash
make neon
bash experiments/sblock_map/run_gov2_sblock_map.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --datasets both --threads 10 --progress 1000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-11-01/runs/inv_list/gov2-sblock-maps-expanded-$(date +%Y%m%d-%H%M%S)"
```

No max-vectors limit means full collections. Add `--max-vectors 10000` for
10K per set. Optional `--reference RUN_ROOT` compares scope and SHA-256 with
an existing matching run's sorted/profile and unsorted/profile. Baseline options
are all-current selection: adaptive superblocks plus opcode 90, GAP-XOR,
shared original-polarity BWT 72/96, 16-byte cutoff, opcode-65 endpoint selection.

Harness checks the USB mount and free space, compresses, verifies against raw
IDs and independently decodes. Verified streams are deleted by default; use
`--keep-streams` to retain. Inputs/sources remain on the system drive. Each set
gets `sblock-map-report.txt` and `profile/sblock-map-measurement.json`.
This is a diagnostic run, not a candidate wire-format or decoder-speed test;
encoding time includes substantial trial work (up to 144 record layouts).

Focused test:

```bash
clang++ -std=c++17 -O1 -g -fsanitize=address,undefined -DBMNEONOPT -Isrc \
  utils/inv_list/experiments/sblock_map/test_sblock_map.cpp -o /tmp/test_sblock_map
/tmp/test_sblock_map
```

Tests materialize candidate maps/records and independently restore positions
across widths, polarities and boundary modes, comparing full record costs.
