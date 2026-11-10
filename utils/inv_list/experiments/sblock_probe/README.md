> Historical experiment notes. The exploratory measurement code and old
> runners described below were retired on 2026-11-01. Use
> [current highest-compression verification](../current/README.md).
> Retained codec correctness tests remain where applicable.

# Ordinary superblock BIC and corrected-run survey

`-sblock-probe` observes only emitted opcode-68 records. The serializer counts
bits written for the opcode and other header fields, the global minimum-delta
correction metadata, BIC of the interior position list, and the final bit-writer
padding. Their sum is checked against the record bytes. The correction flag
is one bit within the existing header flag byte; the separate correction
field counts only the code and value. A dry BIC pass on the same uncorrected
position list measures the actual BIC gain or loss from the correction,
without changing the emitted record. It also counts
consecutive positions in the original sorted array and in the actual array
after the opcode-68 global correction. Run lengths use log2 buckets; coverage
for runs of at least 2, 4 and 8 positions is exact. It counts empty 65,536-bit
blocks, empty 1,024-bit waves, and interior empty-wave stretches before and
after correction. Leading and trailing empty waves are separate, so moving
positions into a dense prefix is visible. Rejected opcode-68
candidates replaced by adaptive opcode 92 do not enter these totals.

The harness also enables the emitted-opcode profiler. Its report checks that
the probe's count and bytes equal the opcode-68 row in the resulting stream,
then verifies the stream against raw Gov2 IDs and independently decodes it.
All reporting stays on USB_1TB; inputs, scripts, and code stay in this repo.
The verified `.bv` is removed unless `--keep-streams` is passed.

Build from `utils/inv_list` with `make neon`, then run:

```bash
bash experiments/sblock_probe/run_gov2_sblock_probe.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --datasets sorted --threads 10 --progress 10000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-sorted-sblock-probe-$(date +%Y%m%d-%H%M%S)"
```

The runner uses the established full-collection level-6 selection:
`-block-select-step2 -periodic-step1-interleaved 4 -bwt-min-gain 16 -gap-xor
-sblock-gap-tolerance 100 -sblock-adaptive`. `--max-vectors N` limits a
preliminary run. The final analysis is in
`sorted/sblock-probe-report.txt`; the opcode table is in `report.txt`.

## Dry run-endpoint candidate

Add `--sblock-run-trial` to the same command for a complete-record candidate
trial and two distributions. This option groups emitted opcode-68 records by
the percentage of corrected positions in runs of length at least four and by
the emitted BIC bits per position. A trial is eligible when
`2 * run_count < position_count` (fewer endpoints than positions).

Each trial uses a one-byte placeholder for a future opcode. The candidate
stores the superblock index, run count, global minimum correction, first run
start, reverse distance from the superblock end to the final exclusive end,
and BIC of the interior interleaved half-open run endpoints. The first and
last endpoints use 24-bit coordinates. The trial measures the actual header,
BIC, and final padding bytes. A decoder could expand runs and reverse the
global correction by position rank; this experiment does not add that decoder
or emit the candidate. `sorted/sblock-run-trial-report.txt` shows forced-trial
cost and projected per-record wins with ordinary fallback. These are
projections, not a verified compressed collection. The existing stream is
still verified and decoded, and the run trial leaves its bytes unchanged.

Full Gov2.sorted trial:

```bash
bash experiments/sblock_probe/run_gov2_sblock_probe.sh \
  --sblock-run-trial --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --datasets sorted --threads 10 --progress 10000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-sorted-sblock-run-trial-$(date +%Y%m%d-%H%M%S)"
```

Future experiment: compare constituent blocks pairwise after the global
minimum-delta correction, looking at intersections and residual positions.
The corrected array should be divided by original 65,536-bit block membership,
and local coordinates must be normalized before comparison. This is only a
diagnostic idea; no pairwise search or codec is added here.

## Dry empty-wave compaction trial

`--sblock-wave-trial` measures emitted opcode-68 records with at least one
empty 1,024-bit wave. Wave occupancy is computed from the original position
list, before ordinary minimum correction. Empty waves are removed and the
packed positions then receive global minimum-delta correction (minimum
adjacent difference minus one, including zero; no effort gate).

The scratch wire model uses a one-byte opcode field, two mode bits (map
polarity and boundary policy), a two-bit superblock-index mode and its index,
gamma position count, gamma selected-wave-list count, and gamma(correction+1).
The selected wave list is BIC in [0,16383]. Both empty and occupied polarities
are trialled. K occupied waves and n positions give the inclusive upper bound
`K*1024 - 1 - (n-1)*correction`. Derived-boundary mode BIC-encodes all positions
in [0,upper]. Explicit-boundary mode gamma-encodes first+1 and upper-last+1,
then BIC-encodes only interior positions. All fields share one bit accumulator
and one final flush. Four complete records are measured; ties prefer occupied
polarity and derived bounds. Dry-trial mode does not emit a new opcode. The
emitting experiment below uses opcode 90 with this exact wire model.

The report compares each mode separately and their per-record minimum against
actual emitted opcode-68 bytes. It includes metadata/map/BIC/padding costs,
projected savings with ordinary fallback, forced net savings, chosen modes,
empty-wave location, and log2 buckets of empty waves before the last occupied
wave. Trailing-only emptiness is kept separate because it shifts no positions.
Adaptive opcode-92 winners remain outside this trial, as in the other probes.

Build separately with `make neon`, then run the full sorted collection:

```bash
bash experiments/sblock_probe/run_gov2_sblock_probe.sh \
  --sblock-wave-trial --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --datasets sorted --threads 10 --progress 10000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-sorted-sblock-wave-$(date +%Y%m%d-%H%M%S)"
```

Read `sorted/sblock-wave-trial-report.txt`. The full harness verifies raw IDs,
decodes the ordinary emitted stream, and checks probe accounting against the
opcode profiler. Candidate sizes are projections, not a newly encoded collection.
The focused `test_sblock_wave.cpp` independently decodes all four scratch
formats, checks adversarial and random lists, and checks unchanged serializer
bytes. Example focused test command from the repository root:

```bash
c++ -std=c++17 -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer \
  -I src utils/inv_list/experiments/sblock_probe/test_sblock_wave.cpp \
  -o /tmp/test_sblock_wave
/tmp/test_sblock_wave
```

Deferred map alternatives: grouped/hierarchical occupancy flags, gamma-coded
occupancy runs, and GAP encoding. The 16,384 wave positions fit 16-bit GAP
coordinates; reuse would need to avoid a nested record or extra accumulator
flush. These alternatives are not part of this trial.

Prior full sorted run-endpoint trial (20261030-132935): 86,462 eligible
opcode-68 records, 18,840 wins, 443,772 bytes projected savings (0.01657% of
collection), versus 23,104,976 bytes growth when forced on every eligible
record. This does not justify adding a run-endpoint opcode.

## Emitted wave opcode experiment

Opt-in `-sblock-wave` enables opcode 90 for a complete wave-compacted record
only if it is strictly smaller than the ordinary opcode-68 record and any
adaptive opcode-92 candidate. Ties retain the existing codec. The four
polarity/boundary modes share the dry-trial bitstream, and the decoder
reconstructs the original sorted positions before setting the superblock.
Historical opcode 90 was a retired experiment; old BM opcodes and their
readers remain unchanged. Run this experiment with a fresh output path on
USB_1TB. The runner makes full baseline and candidate streams in sequence,
verifies each against raw IDs, decodes each, checks stream bytes and opcode
counts, then reports size and elapsed times. It deletes only verified `.bv`
files from its own run unless `--keep-streams` is given.

```bash
bash experiments/sblock_probe/run_gov2_sblock_wave_encode.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --threads 10 --progress 10000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-sorted-sblock-wave-op90-$(date +%Y%m%d-%H%M%S)"
```

Use `--max-vectors 10000` for a preliminary sample. The final comparison is
in `report.txt`; each phase log stays under `baseline/` or `wave/`.
Use `--datasets unsorted` to run the same paired experiment on Gov2.unsorted.

## Occupancy-map follow-up

The final-opcode-68/90/92 wave-size and map diagnostic is documented in
[../sblock_map/README.md](../sblock_map/README.md). It compares current BIC,
minimum-delta BIC, hierarchical 8/32/128 flag groups and gamma occupancy runs
across 128/256/512/1024/2048/4096-bit waves, with separate size-only,
map-only and joint arms. It accounts for metadata/padding without changing emitted records.
