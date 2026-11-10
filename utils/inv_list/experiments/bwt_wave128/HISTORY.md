# Experimental BWT wave packing (opcode 95)

Build separately in `utils/inv_list` with `make neon`. The paired runner
compares the existing level-6 step-2/BWT selector, 16-byte BWT cutoff,
GAP-XOR, adaptive and wave superblocks, and opcode-65 endpoint selection
against the identical configuration with `-bwt-wave128` enabled.

The new option normalizes a BIT block's polarity, applies whole-block BWT,
records occupied 128-bit waves, and packs them stably. Its residual BIC is
bounded by `128 * occupied_waves`. It trials the occupied and empty map
polarities, then uses the existing periodic chain finder for steps 1 and 2.
Chain exceptions remain in the XOR residual. The opcode, wave map, residual,
and chain formula are one record with a single bit accumulator. The decoder
undoes chains, scatters waves, inverts BWT, and restores polarity. Opcode 95
is experimental; the existing opcode 72 decoder is unchanged.

Run from `utils/inv_list`:

```bash
bash experiments/bwt_wave128/run_gov2_bwt_wave128.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-31/runs/inv_list/gov2-sorted-bwt-wave128-$(date +%Y%m%d-%H%M%S)" \
  --datasets sorted --threads 10 --progress 1000
```

`--max-vectors N` permits a sampled prefix; omit it for the full collection.
The script requires a mounted USB_1TB drive and at least 4 GiB free before
each compression. It writes the paired reports and logs to the output directory,
verifies each stream against raw IDs, independently decodes both, and removes
only its own verified `.bv` streams by default. Use `--keep-streams` if the
serialized output needs to be retained.

The final `report.txt` checks that actual collection savings equal the encoder's
selected-record savings and that opcode-95 field bits match its emitted bytes.
It shows codec selections, metadata/residual/formula costs, chain contribution,
and how many winning blocks the proposed empty-wave gate would have skipped.
The gate is measured but does not reject candidates.

## Framing / packing / chain ablation

The follow-up runner compares four complete collections under the same best
baseline options:

1. Existing selector without a new BWT opcode.
2. Opcode 96: normalized whole-block BWT with the same one-accumulator BIC
   residual framing, but no wave map or chain search.
3. Opcode 95: packed 128-bit waves and the bounded BIC residual, without chain
   search.
4. Opcode 95: packed waves with the existing step-1/step-2 chain search.

Opcode 95's wire layout and decoder remain compatible with the preceding
experiment. Opcode 96 is a separate experimental control. Every candidate
still faces complete-record comparison with the same 16-byte BWT cutoff.
The output report reconciles stream savings and opcode field bits; all four
streams must verify against raw IDs and decode independently.

```bash
bash experiments/bwt_wave128/run_gov2_bwt_wave_ablation.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-31/runs/inv_list/gov2-sorted-bwt-wave-ablation-$(date +%Y%m%d-%H%M%S)" \
  --datasets sorted --threads 10 --progress 1000
```

The full sorted ablation verified and decoded all four streams (1,107,205
vectors, 5,979,715,441 IDs). Baseline was 2,664,780,074 bytes; unpacked was
2,662,268,168; packed without chains was 2,662,779,156; packed with chains was
2,662,754,984. Unpacked won as an addition to the baseline, which still retained
opcode 72. It has not been established as a replacement for opcode 72.
The report is in `gov2-sorted-bwt-wave-ablation-20261031-122435` on USB_1TB.

## Shared BWT policy

The earlier comparison computed original-polarity BWT for opcode 72 and a
second, normalized transform for dense blocks in the opcode-95/96 candidate.
The full sorted polarity comparison found original-polarity 96 smaller by
252 bytes. The shared selector therefore uses original polarity throughout:
one whole-block BWT supplies opcode 72, direct-BIC opcode 96, and optionally
packed/chain opcode 95. Original blocks remain available to non-BWT candidates.

This policy applies to the new shared selector. Previous experimental controls
remain available for reproducing the earlier comparison. Existing opcode
layouts and decoders, including streams carrying an inversion flag, retain
their meanings.

## Opcode-72 replacement and polarity experiment

Run `run_gov2_bwt_compare.sh` to compare five complete streams:

| Mode | Candidates |
| --- | --- |
| baseline | Existing best selector, including opcode 72 |
| normalized | Baseline plus opcode 96 with pre-BWT polarity normalization |
| raw | Baseline plus opcode 96 using original input polarity |
| without72 | Normalized selection with every opcode-72 frame excluded |
| probe | Normalized selection plus scratch comparisons of raw 96, normalized 96 and packed/chains 95 |

All passes use adaptive/wave superblocks, GAP-XOR, opcode-65 endpoint selection,
and the established BWT/step-2 selector with a 16-byte cutoff. Each pass verifies
against the raw IDs and decodes independently. The report requires the probe
and normalized streams to have identical SHA-256 hashes. Streams are removed
only after verification and accounting pass; hashes, indices, logs and reports
remain. Build is separate:

```bash
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make neon
bash experiments/bwt_wave128/run_gov2_bwt_compare.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-31/runs/inv_list/gov2-sorted-bwt-compare-$(date +%Y%m%d-%H%M%S)" \
  --datasets sorted --threads 10 --progress 1000
```

Omit `--max-vectors` for the full collection; use `--max-vectors 1000` for a
prefix. `--threads 1` is available for controlled serial timing. The harness
requires the USB mount and checks available space before every compression.

The final table includes absolute bytes, weighted bits per original ID,
encode/verify/decode times and opcode-72/96 counts. It reports:

- Actual stream-size cost of removing all opcode-72 frames.
- Net collection-size effect of normalization versus original polarity.
- Direct complete-record 72-versus-95/96 comparisons: strict wins, ties,
  unavailable packed candidates, and otherwise-smaller candidates rejected by
  the BWT cutoff. Payload opcode and post-marker distributions identify which
  representations survive.
- Separate BWT-alone, BWT-plus-chain and chain-only counts. Opcode 72 is also a
  wrapper for periodic formulas without BWT, so its top-level count is not a
  BWT count. Both baseline winners and winners retained after normalized 96
  are reported.
- Packed95-versus-normalized96 comparisons and the estimated complementary
  savings from adding packed95 to normalized selection. This last estimate is
  from scratch trials; no combined 72/95/96 stream is emitted in this experiment.

The exclusion pass still computes baseline candidates, then removes opcode 72
from final selection and falls back to ordinary encoding or opcode 96. It keeps
the original non-BWT record as the 16-byte cutoff comparator, even when that
record itself uses a chain-only opcode-72 frame. This isolates removal under
the existing acceptance policy. Its runtime is not an estimate of a cleaned-up
encoder with opcode 72 retired. Probe runtime includes extra trial work and is
also reported separately from ordinary encoding timings.

The comparison implementation is temporary: `bmbwtcompare.h` contains counters,
serializer controls gate the trial/exclusion paths, and `inv_list` prints their
aggregates. No decoder or opcode layout is changed. Transform consolidation and
a common polarity policy remain deferred until the results justify retention.

Validation for the staged comparison: NEON build and focused ASan/UBSan tests
passed (dense/sparse round trips, diagnostic byte identity, opcode-72 exclusion,
and existing packed/unpacked direct and forced-stream round trips). The final
two-vector smoke passed all five compress/verify/decode/report phases, including
SHA-256 identity and exact field accounting. This is focused validation, not a
full-collection result. Its report is under
`gov2-bwt-compare-smoke-final-20261031` on USB_1TB.

## Shared original-polarity selection: sorted and unsorted

`run_gov2_bwt_shared.sh` runs three configurations on each requested dataset:

1. `baseline`: existing best selector, adaptive and wave superblocks, GAP-XOR,
   opcode-65 endpoint choice, and the 16-byte BWT cutoff.
2. `shared`: baseline plus opcode 96, using the same original-polarity BWT as
   opcode 72 and its BWT/chain candidates.
3. `packed`: shared selection plus opcode 95, also using the original-polarity
   transform. This is a verified combined stream, measuring whether 95 merits
   retention on either dataset before removing that experimental codec.

The direct-BIC path reads the shared BWT output without wave analysis, polarity
normalization, or copying the whole block. Packed/chain trials use scratch
storage. Actual transform-call counters must equal eligible block counts in
both shared modes. Decoders and wire layouts are unchanged. The preceding raw
control also reused the original BWT, so comparisons against it should not be
interpreted as elimination of a transform that it already shared.

```bash
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make neon
bash experiments/bwt_wave128/run_gov2_bwt_shared.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-31/runs/inv_list/gov2-bwt-shared-$(date +%Y%m%d-%H%M%S)" \
  --datasets both --threads 10 --progress 1000 \
  --reference /Volumes/USB_1TB/BitMagic-experiments/2026-10-31/runs/inv_list/gov2-sorted-bwt-compare-20261031-135811
```

No `--max-vectors` means full collections. For a short prefix add
`--max-vectors 1000`; use `--threads 1` when collecting controlled serial timing.
`--reference` is optional: where input counts match, the final report requires
baseline and shared output hashes to match the previous baseline and raw
streams. Missing datasets/different prefixes are reported explicitly.

Every stream verifies against raw IDs and independently decodes. The final
summary reports absolute bytes, weighted bits/ID, encoding/verification/decoding
times, selections for 72/96/95, and the measured incremental saving from 95.
Raw profiling tables stay in the logs. One formatted opcode breakdown per
stream is printed to the terminal; intermediate saved reports are not replayed.
`Opcode ID` identifies the format; `Marker bytes` is the accumulated cost of its
one-byte top-level markers, equal to its record count. `Payload bytes` includes
metadata, embedded markers, formulas, encoded arrays and padding. `Total bytes`
is marker plus payload bytes, with collection headers accounted separately.

The USB mount/free-space checks and automatic deletion of successfully verified
streams remain in place. Logs, reports, indices and SHA-256 files are retained.

Staging validation passed: NEON build; focused ASan/UBSan dense/sparse tests;
six compress/verify/decode smoke streams (two vectors from each dataset).
Shared and packed modes each made exactly 517 BWT calls for 517 eligible sorted
blocks and 770 calls for 770 eligible unsorted blocks. Sorted baseline and
shared hashes matched the preceding baseline/raw smoke streams. Console checks
found exactly six formatted breakdowns and zero raw duplicate tables. The
full sorted/unsorted collection experiment remains for the user to run.

## Wave-size due diligence (scratch diagnostic)

`run_gov2_bwt_waves.sh` defaults to 10,000 vectors from both datasets. It emits
shared original-polarity 72/96 baseline streams, then identical streams with
six wave-size trials (32, 64, 128, 256, 512, 1024 bits). Wave candidates have no
chain search; existing baseline chain choices remain enabled. Both map
polarities are trialled. Fixed-size records use wave-dependent occupied-count
width and BIC bounds; adaptive choice adds three selector bits before the final
32-bit flush. Its size is calculated from the exact unpadded field lengths,
with focused tests comparing this result against an actual selector write.
The fixed 128-bit scratch record is byte-identical to the existing opcode-95
writer without chains. No new-size records are emitted or added to the decoder.

```bash
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make neon
bash experiments/bwt_wave128/run_gov2_bwt_waves.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-31/runs/inv_list/gov2-bwt-waves-10k-$(date +%Y%m%d-%H%M%S)" \
  --datasets both --threads 10 --max-vectors 10000 --progress 1000
```

The report requires identical baseline/probe hashes, raw-ID verification,
independent decoding, one BWT call per eligible block, and exact accounting of
accepted fixed-candidate fields. Fixed and adaptive gains are calculated
complete-record opportunities against the shared baseline, with the existing
16-byte cutoff against the non-BWT alternative. The new formats would still
need emitted-stream validation if retained. Trial times per size are accumulated
worker time (including mapping, scanning and encoding), not wall-time estimates
for six separate encoders.

The completed paired 10K measurements and interpretation are in
[WAVE_SIZE_FINDINGS.md](WAVE_SIZE_FINDINGS.md). Fixed 1024-bit waves gave the
largest sorted gain (15,266 bytes, 0.00105%); unsorted gains were negligible.
