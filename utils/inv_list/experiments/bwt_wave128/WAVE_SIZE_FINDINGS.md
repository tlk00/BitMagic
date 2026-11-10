# BWT wave-size due diligence, 10K sorted and unsorted Gov2

Run: `gov2-bwt-waves-10k-20261031-174949` under
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-31/runs/inv_list/`.

## Scope and method

The first 10,000 vectors of each dataset are used, not a random sample. Each
prefix contains 4,818,820,706 IDs. The baseline retains the current shared
original-polarity BWT selector (72/96), chain choices, superblock choices,
GAP-XOR, opcode-65 endpoint selection, and the 16-byte BWT cutoff. Compression,
verification, and decoding use 10 workers.

The diagnostic reuses the same single BWT result per eligible block and trials
32, 64, 128, 256, 512, and 1024-bit waves. It stably packs occupied waves, tries
both empty and occupied BIC maps, and BIC-encodes the compacted set-bit list.
There is no additional chain search for these wave candidates. Any block with
both empty and occupied waves can compete; there is no minimum-three-waves gate.

Fixed-width candidates charge the entire record, including its marker,
primary index, occupancy metadata, map, residual, empty-chain marker and final
32-bit accumulator padding. Adaptive selection additionally charges a 3-bit
wave-width selector and recalculates padding. Candidates must beat the selected
baseline record and clear the existing 16-byte threshold against the non-BWT
alternative. Fixed-width gains are alternatives, not additive.

These are scratch trial records. The probe emits the unchanged baseline
stream, and the harness checks matching SHA-256 hashes as well as raw-ID
verification and independent decoding. Gains therefore describe calculated
complete-record opportunities, not validated variable-wave output streams.

## Validation

Before the dataset runs, native ARM NEON build and focused ASan/UBSan tests
passed. Tests cover scratch map/position round trips for all six widths and
both map polarities, actual versus calculated adaptive-selector record sizes,
128-bit scratch byte identity with the existing chain-free opcode-95 writer,
and unchanged serializer output with probing enabled. A two-vector paired
smoke run on both datasets passed stream verification, decoding and hashes.

## Results

All four 10K streams completed successfully. Baseline/probe SHA-256 hashes
match on each dataset; raw-ID verification, independent decoding, one-BWT-per-
eligible-block accounting, and complete fixed-candidate field accounting pass.

Sorted baseline: **1,457,318,527 bytes**, or **2.419378 bits/ID**. The probe
stream has the same SHA-256 hash. Both streams passed raw-ID verification and
independent decoding. All 634,155 eligible blocks were probed and each used
exactly one BWT transform. Selected-candidate field totals reconcile exactly
with their complete record sizes.

| Wave bits | Sorted winning blocks | Sorted saved bytes | Sorted gain | Unsorted winning blocks | Unsorted saved bytes |
|---:|---:|---:|---:|---:|---:|
| 32 | 48 | 246 | 0.000017% | 0 | 0 |
| 64 | 173 | 972 | 0.000067% | 0 | 0 |
| 128 | 518 | 2,990 | 0.000205% | 3 | 48 |
| 256 | 943 | 5,018 | 0.000344% | 2 | 32 |
| 512 | 1,507 | 6,962 | 0.000478% | 3 | 36 |
| 1024 | 3,294 | 15,266 | 0.001048% | 3 | 36 |
| Adaptive six-size selector | 2,320 | 11,500 | 0.000789% | 4 | 52 |

The best fixed width, 1024, wins on 0.5194% of eligible sorted blocks, averaging
4.63 bytes saved per winner. Its projected collection size is 1,457,303,261
bytes. Adaptive selection is worse than that fixed choice: its 3-bit selector
can push a record across a 32-bit flush boundary, erasing a small saving.

Larger waves have fewer packable blocks (500,276 at 1024 versus 633,373 at 32),
but cheaper maps. Among winning candidates, mean map size falls from 42.29
bytes at 32 bits to 3.36 bytes at 1024. These are different sets of winning
blocks, so this is descriptive rather than a paired causal measurement.

Unsorted baseline: **3,561,956,968 bytes**, or **5.913409 bits/ID**. All
1,362,922 BWT-eligible blocks were probed. The best fixed width saves just
48 bytes (0.00000135%); adaptive selection saves 52 bytes (0.00000146%).
Projected sizes are 3,561,956,920 and 3,561,956,916 bytes respectively.
At 32 bits, 1,336,198 blocks have empty waves; at 1024 bits only 22,108 do.
Finding empty waves is plainly insufficient to predict a useful complete
record. Existing encoding and map/framing costs absorb virtually all benefit.

| Dataset / mode | Encode seconds | Verify seconds | Decode seconds |
|---|---:|---:|---:|
| Sorted baseline | 694.197 | 10.782 | 7.337 |
| Sorted six-width probe | 718.555 | 10.671 | 7.429 |
| Unsorted baseline | 957.301 | 11.850 | 7.393 |
| Unsorted six-width probe | 986.308 | 11.961 | 7.382 |

Sorted probe overhead was 24.358 seconds (3.51%). Its accumulated trial time
across workers was 217.07 seconds.
Unsorted overhead was 29.007 seconds (3.03%), with 264.58 accumulated worker
seconds in the trials.

## Interpretation and recommendation

- 128 bits was not the best fixed width for sorted: 1024 improves its saving
  by 12,276 bytes. Nevertheless, even that result is only 0.00105% of the
  collection. The sweep finds no material size opportunity on unsorted.
- Adaptive width selection is not justified by these results. On sorted,
  its selector/padding penalty costs more than its additional flexibility
  recovers compared with fixed 1024. On unsorted its total gain is 52 bytes.
- The retained 16-byte cutoff is measured against the non-BWT alternative,
  not against the already selected BWT record. Thus a wave candidate can
  legitimately contribute only 4 bytes beyond the baseline while still
  satisfying the cutoff. The report does not confuse those two comparisons.
- Keep the existing shared 72/96 path. This experiment provides no reason to
  add a production variable-wave decoder or enable adaptive wave search.
  No existing opcode was retired by this task.
- This is not an exhaustive rejection of wave-plus-chain coding. Wave
  candidates here have no new chain search, and 1024 is the largest tested
  width, not a proven global optimum. If one final targeted check is desired
  before retiring opcode 95, fixed 1024 with the existing chain search is a
  more focused next test than retaining six adaptive widths.

The previous full shared-BWT report (`gov2-bwt-shared-20261031-153826`) counted
634,586 sorted and 1,362,922 unsorted eligible blocks. These prefixes therefore
cover 99.93% and 100% of those counts respectively. That supports the relevance
of the measurement to BWT opportunities, without making this a full-collection
verification or proving identical per-block membership across runs. The prior
full 128-bit packed-with-chains result saved 16,628 sorted / 912 unsorted bytes;
it is not directly interchangeable with this chain-free width sweep.

## Artifacts and reproduction

- Full counters and size accounting: `report.txt`, `sorted.waves.json`, and
  `unsorted.waves.json` in the run directory above.
- Per-mode/per-dataset logs and hashes: `baseline/{sorted,unsorted}/profile/`
  and `probe/{sorted,unsorted}/profile/`. Verified `.bv` streams were removed
  by the harness; reports, logs, indexes, and hashes remain.
- Sorted SHA-256: `0bb703a6ca07503ef77d99b3329eb9af785bfedcc54de7c5b870487caa7d81b2`.
- Unsorted SHA-256: `d016ffc7feb7b72463dd2794a3d370f87e3588e73c3139591593e8fec2001b04`.
- Source diagnostic: `src/bmbwtwavesprobe.h`; focused test:
  `utils/inv_list/experiments/bwt_wave128/test_bwt_waves.cpp`.

From `utils/inv_list`, build separately with `make neon`, then run:

```bash
bash experiments/bwt_wave128/run_gov2_bwt_waves.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-31/runs/inv_list/gov2-bwt-waves-10k-$(date +%Y%m%d-%H%M%S)" \
  --datasets both --threads 10 --max-vectors 10000 --progress 1000
```

## Timing interpretation

Baseline/probe encode times are single-pass wall measurements with 10 workers.
Per-width trial times sum elapsed time across workers, including wave scanning,
packing, encoding, and scheduling effects. They are not six independent encoder
benchmarks. Unchanged streams should have equivalent decode work; small timing
differences do not establish a decoding speed change.
