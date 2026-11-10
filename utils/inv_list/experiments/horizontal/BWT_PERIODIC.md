# Whole-block BWT followed by periodic step 1

Build `bminv` separately, then run the first 1K Gov2.sorted vectors:

```bash
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make
caffeinate -i bash experiments/horizontal/run_gov2_bwt_periodic.sh \
  --max-vectors 1000 --threads 10 --batch-vectors 512 --progress 25
```

For the same first 10K vectors, use `--max-vectors 10000`. Each run creates a
new timestamped output directory under
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/`.
The input and source stay in the BitMagic repository. The script checks the USB
mount and rejects an existing output directory. `--out DIR` selects another
nonexistent output directory.

The harness produces three collections from the same input: ordinary,
independent three-way selection, and four-way selection with BWT followed by
periodic step 1. The fourth candidate reuses the BWT-transformed block and its
primary index. The candidate record contains the periodic residual, its formula,
and the BWT primary index. Decoding restores the transformed block from the
periodic formula and then inverts BWT. All candidates are compared by their
complete serialized record size; the previous three-way winner wins ties.
GAP blocks retain ordinary encoding.

Periodic settings match the previous experiment: independent order, step 1,
full sample, minimum chain length 4, one hole permitted, 256 search and formula
chains, trial budget 2304. The experiment is opt-in via `-block-select-compose`.
The previous `-block-select` path remains available.

The harness runs compression, verification and decode with 10 workers by default.
Compression writes a small `.idx` sidecar for each collection with exact BLOB
boundaries. The read workers use it to seek to their assigned vectors and reject
any BLOB that does not consume its indexed byte range. The sidecar is required
for threaded read commands and is checked against the BLOB size. Existing BLOB
bytes are unchanged. The terminal displays the latest progress on one line and
hides per-block diagnostic rows; complete rows remain in the phase logs.

The harness verifies all three collections against the raw input, times decoding
and writes `report.txt`, `metrics.tsv`, per-block `blocks.tsv`, the BLOBs, logs and
timing files. The report reconciles each block's four candidate sizes and ties
with both collection byte differences. It separates composed wins by whether BWT
alone helped and reports the extra periodic search time. The worker times include
ordinary evaluation inside the original search paths; the transformed periodic
time covers the added search on the already transformed block.

A short focused round-trip test (including 32- and 64-bit builds under
AddressSanitizer) can be built and run separately:

```bash
bash experiments/horizontal/build_block_selection_tests.sh ../.. /tmp/bm-block-selection-tests
/tmp/bm-block-selection-tests/block-selection-test32
/tmp/bm-block-selection-tests/block-selection-test64
```
