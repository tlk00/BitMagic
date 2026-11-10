# Emitted bvector opcode sizes

Build `bminv` separately from `utils/inv_list` using `make neon`.
The optional `-profile-opcodes` switch profiles each already-serialized BLOB
before the main thread releases it. The collector in `src/bmserial_profile.h`
uses the bounded deserialization-index pass to obtain actual top-level record
boundaries. Size-only collection avoids restoring bit blocks; the inv_list
profile also reconstructs one vector at a time for exact original-ID counts. No scoring or candidate estimates enter the counts. The byte input
view does not copy the BLOB. Collection data has no printing code.

`src/bmdbg.h` provides opcode names, families and a human-readable table.
For each top-level record, count one opcode byte and every remaining byte as
payload. Payload includes codec metadata, bit-coded data and padding. Packed
zero runs are one marker byte with zero payload. Headers and ID-list BLOBs
have separate totals. Embedded residual and formula opcodes inside BWT,
horizontal or GAP-XOR frames are counted in their parent record's payload;
this is an outer-record accounting pass, not a split of BIC/gamma bit fields.
All row, family and whole-collection totals are checked against the emitted
BLOB byte lengths. Profiling adds a parsing pass and is disabled by default.

From `utils/inv_list`, run the full Gov2.sorted collection:

```bash
bash experiments/opcode_profile/run_gov2_opcode_profile.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --datasets sorted --threads 10 --progress 10000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-sorted-opcodes-$(date +%Y%m%d-%H%M%S)"
```

The harness uses level 6, `-block-select-step2`, fixed interleaved step 1,
16-byte BWT cutoff, GAP-XOR, 100% superblock grouping tolerance, and
adaptive superblock corrections. It verifies against raw IDs, independently
decodes the stream, and checks the profile's byte accounting. The generated
`.bv` is removed after successful validation unless `--keep-streams` is set;
reports, logs, index and JSON measurements remain on the USB. Raw Gov2 input
stays in the repository. `--max-vectors N` supports a limited run, with zero
meaning the full collection. `--datasets both` adds unsorted Gov2. The USB
mount and free space are checked before compression.

Focused test: `c++ -std=c++17 -O2 -DBMNEONOPT -Isrc
utils/inv_list/experiments/opcode_profile/test_opcode_profile.cpp -o /tmp/test-opcode-profile
&& /tmp/test-opcode-profile` from repository root. It covers ordinary BLOBs,
packed zero runs, legacy ID lists, repeated accumulation and a truncated
BLOB rejection. A 20-vector threaded smoke run checked full harness reporting
and byte-identical serialization with profiling enabled.

## Per-opcode bits per integer

`-profile-opcodes` now collects original set-bit cardinalities as well as byte
spans. `Bits/int` is `8 * total record bytes / represented original IDs`;
`Payload b/int` excludes only the one-byte top-level marker. Both are weighted
ratios of totals, not averages of per-record ratios. Inverted, BWT and GAP-XOR
records use restored set-bit counts; superblocks count all original IDs they
represent. Empty/framing records show `n/a`. Collection bits/ID also includes
collection-vector headers and legacy plain-ID-list bytes.

The optional collector performs a bounded full decode in addition to the index
pass. Profiling time is included in encode wall time, so these runs should not
be compared as pure encoder timings with older size-only profile runs. Ordinary
serialization and all wire formats are unchanged. Size-only library collectors
remain dry unless `set_collect_ids()` is requested. Counts reconcile against
each decoded vector and against the raw-ID verification total in the report.
Reporting lives in `bmdbg.h` and the harness; collection lives in BM headers.
The final table and `profile/measurement.json` retain absolute bytes and IDs.

Full unsorted measurement, from `utils/inv_list` (build separately):

```bash
make neon
bash experiments/current/run_gov2_highest.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --datasets unsorted --threads 10 --progress 1000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-unsorted-opcode-bpi-$(date +%Y%m%d-%H%M%S)"
```

This enables all retained highest-compression options, raw verification and
independent decode. Generated streams are removed after successful reporting;
logs, hashes, index and reports remain. Add `--keep-streams` to retain them.

Validation of this addition: focused ASan/UBSan tests cover sparse/mixed/full
blocks, bookmarks, inverted/BWT inputs, retained superblock selectors, empty
records and legacy ID-list framing. A two-vector unsorted end-to-end smoke
verified and independently decoded 25,120,402 IDs; per-opcode counts reconcile.
The 6,173,674-byte stream is byte-identical to the earlier retained two-vector
small-wave baseline. This is focused validation, not a full Gov2 run.

Full Gov2.unsorted measurement and interpretation: [results](RESULTS.md).
