# Current highest-compression Gov2 verification

Build separately with `make neon` in utils/inv_list, then run:

```bash
bash experiments/current/run_gov2_highest.sh \
  --bin ./bminv --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --datasets both --threads 10 --progress 10000 \
  --level 10 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-11-03/runs/inv_list/gov2-level10-$(date +%Y%m%d-%H%M%S)"
```

No max-vectors argument means full collections. For a short run add
`--max-vectors 10000` and omit `--reference` (the full reference has a different scope).

The harness now uses the library API preset `set_compression_level(10)`;
codec flags are not assembled by the harness. Level 10 enables retained
non-BWT selectors and shared whole-block BWT/chain competition, accepting
one byte of complete-record BWT savings. Level 9 uses a 16-byte threshold;
level 8 excludes BWT/chains, and level 7 selects only released formats.
The default remains level 6. See the in-place serializer Doxygen table.
Use `--bwt-min-gain N` with level 9/10 for a custom BWT threshold. The inverse
search cap defaults to 196,607 zero positions. Reapplying a library level
resets its preset; apply custom option setters afterwards.

The previous level-6-plus-explicit-options reference used a 16-byte threshold,
so do not request byte identity for the new level-10 run. To reproduce its
policy with the API, use `--level 9` and a matching `--reference`.

Experiment archives and new output live under
`/Volumes/USB_1TB/BitMagic-experiments`. The default input directory is
`/Volumes/USB_1TB/int_set_bench/Gov2Flat`; sources and scripts stay in this
repository. The previous DATAFAT32 archive was migrated with per-file SHA-256
verification; its receipt is `migration-from-DATAFAT32-20261103.json` at the
new archive root. Historical run logs retain the paths used at execution time.

Progress is emitted after `--progress N` vectors or 30 seconds, whichever
comes first, including during ordered worker waits in compression, verification,
and decoding. `vectors` counts records written/checked; `prepared` counts
dispatched jobs, and `completed` includes later jobs finished ahead of the
ordered writer. Terminal output replaces one row; captured output prints each
update immediately. Deadlines are best effort during blocking filesystem I/O.
Direct bminv use without `-progress-vectors` retains disabled progress.

Focused reporting checks (no Gov2 input required), from utils/inv_list:

```bash
c++ -std=c++17 -O2 -I../../src experiments/current/test_progress.cpp -o /tmp/bm-test-progress
/tmp/bm-test-progress
# Optional wall-clock check: one slow job, report at 30 seconds.
/tmp/bm-test-progress --realtime
python3 experiments/current/test_progress_console.py
```

The harness checks the USB mount/free space, compresses both datasets,
verifies against original raw IDs, independently decodes and profiles all
emitted records. Ten workers are used for compression and read checks.
An optional `--reference VERIFIED_RUN_ROOT` requires identical vectors, IDs,
byte count and SHA-256; its inputs/options must match the current run.

`verification-report.txt` prints absolute bytes, weighted bits/ID and phase
timings. `report.txt` and per-set `profile/report.txt` retain the detailed
opcode breakdown. Generic stream validation has no dependency on codec-specific
statistics lines.

## November 3 cleanup

Removed dedicated counters/getters/aggregation for GAP-XOR, shared BWT,
adaptive/fixed-wave/high-count superblocks and opcode-68 correction trials.
Their experimental reports are removed from inv_list. Generic opcode byte/ID
profiling, progress, stream reconciliation and internal dry-score/write-size
checks remain. `--diagnostics` now adds worker timing diagnostics. Historical
comparison harnesses use verified complete-stream sizes and generic profiles.

The BM deserializer normalizes restored blocks when the destination is empty
at decode entry, including positive/inverse superblocks. Finalization happens
after transforms and pending reference XOR; uniform blocks are canonicalized,
GAP allocations can shrink, and eligible BIT blocks become GAPs. Nonempty
OR destinations retain their existing merge behavior. Memory and streaming
readers share the implementation; serialized formats have not changed.

Focused regression source: `tests/regression/serial_gap_normalization.cpp` (at repository
root). It checks block representation as well as bits, sparse/dense/solid/BIT
cases, opcodes 68/70/72/90/92/94/96/97/98, both opcode-98 polarities, reference
XOR, empty-only behavior, reader reuse, and streaming/index/range reads.

Generated .bv streams are deleted after successful verification/report;
logs, indices and hashes remain. `--keep-streams` retains streams. Inputs and
sources stay on the system drive. Full-collection validation is the next user run.

## Opcode registry and decoder retirement (November 3)

All released opcode IDs 0..69 and all 20 retained new assignments are unchanged.
New assignments are in numerical order with compile-time checks for duplicates,
legacy collisions and the ordinary-opcode limit of 128. Vacant values are
82..83, 86..89, 91, 93, 95 and 99..127 (38 slots). Bytes 128..255 belong to
packed zero-run framing and are not available for ordinary records.

Optional writers still use horizontal/per-wave BWT and the older grouped-chain
formula markers. Those formats remain active. The retired decoder-only parts
are opcode 92's earlier fixed-step header/global-relative predictor and opcode
96's earlier inverse-polarity/chain-formula variants. Their unreleased streams
now require the older binary. Retained current writer output is byte-identical.
Opcode 96 keeps its fixed zero header fields and empty-formula marker to avoid
changing emitted bytes. Released decoder paths remain present.

The full command above uses the November 3 verified run as a byte-identity
reference, checking full sorted and unsorted sizes and SHA-256 in addition to
raw verification and independent decode.

Focused checks added for this retirement:

- `tests/regression/serial_opcode_registry.cpp`: freezes the 70 released assignments and
  checks all retained/vacant slots.
- `tests/regression/serial_released_compat.cpp`: generates fixtures against a released
  checkout with `BM_WRITE_RELEASE_FIXTURES`, then reads those fixtures against
  the current headers. The local v9.3.1 test covered 184 files, levels 0..6,
  bookmarks, reference XOR, empty destinations, OR merging and bounded streams.
- Existing adaptive-delta and shared-BWT tests now reject the removed variants.

Native NEON build, ASan/UBSan codec and 32-/64-bit GAP-normalization checks,
released-stream compatibility, and the 10-worker pre/post byte-identity smoke
passed. Full-collection validation of this retirement is the next user run.

Focused CLI/configuration cleanup checks (from utils/inv_list):

```bash
python3 experiments/current/test_cli_cleanup.py --bin ./bminv
# Optional: compare every fixture stream against a saved pre-cleanup binary.
python3 experiments/current/test_cli_cleanup.py --bin ./bminv --reference-bin /path/to/previous/bminv
```

The check exercises levels 4 and 6–10, serial and 1/10-worker encoding,
verification/decode, final-position mode switches, and invalid arguments.
The retired `-svin`/`-svout` options fail explicitly; use the separate
`experiments/sparse_vector/svexp.cpp` utility for sparse-vector experiments.

Pair-search refactoring parity check (small synthetic collection; no Gov2 I/O):

```bash
python3 experiments/current/test_xbv_pair_search.py --bin ./bminv --reference-bin /path/to/previous/bminv
```

It compares candidate, best and measured CSVs plus summaries (excluding elapsed
time), with inspection on/off, ratio filtering, prefixes, empty results and
one/four evaluation workers.

Recover an existing profiled compression after a verify/decode/report harness
failure without repeating compression:

```bash
bash experiments/current/recover_gov2_profile.sh --run EXISTING_RUN \
  --dataset sorted --bin ./bminv --threads 8 --progress 1000
```

The recovery script reads the recorded vector limit, preserves prior verify/
decode logs and timings, verifies and decodes the saved stream, and regenerates
the profile report and SHA-256. It retains the compressed stream and original
failure marker. The temporary `horizontal/progress_console.py` relocation shim was retired;
all supported harnesses use `current/progress_console.py`.

## Portable collection checks

`../../run_collection.sh` (relative to this directory) replaces the old top-level
Gov2/ClueWeb/ARM runners. From utils/inv_list:

```bash
bash run_collection.sh --bin ./bminv --input /data/gov2.sorted \
  --out /results/gov2-sorted --levels 6,7,8,9,10 --threads 10
python3 experiments/current/test_collection_runner.py --bin ./bminv
```

The output directory must be new, with an existing parent. This runner preserves
streams and performs raw verification and independent decoding per preset. The
profiled USB workflow above remains supported. Old standalone BWT/periodic and
compose/multistep CLI modes are retired; see ../RETIREMENT_CANDIDATES.md.
