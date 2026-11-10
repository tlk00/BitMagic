# Independent ordinary / periodic step 1 / whole-block BWT experiment

## Formula-scoring policy trial (2026-10-27)

The first 10,000 Gov2.sorted vectors were encoded with one compression worker
and the current level-6, step-2 block selector (16-byte BWT threshold,
gamma-versus-fixed step-1 formula choice). Four temporary implementations
produced byte-identical 1,466,801,846-byte streams and indexes; all verified
against the input and decoded. Compression wall times were 4257.27 seconds
for the existing scorer, 4240.77 for bounded gamma-first, 4247.37 for bounded
fixed-first, and 4236.37 for bounded gamma-first with deferred emission.
The last variant skipped 18,220,891 of 27,587,059 formula writes, but saved
only 20.90 seconds (0.49%) against the control. This single-run difference is
too small to justify maintaining the extra scoring and emission paths. The
existing scorer remains; the temporary policies and CLI option were removed.
Full phase logs and the timing table are on the USB under
`gov2-sorted-formula-policy-20261026-170514`.

Build separately, then run from `BitMagic/utils/inv_list`:

```bash
make
caffeinate -i bash experiments/horizontal/run_gov2_block_select.sh \
  --max-vectors 10000 --threads 8 --batch-vectors 512 --progress 25
```

The default input is `/Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.sorted`. The timestamped
output directory is under
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/`.
The harness checks that the USB volume is mounted before creating output and
before each phase. Both serialized collections are retained. Use `--out DIR`
to choose a different, nonexistent destination, or `--max-vectors 1000` for a
shorter run. Sources and input remain on the system drive.

The harness compresses the same vectors twice: ordinary baseline, then
`-block-select`. Both use the requested thread count. It verifies both against
the original collection and times complete decoding separately. Progress is
reported at batch completion; `--progress` does not reduce the batch size.

For each eligible BIT block, the selector runs periodic step 1 and whole-block
BWT independently on the original block. It compares complete serialized
records, including formula, wrapper and padding bytes, against ordinary
serialization. Ordinary wins ties; periodic wins equal profitable candidates.
GAP blocks retain ordinary encoding. BWT is evaluated even if periodic wins or
BWT increases GAP count. There is no BWT-plus-periodic composition.

The fixed periodic settings match the previous step-1 experiment: independent
order, steps 1..1, sample 1, minimum length 4, at most one hole, 256 search and
formula chains, and trial budget 2304. BWT transforms the whole 65,536-bit block,
with sample 1. No additional selector opcode is needed: the selected record uses
its existing opcode. Search is disabled by default and mutually exclusive with
other experiment modes. The direct CLI option is `-block-select`, which now
uses a 16-byte minimum BWT gain; this original harness passes
`-bwt-min-gain 0` to reproduce unrestricted selection.

Outputs:

- `report.txt`: actual collection savings, compression/verification/decode wall
  times, chosen representation counts and candidate evaluation times.
- `blocks.tsv`: vector/block coordinates, complete candidate sizes, selected
  kind (0 ordinary, 1 periodic, 2 BWT), trials and raw BWT framing/payload sizes.
- `metrics.tsv`: summary values suitable for comparison between runs.
- `run.txt`, phase logs and timing files: commands and detailed diagnostics.
- `baseline.bv`, `selection.bv`: verified serialized collections.

Periodic and BWT candidate sizes in the comparison include each path's ordinary
fallback. Their independent savings are diagnostic estimates from the candidate
records, not extra output collections. The report checks every block's minimum,
tie breaking, BWT framing, aggregate counters and collection byte difference.
Each search path currently performs its own ordinary encoding; worker timing
includes that cost. Summed worker elapsed time is not wall time or CPU time.

## Focused regression tests

```bash
bash experiments/horizontal/build_block_selection_tests.sh ../.. /tmp/bm-block-selection-tests
/tmp/bm-block-selection-tests/block-selection-test32
/tmp/bm-block-selection-tests/block-selection-test64
```

These opt-in AddressSanitizer tests cover empty/full, sparse GAP, dense,
structured, random and multi-block vectors; round trips for the independent,
composed, multistep, step-2 and thresholded selection modes;
minimum size, tie breaking and savings accounting. They do not enable the
experiment in the library stress suite.

## BWT gain threshold design note

The experimental `-bwt-min-gain N` option requires a BWT-bearing winner to
beat the best non-BWT record by at least `N` complete serialized bytes. A
rejected BWT record falls back to the better of ordinary and periodic step 1.
The opt-in selector now defaults to 16 bytes; `-bwt-min-gain 0` explicitly
restores unrestricted size-first selection. Ordinary serialization remains
unchanged unless a block-selection mode is enabled. The older comparison
harnesses pass zero explicitly to preserve their unrestricted controls.
Future compression levels may choose different thresholds once measured
decode cost and size loss justify them. The unsorted Gov2 distribution is a
useful control for setting a threshold because its BWT gains are often small.

## Step-1 chain-start encoding experiments

`-periodic-step1-start-reduce` is opt-in. For a step-1 group of `n` sorted,
nonoverlapping chains with minimum length four, it encodes start `i` as
`start[i] - 4*i` and narrows the BIC upper bound to `65536 - 4*n`.
The decoder adds `4*i` back. Other steps retain the original start encoding.
New post-block markers distinguish this formula from existing serialized
records, so the decoder accepts both versions. The existing formula is emitted
unchanged when the option is absent. Complete encoded block sizes, including
formula and padding, still control candidate selection and ordinary fallback.

The next experiment is to emit the first start explicitly as metadata and use
it as the lower boundary for BIC of the remaining reduced starts. Compare the
metadata cost plus bit-writer padding with the current shared-accumulator BIC
using complete-record bytes. Further candidates to measure are a
stronger final-start upper bound derived from the full lengths, a dynamic
minimum separation, and a conditional formula mode when the fixed reduction
loses. Chain-count and first-start bounds must be reconstructed identically by
the decoder. Run sorted and unsorted distributions before promoting any mode.

From `utils/inv_list`, build with `make neon`, then run
`experiments/horizontal/run_gov2_step1_start_reduce.sh` for a two-stream
comparison. Its default is the complete Gov2.sorted collection with ten
workers, level 6, independent BWT-residual step 2, and a 16-byte BWT gain
threshold. Outputs remain on the mounted USB_1TB drive.

The opt-in `-periodic-step1-first-start` experiment implements the first-start
variant with a separate marker. A one-bit mode chooses gamma(start+1) or a
16-bit literal, whichever takes fewer bits; later starts retain the four-bit
minimum-length reduction and use first+1 as their BIC lower bound. It shares
the formula's bit accumulator and terminal flush. The search still compares
complete records with ordinary and BWT candidates; this experiment compares
the resulting stream against the reduced-start mode rather than silently
mixing both formulas in one run.

`run_gov2_step1_first_start.sh` can build a fresh two-stream comparison with
`--max-vectors N`. For the full sorted collection, pass `--control-dir` pointing
to the verified `run_gov2_step1_start_reduce.sh` output. The runner reuses its
`reduced.bv` and logs read-only, writes only the new stream to the USB drive,
verifies and decodes it, and checks complete-record size accounting.

The separate BIC-ends formula was retired after its full Gov2.sorted run lost
33,361 bytes against explicit-first gamma lengths. Its experimental opcodes
82 and 83 are no longer decoded; preserve the older binary if those streams
must be read.

The opt-in `-periodic-step1-interleaved N` experiment BIC-encodes the
interleaved inclusive coordinates `[start0,end0,start1,end1,...]` after an
explicit first start. Mode 1 uses the guaranteed within-chain and between-chain
minima (3 and 2). Subtracting one less than each minimum leaves strictly
increasing corrected coordinates. Mode 4 compares exact formula bit costs for
explicit-first gamma lengths and fixed interleaving, including the final
32-bit flush, and emits the smaller formula. The residual payload and framing
are identical across these choices; the enclosing block search still compares
complete records with ordinary and BWT candidates. It applies only to step 1;
other steps retain the previous formula. Both modes are disabled by default.
Run `run_gov2_step1_interleaved.sh --max-vectors 1000 --mode auto` from
`utils/inv_list` after building with `make neon`. The runner keeps raw input
local, writes streams to USB_1TB, verifies and decodes each stream, and
reports complete-record size changes. For a full pass with limited USB space,
use `--mode auto --control-dir DIR` with a verified explicit-first control.
The measured-minimum shared and pair formats lost the 10K Gov2.sorted
comparison and were removed, including their decoders. Their experimental
opcodes 86 through 89 remain reserved; use the older binary to read old streams
containing them. See `STEP1_INTERLEAVED_10K.md` for the measurements.

To measure the threshold on emitted streams, run the separate
`run_gov2_bwt_threshold.sh` harness. It compares ordinary, unrestricted
step-2 selection, and thresholded step-2 selection; verifies and decodes all
three; and reports BWT-bearing block counts, extra bytes, and decode times.
Omitting `--max-vectors` runs the entire input collection; pass it only for a
prefix sample. The report uses the actual completed vector count.
This experiment does not change the default production compression levels.

For a decode-only timing comparison of the completed full Gov2.sorted
unrestricted, 16-byte, and 32-byte streams, run from `utils/inv_list`:

```bash
caffeinate -i python3 experiments/horizontal/run_gov2_bwt_decode_compare.py
```

The script checks the existing streams and indexes, warms each once, then
rotates their order across six measured rounds. It writes per-run logs,
wall/worker times, and paired differences to a new small USB output directory.
It does not recompress or produce another BLOB collection. Use `--run16` and
`--run32` if the source run directories change.
