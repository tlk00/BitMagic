> Historical experiment: horizontal search and BWT+XOR search have been retired.
> Their command-line options and dedicated runners are removed; opcode slots 70–71 are vacant.
> The results below are retained for reference. Other BWT and periodic experiments remain available.

# Build separately; run Gov2 with the reporting harness

`run_gov2_unsorted.sh` is a Bash script that uses an already-built `inv_list`
binary. It never invokes make or a compiler. The app remains disabled by default;
the script enables horizontal search only for the experimental run.

## 1. Build the app using your normal make configuration

```sh
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make
```

Use your usual make flags/target for the architecture you want to measure.
If your executable has a different name, supply it with `--bin` below.

## 2. Execute the Gov2 unsorted comparison

Run from `utils/inv_list`:

```sh
HARNESS=/Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list/experiments/horizontal

bash "$HARNESS/run_gov2_unsorted.sh" \
  --bin ./bminv \
  --input /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.unsorted \
  --out ./gov2-unsorted-horizontal-d2-e8-v64 \
  --max-vectors 1000 --progress 10 --diagnostics \
  --level 6 --depth 2 --effort 8 --variants 64 --functions 7
```

The result directory must not exist. For each further configuration, use a new
result directory. Depth 1 evaluates pair XOR; depth 2 adds compound predictors.
No rotation search is implemented yet. Full Gov2 runs are left to you.

| Option | Default | Meaning |
|---|---:|---|
| `--level` | 6 | Existing codec compression level, 0–6, used by both runs |
| `--max-vectors` | Full set | Compress and verify only the first N input vectors |
| `--progress` | 10 | Print live counters every N vectors; 0 disables updates |
| `--diagnostics` | Off | Report search funnel, trial/accepted graph families, and ordinary block-size buckets |
| `--depth` | 2 | Maximum predictor arguments, 1–8 |
| `--effort` | 8 | Compound-argument beam width, 1–64 |
| `--variants` | 64 | Maximum complete trial encodings per block, excluding baseline |
| `--functions` | 7 | Mask: XOR=1, XOR-of-OR=2, XOR-of-AND=4 |

Unary XOR is available in each family as the one-argument case. XOR2 means
`A ^ B ^ C`. The trial limit is a ceiling: duplicate graphs are skipped, and
this first generator has at most 252 profile/prefix slots. Search sees existing
BIT blocks reaching individual encoding; GAP inputs are never expanded and
existing superblock decisions are retained. Reference-vector XOR is excluded.

## What the script executes

1. Ordinary compression, without horizontal search.
2. Compression with the configured horizontal search.
3. Read/decode and compare the baseline collection against the original input.
4. Read/decode and compare the horizontal collection against the original input.
5. Read/decode the baseline without rereading/comparing the original collection.
6. Read/decode the horizontal output without the original-input comparison.

Both verification passes use the existing `inv_list -verify` implementation.
Compressed BLOBs need not match each other; decoded vectors must match the source.
Any failed command stops the script with its exit status and writes `FAILURE.txt`.
A stale executable without horizontal options is rejected before the long run.

## Reports

The script streams all six phases to the terminal while logging them, then prints
and saves `report.txt`, and saves machine-readable
`metrics.tsv`. It retains `baseline.bv`, `horizontal.bv`, all six command logs,
all six timing files, and `run.txt` with configuration and exact commands.
For a limited run, `report.txt` shows vectors and raw bytes actually consumed,
so throughput and bits per integer refer to the same prefix. Compression
progress shows elapsed time, vectors per second, and most recent vector time.
The horizontal phase also shows cumulative blocks, accepted transforms,
trial encodings, predictor evaluations, saved bytes, and the most recent
vector's block/evaluation counts. `--progress 1` provides per-vector timing if
the initial run warrants closer inspection.

The report includes:

- Actual file sizes, bits/integer, bytes and percentage saved, including an explicit
  **NO SAVINGS** result when sizes match.
- Wall/user/system compression time, verification wall time, and decode-only
  wall/user/system time. All six phases have full timing data in `metrics.tsv`.
- Compression throughput, compressed-input decode throughput, and time ratios.
- BIT blocks searched, accepted blocks, fallback blocks, and acceptance percentage.
- Complete trial encodings and predictor evaluations, including per-block and
  per-vector averages and evaluations per compression wall second.
- Baseline bytes attributable to searched blocks; accepted residual and formula/
  framing bytes; formula overhead percentage and net savings per accepted block.
- A cross-check that the actual file-size difference equals reported net savings.
- With `--diagnostics`, the number of heuristic candidate blocks, blocks whose
  trials saved no bytes, reached trial caps, saturated target-wave scores,
  accepted graph families, and savings by ordinary block-size bucket. The
  report checks that family and size-bucket totals match serializer totals.

The block statistics are the blocks searched, not a count of all GAP, zero, full,
or superblock records. Timings include I/O and are a single sequential run;
filesystem cache state is not controlled. Verification time includes reading the
raw input and comparing values, so use decode-only time to assess reconstruction.
Experimental BLOBs require the updated decoder.

## Larger diagnostic run with limited disk space

`run_gov2_diagnostic_stream.sh` sends both encoded collections to `/dev/null`
while retaining their serializer byte counts, timings, diagnostics, and live
progress in a new output directory. It does not retain BLOBs and therefore
does not perform decode verification. Use the full comparison harness above
when retained output and verification are needed. The 1,000-vector diagnostic
run has already passed that full verification cycle.

```sh
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
caffeinate -i bash "$HARNESS/run_gov2_diagnostic_stream.sh" \
  --bin ./bminv --input /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.unsorted \
  --out ./gov2-unsorted-horizontal-diag-d2-first10000-stream \
  --max-vectors 10000 --progress 100 \
  --level 6 --depth 2 --effort 8 --variants 64 --functions 7
```

The 10,000-vector prefix contains 19,275,322,824 raw bytes. The output
directory above must not exist. `caffeinate -i` keeps the Mac awake while the
command runs; omit it if not needed. The stream-only runner refuses a stale
binary and checks count, size-difference, and diagnostic accounting before
writing `report.txt` and `metrics.tsv`.

## Binary BWT wave experiment

Build `bminv` with the updated source, then use `run_gov2_bwt.sh` from
`utils/inv_list`. The script retains baseline and BWT collections, verifies
both against the raw input, times decode, and reports complete byte savings
and GAP complexity before and after BWT. BWT is applied to all 64 waves of
each sampled BIT block, with 80 packed primary-index bytes per accepted block.
The default `--sample 1` evaluates every eligible BIT block. Increase it for
a shorter exploratory run; the report always states the actual sample count.

```sh
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make
caffeinate -i bash "$HARNESS/run_gov2_bwt.sh" \
  --bin ./bminv --input /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.unsorted \
  --out ./gov2-unsorted-bwt-first10000-s1 \
  --max-vectors 10000 --progress 100 --level 6 --sample 1
```

The output directory must be new. The two retained BLOBs require several GiB
of free space. Experimental BWT records require the updated decoder.

## BWT followed by horizontal XOR

`run_gov2_bwt_xor.sh` compares ordinary encoding, XOR search on the original
waves, and BWT followed by XOR search on transformed waves. The combined
encoder accepts BWT+XOR only when its complete record, including 80 packed
BWT index bytes, beats the best original-wave XOR record. All three modes use
the same thread count and are independently verified and timed. Start with
1,000 Gov2 vectors; increase `--max-vectors` if the candidate funnel warrants
it. The output directory must be unused.

```sh
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make
caffeinate -i bash "$HARNESS/run_gov2_bwt_xor.sh" \
  --bin ./bminv --input /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.unsorted \
  --out ./gov2-unsorted-bwt-xor-first1000 \
  --max-vectors 1000 --progress 50 --threads 8 \
  --depth 1 --effort 8 --variants 64 --level 6
```

The report shows original-domain XOR wins, BWT-domain XOR candidates, complete
BWT+XOR wins over original XOR, actual bytes saved, GAP counts, and search
effort. The search report counts original and BWT-domain trials together.

To evaluate every eligible BIT block with the generic threaded engine, add
`--threads` to the same command. Both baseline and BWT compression then use
the same pool size; it selects 8 workers by default. Use
`--threads N` to override and `--batch-vectors N` to reduce the input window
from 512 vectors. Threaded mode requires `--sample 1`. The main thread retains
input and output order, and the report adds worker-seconds, ordered-wait time,
and peak raw input bytes for each compression mode. Build remains the separate `make`
step shown above.

```sh
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
caffeinate -i bash "$HARNESS/run_gov2_bwt.sh" \
  --bin ./bminv --input /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.unsorted \
  --out ./gov2-unsorted-bwt-mt-first10000 \
  --max-vectors 10000 --progress 100 --level 6 --sample 1 --threads
```

## Optional standalone correctness test (separate build and execution)

These tests do not alter or enable existing stress tests. They cover 32/64-bit
address modes, formula encodings and padding, truncated formulas, DAG reversal,
serialization levels 0–6, no-expansion fallback, disabled equivalence, GAP
exclusion, reference-XOR exclusion, set/count/range operations, concatenated
streams, mixed blocks, bookmarks and index scanning.

Build only the test executables (the app is still built with make):

```sh
REPO=/Users/anatoliykuznetsov/dev/BitMagic
TEST_BUILD=/Users/anatoliykuznetsov/Documents/Codex/horizontal-test-02
CXX=clang++ sh "$HARNESS/build_horizontal_tests.sh" "$REPO" "$TEST_BUILD"
```

Then execute only, using your already-built app:

```sh
sh "$HARNESS/run_horizontal_tests.sh" "$REPO" "$TEST_BUILD" "$REPO/utils/inv_list/bminv"
```

The optional test build defaults to AddressSanitizer; set `SANITIZE=none` if
unavailable. `EXTRA_CXXFLAGS` can carry architecture flags. An optional
`BASELINE_BIN=/path/to/pre-change/bminv` on the execution step adds a historical
byte-for-byte comparison of disabled output.
