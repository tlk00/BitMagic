# Test layout

The root contains build files, scripts and documentation. Keep C/C++ sources in
the appropriate subdirectory:

- `common/`: shared internal test headers, fixtures, allocators and generators.
- `regression/`: standalone focused regression programs, each with its own main.
- `stress/` and `stress64/`: 32-/64-bit addressing stress drivers and local helpers.
- `perf/` and `perf64/`: performance drivers and local benchmark helpers.
- `test_threads/` and `bench-tpch/`: specialized test programs.

Common headers may depend on other common headers and library headers. Keep
helpers used by just one driver with that driver. Some existing fixture headers
expect types declared by their including driver; they are not standalone APIs.

## Stress workload modes

Both `bmtest` and `bmtest64` accept `-fast` (or `--fast`) in debug and optimized
builds. It keeps the selected sections, assertions, fixture sizes, and large
indices, while reducing audited stress repetitions and sampling SELECT sweeps.
Without section selectors, every normally enabled section still runs:

```sh
bmtest -fast
bmtest64 --fast --bvops
```

Both `stress/run_all_pl.sh` and `stress64/run_all_pl.sh` accept `-fast` or
`--fast` and forward it to every worker, for example
`./run_all_pl.sh --build stress_debug_neon --fast -j 4` from `stress/`.
The runners reject combining fast and pedantic modes before launching workers.
Fast runs may read existing scheduling weights, but never create or update the
weight history. Per-run elapsed times are still displayed.

All stress-binary options and existing aliases accept either one or two leading hyphens.
`-fast` and `-pedantic` cannot be combined, in either order. Omitting both retains
the normal workload; pedantic workloads are unchanged. Random generators retain
their existing seeding and distributions; fewer trials consume fewer random
numbers, so subsequent cases need not reproduce a normal run.

Fast budgets cap 32-bit operation/aggregator trials at eight. The 64-bit
aggregator uses two complete random fixtures per operation, retaining every
argument-ladder query and all fixed full-universe unit cases. The 64-bit
operation test uses six trials, retaining both allocation strategies and the
four-trial clear cycle followed by reuse. Its generated-vector checks still
cycle through all compression levels. The 64-bit operation test samples one
wide-vector serialization pair per trial, rotating through both sources and both wide patterns across four trials;
normal/pedantic runs retain all eight calls per trial, including duplicates.
Shift/AND aggregation uses two collections and eight repeated checks per collection in the 32-bit driver; the 64-bit driver keeps one complete 20-vector
collection and two identical queries, preserving reset/reuse and both count modes.
Deserialization keeps six trials to exercise its five-trial clear cycle and subsequent reuse.
Compressed collection tests keep eight trials. The 32-bit driver also reduces
Kleene logic trials to eight, randomized block splitting to 256 trials, randomized
delta tests to 16 passes, and interpolation coding to 256 trials per pattern.
Interpolation spacing factors span the original range, retaining the special
first two cases and the original buffer capacities.

SELECT retains its unit cases and full block/vector fixtures. Fast query sweeps
keep positions near both ends and around each 4096-position boundary, including
the 64-bit driver's rank/select queries. Normal and pedantic SELECT sweeps remain
exhaustive.

Interval-enumerator stress keeps the full fixture and both unoptimized/optimized
passes. Fast mode samples spacings around word/block boundaries and the spacing
limit; the 32-bit driver's range-copy checks are capped at 16 per pass.

Additional audited fast workloads:

| Workload | Fast policy |
| --- | --- |
| Random word comparison / block-to-GAP conversion / GAP range checks (32-bit driver) | 4096 / 256 / 64 trials; full blocks and validation retained |
| GAP stress densities / array encoder inputs (32-bit driver) | Sample the original domain; retain full validation and original buffer capacities |
| Bulk insertion spacings | Sample boundary spacings; retain the original range, including the 64-bit driver's indices above 2^32, and the fully populated spacing-1 case |
| Repeated aggregation on a shift collection | 8 checks of the original collection |
| Large-vector insertion / erasure | 8 operations in the 32-bit driver; 2 in the 64-bit driver, retaining both inserted bit values and post-optimization erasure; original generated vectors retained |
| Repeated large-vector shifts (32-bit driver) | 65 shifts, crossing both word widths; the 64-bit driver's existing 33-shift case stays unchanged |
| Block-change patterns, 32-/64-bit bitscan patterns, long log2 sweep, AVX2 unsigned comparison (32-bit driver) | Sample the original input domain, retaining endpoints and power-of-two neighbors |
| Sparse-vector increments | Sample counter values across the full original range; restore skipped pre-increment values at every original index, including indices above 2^32 in the 64-bit driver |
| Sparse EQ scans / filter factors | Sample input values; keep each complete fixture and its checks |
| Compressed gather | One sorted and one shuffled/unsorted gather per sampled NULL spacing; keep the original fixture size, maximum spacing, and full query array |
| Repeated identical string serialization (32-bit driver) | 2 passes, retaining the comparison with the previous encoded buffer |
| Gamma encoder buffers (32-bit driver) | 66 complete buffers, covering all 65535 legal values and wraparound |
| Random bit swaps (32-bit driver) | 64 swaps per original full fixture |
| Sparse interval-shift/join stress | Sample original offsets through the upper endpoint; keep full fixtures, the detailed first pass and an ordinary pass |
| Simple random fill (32-bit driver) | Keep every fill/clear operation and original index range; sample expensive index/range checks, always including the last mutation |
| Range copy | Sample repeated ascending/descending and moving-endpoint queries; preserve source vectors, both optimization passes, fixed cases, sparse-fixture transitions, and high-index word boundaries |
| Random subsets of the full inverted universe | The 32-bit driver keeps empty, small, and full-copy trials and omits huge fractional trials in fast mode. The 64-bit driver keeps empty, singleton, and 256-bit requests, reducing repeated full-universe index builds. Populated fixtures with at most one million set bits retain all fractional/complement checks; larger populations keep 0/1/2/15/full-copy requests. The 32-bit full-universe fixture keeps all original small requests. AND/OR target cross-checks omit repeated half-population subset requests above one million bits, retaining both full-source targets; smaller populations keep all four target cases. |

The shared input sampler preserves dense low/end values and power-of-two
neighbors in addition to samples across the original domain. Outside fast mode
it advances by one, preserving the original sweep and random-number consumption.

Release fast-mode tuning targets roughly 10–15 minutes per harness group on the
reference machine. This is a tuning target, not a deadline: no timer stops a
section and no harness timeout discards unfinished checks. Fixed fast budgets
apply in debug builds too; machine speed and parallel contention affect runtime.
Float scanner distributions run 8 complete reference queries on the original
fixtures; all original random query bounds are still generated. Rank splitting,
GT/EQ scanner queries, basic rank validation, string comparisons and
deserialization ranges, sparse-vector range partitions, sparse gather-deserialization
ranges (including both bookmark passes of the large 32-bit RSC fixture),
and 32-bit merge split positions sample their original query domains.
Full block-number rank-index checks in the 64-bit driver also sample queries.
The 32-bit singleton shift sweeps restore the equivalent pre-shift singleton at
sampled positions and check each resulting shift, including high positions and
right-shift overflow. The 64-bit singleton left-shift sweep similarly restores
sampled positions throughout its original 100,001-position high-index window.
Its 33 shifts of the complete large vector remain unchanged.
Normal and pedantic traversal remains exhaustive.

Import stress samples original input sizes, including the original maximum;
sparse gather queries and signed EQ values sample their full original domains.
Randomized bvector range-deserialization sweeps keep 32 trials per fixture and
bookmark pass, entering the original random gallop after 8 trials. Compressed
decode validation retains all sweep shapes with 4 full-prefix and up to 8 other
queries per sweep. Full compressed-value comparison uses an indexed copy for
large unsynchronized fixtures in fast mode; small fixtures retain unsynchronized
lookup validation. Large unsigned compressed assignment builds the full expected
fixture by append, then samples synchronized/unsynchronized set/increment pairs.
String verification keeps every value and sequential iterator check while
sampling independent iterator construction and full-vector scanner queries; common-prefix and sorted-search queries are sampled.

Large string/numeric midpoint-erasure helpers perform at most 64 erases in fast
mode; small collections still drain to empty. Large non-nullable float streaming
checks retain complete sign/exponent/mantissa plane equality and sample redundant
bit probes; small gather masks and assigned NULL maps remain fully enumerated.
String indexed-gather tests retain every source, NULL policy and serialization
variant, with first/last present and missing point queries from the original lists,
2 block-island trials, and 3 random-island trials covering 1, 8 and 64 islands.

Dense bvector validation above one million set bits retains complete comparison
against the independent flat reference bitmap. It imports that reference once,
compares every bit, independently checks the last bit from raw words, and samples
redundant iterator/range API probes. Smaller fixtures keep the original validator.
Normal/pedantic validators are unchanged. The reference-only `bvector_mini::size()`
accessor exposes its bit capacity for this full-content check.

Compressed-gather NULL spacings use the original interval sampler, including the
original maximum. For early fixtures, exact assigned-index/value enumeration
checks the complete NULL layout without scanning every NULL position. Sparse
interval-shift and join offsets sample the original range through 9,999,999;
selected fixtures retain their full spans and complete reference comparisons.

The large integer serialization stress fixture keeps every legacy XOR on/off and
optimization pass. Fast mode runs its identical stream companion once per fixture
rather than four times; normal and pedantic companion calls remain unchanged.

The 64-bit filter unit fixture retains its high indices and complete value/NULL
plane comparison against a fixed expected vector; fast mode samples redundant
point queries across the full logical span. Normal/pedantic point sweeps stay
exhaustive.

Rank-compression stress retains both original 64-bit sizes (over 8 and 17 billion
positions), full ordinary compression/decompression and both optimization passes.
Fast mode uses the existing 100-bit random-subset queries for the per-source
compressor on each complete index; full dense per-source queries remain in normal
and pedantic runs. Basic rank-compressor unit tests remain unchanged.
Other passes that grow fixtures or sweep positions remain unchanged. Full debug
runtime therefore still depends on these large workloads; fast mode has no
wall-clock cutoff.

Fast mode also samples independent inner position checks in the 32-bit block
first-difference and insertion/erasure sweeps; cumulative outer mutations remain
intact. Range-visitor steps cover at most roughly 16 full-range query spacings.
Known remaining costs include full matrix/string fixture construction, fixed
full-universe aggregator unit cases, and serialization cases whose pass number
changes the populated fixture.

The 64-bit constant sparse-vector range-partition stress fixture also has an
existing bound mismatch: `sv_from` is near the maximum index while the query-loop
upper bound is 120,000, so those query loops do not run. This pass preserves that
normal/pedantic behavior; the separate linear fixture and high-index unit cases run.

The 64-bit `RangeCopyTest()` has a separate existing coverage gap: its inverted
T4/T5 start expressions subtract `to_max == bm::id_max - 1` and then another
offset, underflowing the unsigned index and skipping those sweeps. Fast-mode
sampling does not change that behavior; the high-index T1/T2/T3 cases still run.

## CMake

From the repository root, configure an out-of-tree build and select the targets
needed for the check. Native Mac NEON builds can use:

```sh
cmake -S . -B /tmp/bm-tests -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS=-DBMNEONOPT
cmake --build /tmp/bm-tests --target bmtest bmtest64 bm_debug_allocator_test
```

Executables are written to the build directory's `bin/`. Existing stress and
performance target names are preserved. Shared headers appear in the generated
IDE projects under `Common Test Headers`.

The standalone targets are `bm_debug_allocator_test`,
`bm_serial_compression_levels`, `bm_serial_gap_normalization`,
`bm_serial_opcode_registry`, and `bm_serial_released_compat`.
The first four are registered with CTest; build them before running:

```sh
ctest --test-dir /tmp/bm-tests --output-on-failure
```

`bm_serial_released_compat` requires a directory of fixtures generated against a
released library, supplied as its command-line argument. Building that target
does not verify released-format compatibility. The allocator test also accepts
`block-size` and `ptr-size` to deliberately abort and check failure diagnostics;
these negative checks are not part of its ordinary CTest invocation.

## Make

Existing Make entry points remain available. After moving headers in an existing
checkout, regenerate local dependency files with `make dep` in each affected
driver directory (`stress`, `stress64`, and `perf64`) before rebuilding. Their
Makefiles also track the common headers explicitly.
