# Native AArch64 NEON backend

`make BMOPTFLAGS=-DBMNEONOPT` selects `src/bmneon.h` on AArch64
(including Apple Silicon). ARM32 builds use scalar mode without `BMNEONOPT`.
The native path includes neither SSE headers nor the translator,
and retains the ordinary 64-bit scalar helpers (`BM64OPT`). Public APIs and
serialization formats are unchanged.

The implementation uses the baseline AArch64 Advanced SIMD intrinsics shared
by GCC and Clang; it does not require SVE, dot-product extensions, or an
Apple-specific instruction set. GCC runtime verification requires that toolchain; an Apple Clang build alone
does not validate it.

## Interface and contracts

The header supplies all 51 unique `VECT_*` hooks exposed by `bmavx2.h`, including
the 42 SSE2 hooks. It also provides native helpers for compressed GAP endpoint
summing and decoder OR/AND of unaligned input, which are dispatched directly
from `bmfunc.h` and `encoding.h`.

* Block Boolean operations preserve their existing return contracts: OR
  reports all-ones, AND/SUB/XOR return the OR reduction of the result, and
  digest operations report emptiness. AND-OR digest reports the emptiness of
  the source intersection, independently of the previous destination.
* Population count uses eight independent widening 16-bit accumulators.
  Each batch is limited to 2,048 words; the sum in any 16-bit lane is at most
  8,192. The final horizontal sum widens, so a full block correctly counts
  65,536 bits. Longer ranges use multiple bounded batches.
* Transition counting compares each bit with its predecessor and treats the
  initial bit as its own predecessor. Uniform zero/one blocks have one run.
  The XOR path uses the initial bit of the XOR result.
* In BM's naming, `SHIFT_R1` moves bits toward larger bit indices (`word << 1`).
  The masked shift returns the carry before masking and updates the digest.
* GAP counting sums interval lengths, rather than applying population count
  to compressed endpoints. GAP conversion retains the caller's existing
  sufficient-capacity precondition.
* Ordinary vector stores implement the stream hooks. The interface does not
  require a non-temporal store policy. Word-pointer kernels need word alignment;
  explicit unaligned copy/decoder paths use byte loads.

## Correctness

Both stress drivers include `tests/neon_test.h`. Native builds run its tests
with `-neon`, `-ll`, or the default all-tests selection. Non-native builds
reject `-neon` with a diagnostic. The historical SSE-specific `TestSIMDUtils()`
remains unchanged; the NEON suite is independently compiled and invoked.

After building fresh binaries with the flag above:

```sh
(cd tests/stress && ./bmtest -neon -complexity -silent)
(cd tests/stress64 && ./bmtest64 -neon -complexity -silent)
```

The shared suite checks all 51 hooks against scalar expectations, including
adversarial patterns, fixed-seed random blocks, carries, flags, selected digest
waves, short tails, unsigned comparisons across the signed boundary, duplicate
scatter destinations, and exact-sized search buffers. It also checks multi-batch
count overflow, unaligned decoder tails, public operations, mixed dense/GAP
storage, and serialization/serialized AND/OR at compression levels 0–6.

Run a sanitizer build of these tests as well. Existing broad low-level,
Boolean, aggregation, shift and serialization stress groups remain useful
for longer regressions. `-neon` is a bounded kernel/integration test, not a
replacement for the entire stress suite. In particular, the historical
`TestFindBlockDiff()` contains billions of nested comparisons.

## Benchmarks

The existing `tests/perf/perf.cpp` driver accepts these optional groups:

```sh
./bmperf -bench-kernels -bench-dense -bench-gap \
    -bench-blocks 64 -bench-repeats 1000 -bench-seed 1337
```

* `-bench-kernels`: population count, Boolean counts, runs, fused runs/count,
  and fused XOR runs/count on raw dense blocks.
* `-bench-dense`: public count/distance, interval counting and XOR materialization.
  Fixtures are checked to contain the requested number of bit blocks and zero
  GAP blocks. Materialization includes result allocation/reuse and counting.
* `-bench-gap`: endpoint summation and bit-to-GAP conversion on regular runs.

Groups can be combined. Without an explicit group, the ordinary benchmark
sequence runs along with the bounded SIMD checks. `-silent` runs benchmark
checks once without timing output.
Defaults for explicit workloads are 64 blocks, 200 repeats, and seed 1337.
Valid block counts are 1–4096. The default full perf run executes one grouped
measurement with 64 blocks and 20,000 sweeps per workload.
Explicit benchmark selections report one line per workload in the existing
`test name; duration` format. Five samples are measured and the median elapsed
time across all selected sweeps is reported. Generation, import and fixture
verification are outside each workload's timing. A compiler memory barrier
prevents invariant sweeps being hoisted. The ordinary no-argument perf run
also executes the bounded kernel, dense and GAP checks and reports them together
as one `SIMD fast kernels; duration` line. A scalar build runs these checks via
the scalar fallback kernels.

Compare one-block hot inputs, 64-block inputs, and 4,096-block inputs (two raw
inputs total 64 MiB). Run native and translated binaries sequentially with the
same compiler, flags, seed, block count and repeat count. Compare checksums
before interpreting speedups. Repeat in reverse order to expose warmup or
frequency effects; the output includes compiler and backend identity.

For a translated baseline, compile this updated performance driver against a
separate copy of the pre-change `src` headers. Do not toggle flags on an existing
object directory and assume make will rebuild objects when flags change.
The existing `tests/perf/build_all_arm.sh` and `tests/stress/build_all_arm.sh`
continue to produce the same release executable names.
