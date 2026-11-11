# Cross-platform inv_list reference tests

Run from `utils/inv_list`, using the binary built for the test platform:

```bash
bash run_inv_test_pl.sh --bin ./bminv_neon \
  --reference /Volumes/USB_1TB/gov2-reference-Nov-2026 \
  -j 8 --out /Volumes/USB_1TB/inv-test-full-$(date +%Y%m%d-%H%M%S)
```

| Mode | Work for every selected sorted/unsorted level |
|---|---|
| Normal | Read golden vectors once for in-memory serialization/deserialization at every selected level, with exact equality and per-level size/timing reports. |
| `--pedantic` | Normal round trips plus saved-level threaded/serial raw-ID verification, independent decode, and full reference integrity audits. |
| `--fast` | Validate golden against all raw IDs once per dataset, then compare selected levels against golden with exact bvector equality. |

Use `--start-level 6` to test levels 6–10 inclusive for both datasets (10 cases).
The default is level 1, covering all 20 cases. For a reduced-I/O fast run:

```bash
bash run_inv_test_pl.sh --bin ./bminv \
  --reference /Volumes/USB_1TB/gov2-reference-Nov-2026 \
  --fast --start-level 6 -j 10 --read-memory-mib 8192 \
  --out "/Volumes/USB_1TB/inv-test-level6-$(date +%Y%m%d-%H%M%S)"
```

`--golden-level 6` is the default. Choose any level from 1 through 10.
For each dataset, the runner first validates that golden stream against every raw
ID using the current executable. Only after it passes do other selected levels
use it as their expected collection. Golden validation runs even when its level
is outside `--start-level`; it is a prerequisite and does not increase the selected
case count. If the golden is selected, its completed raw check supplies its
threaded verification in fast mode. Fast/pedantic historical checks decode the
saved selected streams. Normal round trips use the golden stream once across
the selected encoder levels and compare each fresh decode exactly. A golden failure
stops the run before dependent comparisons for that dataset.

Fast mode omits checksum/history audits, rebuilding and decode-only passes. Summaries report the actual case count; `test.json`
records the level range and count. Normal and pedantic modes also accept
`--start-level`. Pedantic mode performs the final all-file checksum scan.

The raw inputs must be `<reference>/Gov2Flat/gov2.sorted` and `gov2.unsorted`,
Pedantic mode checks their recorded input hashes; uncached runs validate the golden
against the supplied raw inputs. The current bminv verifier imports each raw ID
list into a bvector and calls `equal()` on the decoded vector. Other streams are compared against that validated logical baseline.
Serialized size/checksum changes are reported, never treated as logical failures.
The target executable must support `-golden-bvin` and `-golden-index-in`.

`--bin` accepts an executable name on PATH or a path. `-j N` / `--threads N`
sets the worker budget; default is logical CPUs minus two (at least one).
Cases run sequentially. Build binaries separately; the script snapshots its
chosen executable before running so concurrent rebuilds cannot alter later cases.

`--out` must be fresh and disjoint from the reference directory, including
resolved symlink aliases. Commands only pass reference paths as read inputs;
logs, summaries and process working directories are under `--out`.
The script never chmods, removes, or overwrites reference inputs. Pedantic mode
checks saved checksums before testing and again on success. Normal and fast
modes use the saved TSV metadata and basic matrix/file discovery checks.
This is protection against script mistakes, not sandboxing of an untrusted binary.

The first failure stops the run. A target error retains its positive exit code;
a signal maps to `128 + signal` (SIGSEGV → 139). The console names the case,
phase and log. Core dumps are disabled. `FAILURE.txt`, partial `summary.tsv`,
and per-phase logs/commands/results remain for diagnosis. Success creates
`PASS.txt`. `--reserve-gib 8` sets the output free-space reserve.

Fixture checks (no production reference writes):

```bash
python3 reference/test_inv_reference.py --bin ./bminv
```

The optional `--crash-test` fixture deliberately sends SIGSEGV to a small Python
wrapper to check exit 139 propagation. macOS may display a crash dialog; routine
fixture tests omit this test. The production runner never injects crashes.

Add `--silent` to suppress console progress and per-phase announcements while
retaining errors, completion summaries, and full progress in the phase logs.
Direct `bminv --silent` suppresses its own progress even with
`-progress-vectors N`. Legacy `-s` keeps the earlier behavior: requested
structured progress remains available to existing harnesses.

## Indexed verification I/O

Threaded verification and decode use two sequential adjacent-BLOB batch slots.
The coordinator reads the next batch while workers decode/compare the active
batch, then drains active tasks before exchanging the slots. There is one reader
and the configured worker pool; each file advances sequentially. Normally each
slot can use half the configured budget. The combined staged input bytes across
both slots stay within the budget. A single vector/pair above half the budget
can use up to the full budget after active work drains. A pair above the full
budget is rejected before payload allocation. Buffered reads perform read-ahead
in the coordinator; mmap windows request advisory prefetch and can still fault
in workers. Output reports `read_ahead_batches` and combined peak staged bytes.
Use `--read-memory-mib 8192` for an 8 GiB budget, or `--read-memory-mib 16384`
for 16 GiB of combined staged input data. Raw verification counts raw IDs plus
encoded bytes; collection comparison counts encoded bytes from both files. Positive MiB values are
limited by platform address/stream sizes. The default remains 512 MiB.
The setting is an upper bound; batches allocate only the data they contain. The index, worker
bitvectors/decoder/comparison workspaces and filesystem cache need additional RAM.
The existing 512-vector batch limit can result in smaller batches.

`--read-mmap` selects Unix read-only compressed-file windows and advisory
sequential/read-ahead hints. Raw inputs remain sequential buffered reads.
Windows uses buffered batches. Both backends preserve exact logical comparison and
exact indexed BLOB consumption. These options apply to fast and normal threaded
read phases. Commands and test.json record the selection.

Focused backend/error checks from utils/inv_list:

```bash
python3 experiments/current/test_buffered_read.py --bin ./bminv
```

## Worker failure propagation

Read/compression workers capture exceptions with `std::exception_ptr`. After an
acquire completion wait, the main thread checks both the captured exception and
the pool task error code. Read failures name the verification/decode phase,
source files, zero-based vector index and indexed BLOB offsets. Workers finish
before their batch buffers or mappings are released. Failures exit with status 1;
the reference runner retains the diagnostic in the phase log, marks the summary
failed and writes FAILURE.txt without a success marker.

Induced failures in `test_buffered_read.py` cover first/later ID mismatches and
truncated BLOBs with buffered/mapped reads and quiet/progress/silent execution.
`test_inv_reference.py` exercises real worker failures through the runner's logs,
result JSON, summary and failure markers. A standalone background-pool check also
covers allocation/non-standard exceptions and pool-level error codes:

```bash
c++ -std=c++17 -O1 -pthread -I../../src -I. \
  experiments/current/test_worker_error.cpp -o /tmp/inv-worker-error
/tmp/inv-worker-error
```

Direct indexed collection comparison, after independently validating the golden:

```bash
./bminv -verify -threads 8 -read-memory-mib 8192 \
  -bvin level10.bv -bv-index-in level10.idx \
  -golden-bvin level6.bv -golden-index-in level6.idx
```

Each file has its own index and sequential batch window. The combined staged
bytes obey the read budget; indices, per-worker decoded vectors, codec scratch
space and mapping page alignment need additional memory. Failures include both
file names and indexed BLOB ranges. The focused fixtures induce golden/raw
mismatches, first/later collection mismatches, count differences, malformed BLOBs
and extra bytes, and check runner failure logs and markers.

## Cached SHA-256 validation

Add `--hash` explicitly to use the integrated checksum path. For example, after
rebuilding the platform binary, run from `utils/inv_list`:

```bash
bash run_inv_test_pl.sh --bin ./bminv \
  --reference /Volumes/USB_1TB/gov2-reference-Nov-2026 \
  --fast --hash --start-level 6 --golden-level 6 \
  -j 10 --read-memory-mib 16384 --read-mmap \
  --out "/Volumes/USB_1TB/inv-test-sha-$(date +%Y%m%d-%H%M%S)"
```

The default cache directory is `<reference>/validation-cache`. Use
`--hash-cache-dir PATH` to choose a separate persistent directory. On first use
for each dataset, the runner completes full golden/raw-ID equality verification,
then asks inv_list to generate the level-6 cache. Subsequent runs reuse that
cache and verify each selected candidate against it. `--golden-level` chooses
another baseline if desired. Plain `--fast` continues using exact golden
comparison. Normal mode also accepts `--hash` to anchor its golden stream to the
cache before in-memory round trips. Pedantic mode adds serial raw checks and
integrity audits.

The reference streams, indices and raw files are read inputs. Hash mode creates
its authorized cache subdirectory, or writes into the explicitly selected cache
directory. Logs and summaries remain under fresh `--out`. Cache keys
bind the canonical format, dataset, level, vector/ID counts and frozen TSV source
SHA-256 values. The binary cache format is portable between target platforms.
A JSON sidecar records provenance, origin file metadata and the cache checksum.
Initial generation checks source metadata for concurrent mutation and verifies
the generated header's source SHA-256 against the frozen golden stream hash.
Ready caches are immutable. Unique staging files and a per-cache mkdir lock
protect publication. A failed build is cleaned without publishing a ready pair.
A lock left by an interrupted process is reported with its path for inspection.

On reuse, the runner reads the small cache and its sidecar, and the C++ verifier
also validates the cache footer. Fast hash runs trust the frozen source hashes
in `references.tsv`; pedantic mode performs the original reference integrity audits.
A source hash/scope change selects a different cache key. Golden and raw input
contents are read during initial generation only, apart from a selected golden
level subsequently being tested as a candidate itself.

Each record contains 32 SHA-256 bytes and an 8-byte little-endian population
count. Verification first compares the decoded bvector's population count, then
computes its canonical block hash. An equal-cardinality ID substitution must
still pass the SHA-256 comparison. The cache has a versioned 64-byte header and
a 32-byte SHA-256 footer. Empty and allocated-zero blocks are omitted, GAP blocks
expand into aligned scratch, FULL blocks supply all-one bytes, and each occupied
8 KiB logical block is followed by its uint64 LE absolute block number.

`src/bmbvhash.h` provides the STL-free `bm::bvector_hash<BV, Hash>` class. It
borrows a functor supporting `reset`, `update`, `final` and `digest_size`, and
resets it for each `compute`. It retains two aligned scratch blocks through the
supplied `BV::allocator_type`; callers can pass `bv.get_allocator()`. Its format
and traversal are independent of the application crypto provider and disk cache.
One traversal/functor instance is reused by each worker across its batch jobs.
The application pipeline keeps two slots under the shared staged-input budget;
indices, decoded vectors and codec/hash workspaces require additional memory.

Direct commands support independent raw-derived cache experiments as well:

```bash
./bminv -threads 8 -read-memory-mib 16384 \
  -bvin level6.bv -bv-index-in level6.idx -hash-cache-out level6.sha256
./bminv -threads 8 -read-memory-mib 16384 \
  -u32in Gov2Flat/gov2.sorted -hash-cache-out raw.sha256
./bminv -verify -threads 8 -read-memory-mib 16384 \
  -bvin level10.bv -bv-index-in level10.idx -hash-cache-in level6.sha256
```

Direct cache creation requires a fresh output path. The runner supplies the
initial golden/raw proof and provenance sidecar automatically. Direct raw cache
creation imports sorted, unique uint32 IDs into temporary bvectors with a bounded
4,096-ID conversion chunk, and produces the same logical records. Raw payloads
are staged within the configured budget before import. A single record exceeding
the full budget fails before payload allocation.

### Crypto providers and dependency-free builds

macOS uses CommonCrypto, Windows uses CNG (`bcrypt`), and Linux uses OpenSSL EVP.
The generic traversal header has no crypto or STL dependency. Linux Make builds
find libcrypto with pkg-config, and root CMake builds use `OpenSSL::Crypto`.
When OpenSSL is unavailable these builds select exact-only operation. Explicit
fallback builds can use:

```bash
make BM_NO_SHA256=YES
# Or configure the root CMake project with -DBM_INV_SHA256=OFF.
```

A direct compiler build can define `BM_NO_SHA256`. Hash commands then report that
SHA-256 support is disabled, while exact verification is available. Windows
MinGW links `bcrypt`; MSVC also receives the library directive in the adapter.

A possible embedded fallback is a small reviewed SHA-256 provider implementing the
same functor interface. One candidate is `amosnier/sha-2`, a portable C99
implementation with a streaming API and a 0BSD licensing option. Pin and review
its source before inclusion. [Upstream implementation](https://github.com/amosnier/sha-2),
[upstream license](https://github.com/amosnier/sha-2/blob/master/LICENSE.md).
 Before vendoring one, check its license/provenance and
run standard known-answer vectors, incremental/chunk-boundary tests and the
independent Python canonical oracle on supported compilers. The current
fallback is an exact-only build and uses the installed platform provider when
hashing is enabled.

Focused checks include cardinality-only corruption with a recomputed cache
footer, equal-cardinality substitutions, raw/BV digest equivalence, cached-only
reads, cache corruption/truncation, and failure cleanup:

```bash
python3 reference/test_hash_cache.py --bin ./bminv --inv-bin ./bminv
python3 reference/test_inv_reference.py --bin ./bminv
c++ -std=c++17 -O1 -I../../src experiments/current/test_bvector_hash.cpp \
  -o /tmp/test-bvector-hash
/tmp/test-bvector-hash
```

The generic traversal fixture also supports `-DBM64ADDR` and checks allocator
propagation, two-block workspace release, reuse/reset, GAP/bit/FULL equivalence,
holes and empty vectors. The initial block-SHA experiment's temporary datasets,
caches and builds were removed after its results were recorded in
`experiments/block_sha/RESULTS.md`.

## In-memory normal validation and pedantic coverage

Normal mode now performs serialization verification in RAM. For each dataset it
reads the golden collection once in the round-trip phase. Each worker decodes
one golden vector, optionally checks cached POPCOUNT/SHA-256, then serializes and
freshly deserializes that vector at every selected level. Exact bvector equality
checks the rebuilt contents. Decoders require exact BLOB consumption. Optimizing
the golden representation is checked against an untouched copy before encoding.
The temporary optimization copy is released before rebuilt vectors are created.
No candidate BV/index files are written. Per-level vector counts, encoded byte
totals and encode/decode/comparison worker times appear in `roundtrip-levels.log`;
`summary.tsv` keeps a row for every dataset/level and links to that shared log.

With `--hash`, the first run establishes the golden cache through the complete
raw-ID check. Warm-cache normal runs read the golden stream and its index once,
plus the small checksum cache. Without `--hash`, the runner first checks the
golden against raw IDs on every run, then performs the same round-trip pass.

`--pedantic` adds the saved historical stream checks at every selected level.
Those include threaded logical verification, serial raw-ID verification and an
independent decode, plus the reference checksums and final immutability audit.
This exercises historical encoding choices that a new encoder might stop
emitting. `--fast` and `--pedantic` are mutually exclusive. Plain `--fast` and
`--fast --hash` continue to scan selected saved streams without reserialization.

Normal cached validation from `utils/inv_list`:

```bash
bash run_inv_test_pl.sh --bin ./bminv \
  --reference /Volumes/USB_1TB/gov2-reference-Nov-2026 \
  --hash --start-level 6 --golden-level 6 -j 10 --read-memory-mib 16384 \
  --out "/Volumes/USB_1TB/inv-normal-$(date +%Y%m%d-%H%M%S)"
```

Pedantic validation of the same levels:

```bash
bash run_inv_test_pl.sh --bin ./bminv \
  --reference /Volumes/USB_1TB/gov2-reference-Nov-2026 \
  --pedantic --hash --start-level 6 --golden-level 6 -j 10 --read-memory-mib 16384 \
  --out "/Volumes/USB_1TB/inv-pedantic-$(date +%Y%m%d-%H%M%S)"
```

The normal round-trip phase divides `--read-memory-mib` equally between encoded
read-ahead/record staging and rebuilt BLOB buffers across workers. Each worker's
buffer is limited to the second half divided by the configured worker count.
A conservative `calc_stat().max_serialize_mem` bound is checked before allocation;
exceeding it reports the vector and asks for more memory or fewer threads.
Old storage is released before growing a buffer. Decoded golden/restored vectors,
serializer/decoder scratch and indices require additional RAM. Buffers grow only
as needed and are reused across levels and vectors. At most two decoded vectors
are retained by a worker at once. Buffered reads are the baseline; `--read-mmap`
remains an optional Unix backend.

Direct in-memory round-trip command:

```bash
./bminv -verify -roundtrip -threads 10 -read-memory-mib 16384 \
  -bvin level6.bv -bv-index-in level6.idx \
  -roundtrip-start-level 6 -roundtrip-end-level 10 -hash-cache-in level6.sha256
```

Omit `-hash-cache-in` for independently anchored exact-only runs. The dependency-free
build supports this path too. To induce valid-but-wrong rebuilt contents and
extra indexed bytes in real encoder subclasses, run the focused worker fixture:

```bash
c++ -std=c++17 -O1 -pthread -I../../src -I. \
  experiments/current/test_roundtrip.cpp -o /tmp/test-inv-roundtrip
/tmp/test-inv-roundtrip
```
