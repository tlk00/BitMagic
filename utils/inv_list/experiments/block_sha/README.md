# Canonical-block SHA-256 experiment

This standalone development tool compares the existing exact bvector equality
check with a representation-independent SHA-256 check. It can save vector
fingerprints from a compressed golden collection or from raw integer lists,
then verify other compressed collections using only their BLOBs and the cache.
It uses application-internal block APIs and OpenSSL EVP. Production inv_list
and its reference runner are independent of this experiment.

Build from the repository root with an installed OpenSSL development package.
On ARM64, enable native BitMagic NEON decoding as shown below.

```bash
cmake -S utils/inv_list/experiments/block_sha -B /tmp/block-sha-build \
  -DCMAKE_BUILD_TYPE=Release -DBLOCK_SHA_NEON=ON
cmake --build /tmp/block-sha-build --config Release -j 4
```

If CMake needs the OpenSSL location on Homebrew macOS, supply
`-DOPENSSL_ROOT_DIR=/opt/homebrew/opt/openssl@3`. Use the default
`BLOCK_SHA_NEON=OFF` for other architectures. Windows builds need an OpenSSL
installation discoverable by CMake. With a Visual Studio generator the executable
is under the build directory's `Release` subdirectory. Windows hardware validation
is pending.

Run all experiment phases against an existing reference directory. The output
must be fresh and outside the reference tree. Default scope is the whole dataset.

```bash
python3 utils/inv_list/experiments/block_sha/run_experiment.py \
  --bin /tmp/block-sha-build/block_sha \
  --reference /Volumes/USB_1TB/gov2-reference-Nov-2026 \
  --dataset sorted --threads 8 --memory-mib 16384 --batch-vectors 2048 \
  --out /tmp/block-sha-sorted-run
```

The runner executes exact level-6/level-10 comparison, SHA comparison of the same
pair, golden-BV cache generation, raw-ID cache generation, byte comparison of
both caches' logical records, and cached validation of levels 7, 8, 9 and 10.
`--dataset unsorted` chooses the other raw collection and its associated streams.
`--max-vectors N` explicitly restricts every phase to the first N vectors.
`--mmap` requests Unix read-only mapping windows with advisory read-ahead hints.
Each phase saves its exact command and combined log. `results.json` contains
measurements, and completion creates `PASS.txt`. Failed validation creates
`FAILURE.txt` and exits with failure. Cache output paths must be fresh.

A subsequent candidate scan can reuse the raw-derived sorted cache from the
recorded experiment directly. This command opens the candidate and cache only.

```bash
/tmp/block-sha-build/block_sha --mode cached \
  --candidate /Volumes/USB_1TB/gov2-reference-Nov-2026/gov2.sorted.level-07.bv \
  --candidate-index /Volumes/USB_1TB/gov2-reference-Nov-2026/gov2.sorted.level-07.idx \
  --cache /tmp/inv-block-sha-gov2-sorted-20261111/raw.sha \
  --threads 8 --memory-mib 16384 --batch-vectors 2048
```

Individual modes use these arguments.

| Mode | Required inputs | Work |
|---|---|---|
| `exact` | `--candidate BV --candidate-index IDX --golden BV --golden-index IDX` | Decode both vectors and call `equal()` |
| `pair-sha` | Same pair inputs | Decode both, hash canonical blocks and compare digests |
| `cache-bv` | `--candidate BV --candidate-index IDX --cache NEW` | Decode and persist golden digests |
| `cache-raw` | `--raw RAW --cache NEW` | Import raw IDs into temporary bvectors and persist the same digests |
| `cached` | `--candidate BV --candidate-index IDX --cache EXISTING` | Decode and compare against persisted digest/cardinality records |

All modes accept `--threads N`, `--memory-mib N`, `--batch-vectors N` and
`--max-vectors N`. Pair inputs need matching vector counts even for prefix runs.
Cached runs require the cache count to match the selected candidate scope.
The raw format is little-endian uint32 count followed by that many uint32 IDs.
IDs must be strictly increasing, unique and below BitMagic's reserved maximum ID.
Empty raw vectors are supported. The raw builder imports through a reusable
4,096-ID chunk inside each task, with a fresh temporary bvector for each vector.

## Canonical stream and disk format

Each vector digest starts with the 15 ASCII bytes `BMBLOCKSHA256V1`.
Traverse logical blocks in ascending absolute block order, skipping zero blocks
including allocated bit-blocks whose contents are all zero. A GAP block is
expanded with BitMagic's conversion into allocator-aligned heap scratch. FULL
blocks supply the all-ones block. Every occupied block contributes its 8,192
logical bytes, encoded as little-endian 32-bit words, followed by its absolute
block number as little-endian uint64. On little-endian machines the bit-block
bytes feed EVP directly. No set-bit enumeration is needed. Trailing logical
size and allocation shape do not contribute to this set-content digest.

One streaming SHA-256 context hashes the whole vector. Source compression level,
serialized byte identity, GAP optimization and block allocation choices therefore
produce identical vector digests when their sets are equal.

The cache has a 64-byte header, 40-byte records in vector order and a 32-byte
SHA-256 footer over the preceding header and records. All integer fields are LE.

| Header offset | Bytes | Meaning |
|---|---|---|
| 0 | 8 | ASCII `BMSHA001`, identifying this canonical stream version |
| 8 | 8 | Vector count |
| 16 | 8 | Source kind, 1 for raw and 2 for compressed BV |
| 24 | 8 | Source bytes consumed, prefix bytes for a prefix build |
| 32 | 32 | SHA-256 of consumed source bytes |

Each record stores the canonical vector SHA-256 followed by uint64 cardinality.
Raw and BV headers differ, while their record sections must match. Cached checks
validate the footer, count, exact file size and each decoded vector. They read
only the candidate stream/index and cache. Cache provenance identifies the input
used to construct it; the experiment runner establishes golden/raw agreement.
Publication uses a fresh `.partial` file followed by rename after successful
completion. Failures preserve the partial for diagnosis. Use one writer per fresh
cache path. Digest matching has the probabilistic assurance of SHA-256.

## Scheduling, memory and measurement

A coordinator reads the next batch while asynchronous workers process the active
batch. An atomic job index distributes vector tasks across the configured workers.
Every worker future is joined before the active slot is released or an error is
propagated. Results are written in original vector order. Worker failures name
the vector ordinal and candidate BLOB offsets. Each indexed BLOB must decode to
its exact boundary.

Both slots share one staging budget. Charged bytes include candidate/golden
encoded data, raw payloads, and 40 bytes of record data per staged vector. Normal
slots use half the budget. A single larger vector can use up to the full budget
after active work drains; exceeding the full budget fails before payload
allocation. Raw lists are staged under this bound before chunked import. Indices,
job metadata, decoded bvectors, allocator scratch, EVP contexts, stream buffers,
mapping page alignment and OS cache require additional memory. Per task, bvectors
and decoder workspaces are destroyed on success or exception. This prototype
creates a bounded set of async workers per batch; thread startup and allocator
costs are included in wall time.

Output reports total wall time, summed worker decode/import-plus-cardinality time,
summed equality/hash time, logical ID count, input bytes, and peak charged staged
bytes. These worker durations include scheduling delays and exclude some teardown
costs. Cache input bytes and index reads are additional to `input_bytes`. Cache
creation also hashes source bytes for provenance, and cached checks scan the small
cache once to validate its footer before reading records. The runner's wall time
includes process startup and index parsing. Run order and filesystem caching
influence timings, so single runs are observational measurements.

## Focused checks

```bash
python3 utils/inv_list/experiments/block_sha/test_block_sha.py \
  --bin /tmp/block-sha-build/block_sha --inv-bin /path/to/current/bminv
```

Python independently constructs canonical bytes and SHA-256 digests for sparse,
dense, FULL, high-ID, and holed vectors. Tests compare raw/BV cache records, exercise
buffered/mapped reads and multiple batches, validate after removing the golden
inputs, and induce later equal-cardinality ID substitutions, extra BLOB bytes,
cache corruption/truncation/count mismatches, raw truncation and oversized lengths.
Failed raw cache generation must leave no published cache. The same fixtures can
run with an ASan/UBSan build of the experimental executable.
