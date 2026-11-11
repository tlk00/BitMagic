# Cross-platform inv_list reference tests

Run from `utils/inv_list`, using the binary built for the test platform:

```bash
bash run_inv_test_pl.sh --bin ./bminv_neon \
  --reference /Volumes/USB_1TB/gov2-reference-Nov-2026 \
  -j 8 --out /Volumes/USB_1TB/inv-test-full-$(date +%Y%m%d-%H%M%S)
```

| Mode | Work for every selected sorted/unsorted level |
|---|---|
| Normal | Audit saved reference integrity; validate golden against raw IDs; compare other streams against golden with workers; retain serial raw-ID verification, decode and rebuild checks. |
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
threaded verification. All other threaded checks, including rebuilt candidates,
decode both collections and compare exact bvector equality. A golden failure
stops the run before dependent comparisons for that dataset.

Fast mode omits checksum/history audits, rebuilding and decode-only passes. Summaries report the actual case count; `test.json`
records the level range and count. Normal mode also accepts `--start-level` and
retains its reference integrity audits, including the final all-file checksum scan.

The raw inputs must be `<reference>/Gov2Flat/gov2.sorted` and `gov2.unsorted`,
Normal mode checks their recorded input hashes; both modes validate the golden
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
logs, candidates, indexes and process working directories are under `--out`.
The script never chmods, removes, or overwrites references. Normal mode checks
all saved checksums before testing and reference-file checksums again on success.
Fast mode omits these audits entirely; basic matrix/file discovery checks remain.
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
