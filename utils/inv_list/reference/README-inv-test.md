# Cross-platform inv_list reference tests

Run from `utils/inv_list`, using the binary built for the test platform:

```bash
bash run_inv_test_pl.sh --bin ./bminv_neon \
  --reference /Volumes/USB_1TB/gov2-reference-Nov-2026 \
  -j 8 --out /Volumes/USB_1TB/inv-test-full-$(date +%Y%m%d-%H%M%S)
```

| Mode | Work for every sorted/unsorted level 1–10 |
|---|---|
| Normal | Check saved reference integrity; threaded and serial raw-ID verification; independent decode; rebuild into `--out`; repeat both verifications and decode on rebuilt streams. |
| `--fast` | Full threaded raw-ID verification only. Reads every ID; skips checksum/history audits, rebuilding and decode-only passes. |

The raw inputs must be `<reference>/Gov2Flat/gov2.sorted` and `gov2.unsorted`,
Normal mode checks their recorded input hashes; fast mode directly tests decoding
against the supplied raw inputs. The current bminv verifier imports each raw ID
list into a bvector and calls `equal()` on the decoded vector. Both saved and
rebuilt streams are thus compared to the same logical expected vectors.
Serialized size/checksum changes are reported, never treated as logical failures.
An older target binary may use its own equivalent per-ID verification algorithm.

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
