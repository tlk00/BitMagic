# Gov2 validated reference matrix

Build separately from `utils/inv_list`:

```bash
python3 reference/build_reference.py
```

This force-builds one release executable using native ARM NEON or supported
x86-64 AVX2. It records the build command, compiler, source identity and binary
checksum. It also rebuilds `bminv`; the reference executable has a separate name.
No benchmarks run during building.

Generate full sorted/unsorted collections at levels 1–10:

```bash
bash reference/run_gov2_reference.sh \
  --out /Volumes/USB_1TB/gov2-reference-Nov-2026
```

| Setting/content | Meaning |
|---|---|
| `-j N` | Compression and threaded-reader workers; default `max(1, logical CPUs - 2)`. Serial verification uses one reader. |
| `--levels 1,6,10`, `--datasets sorted` | Optional staged subsets; defaults cover both datasets and all levels 1–10. No vector-count limit: every case covers the full input file. |
| `--input-dir`, `--volume` | Override the raw-input directory and required mounted output volume, including Linux mount paths. |
| `--reserve-gib 8` | Minimum free-space reserve; checked before each case and while child phases run. Next-case planning uses raw input bytes as an estimate, not a guaranteed size bound. |
| Final `.bv` / `.idx` | Published only after threaded AND serial raw-ID comparison, independent decode, and opcode/index accounting. |
| `reports/<dataset.level>/` | Logs, commands, timings, opcode TSV, report, and checksum manifest. |
| `references.tsv` | Validated entries with absolute sizes, weighted bits/ID, ratios, per-phase wall times, workers/platform and build provenance. |
| `.tools/`, `matrix.json` | Frozen optimized executable and input/build identity, reused on later invocations. |
| `.partial/` | Incomplete cases; retained on failure. Inspect and move a failed case aside before retrying; never delete validated entries to make room. |

Rerun the same command to resume. Every committed artifact checksum and validation
record is checked before skipping. A different input or explicitly requested binary
is rejected. A directory lock prevents overlapping writers; after a crash, inspect
`.reference-lock/owner.json` and confirm that the owner process is gone before
removing the stale lock directory. Completed manifests are not overwritten.

The script creates a short directory `readme.md`. Inputs and executable sources
remain unchanged by execution. The frozen executable and reference artifacts are
retained on the output volume. Timings include profiling, validation and I/O;
concurrent system load and thermal conditions affect comparisons.

Small-fixture checks after building (no Gov2 or USB writes):

```bash
python3 reference/test_reference.py ./bminv_reference_neon
# On x86: use ./bminv_reference_avx2 instead.
```

For testing another platform binary against this immutable reference matrix,
see [run_inv_test_pl.sh usage](README-inv-test.md). Fast mode still compares
every raw ID; normal mode also rebuilds and tests candidates outside this tree.
