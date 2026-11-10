# Opcode survey and reusable vector playlists

Survey full Gov2.sorted and Gov2.unsorted using all retained highest-compression
selectors, excluding experimental opcode-68 correction selection. Optimize each
reconstructed input vector before serialization. Existing BM workers handle
compression and raw-vector round-trip checks; the main thread writes per-vector
statistics and aggregates the profiler. Batches remain bounded to 512 vectors /
512 MiB input (one oversized vector can exceed that limit). The complete
collection is never retained in RAM or written as a compressed file: verified
BLOBs are discarded via /dev/null after profiling. Sources and scripts stay local;
only measurements, playlists and replay inputs/outputs go to USB.

## Build and generate playlists

From utils/inv_list:

```bash
make neon gov2-playlist-neon
SURVEY_DIR="/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-op68-survey-$(date +%Y%m%d-%H%M%S)"
python3 experiments/opcode_playlist/run_survey.py \
  --bin ./bminv --tool ./gov2_playlist \
  --input-dir /Volumes/USB_1TB/int_set_bench/Gov2Flat \
  --threads 10 --progress 1000 --opcode 68 --top 10000 --min-bytes 256 \
  --out "$SURVEY_DIR"
```

Both full surveys run in sequence. Every serialized BLOB is independently
restored and compared against its original optimized input vector. Profiling
also reconciles original-ID counts, marker/payload/header bytes and full stream
size. Source size/mtime must remain stable. A COMPLETE marker is written only
after successful process completion. No separate full serialized collections or
full verification input copies consume USB space.

Each `DATASET/opcodes.csv` has sparse per-vector/opcode rows with original
zero-based vector index, raw byte offset, vector cardinality and total serialized
size, opcode count, complete record bytes and original IDs represented. Header
and legacy ID-list bytes use pseudo-opcodes -1/-2. Keep this full table: other
opcode playlists and thresholds can be produced without rerunning compression.
The target-eligible tuples are held in memory during ranking, not the vectors.

Eligible vectors have at least 256 complete target-opcode bytes and nonzero
target IDs. Rank by descending exact `target_bytes / target_IDs`, then descending
target bytes, then ascending original index. This is a ratio of totals across
all target records in the vector. It is not a mean of per-record ratios.
The default filter tests target-opcode contribution, not whole-vector bytes.
If fewer than N qualify, select all eligible vectors and report the actual count.

`playlists/gov2.DATASET.op68.high-bpi.min256B.top10000.bv` is an ordinary BM
serialized bvector containing original vector indices, verified by round trip.
The JSON records the ranking, baseline binary hash/options, full survey totals,
source identity, counts and artifact hashes. The ranked CSV preserves score
order and exact source-record SHA-256 hashes. The bvector supplies membership;
replay uses ascending original index order. Playlists are ranked independently
for sorted/unsorted. For other opcodes or sizes, rerun `playlist.py select` with
new parameters (fresh semantic names; existing playlists are not overwritten).

Survey report.txt shows full opcode sizes/efficiencies and the playlists' byte
contributions, cardinalities and weighted bits/ID. survey_summary.json preserves
all per-opcode full-collection efficiencies, including sorted Gov2.

## Evaluate both playlists

To complement the sparse high-bits/ID cases, select by largest total target
opcode byte contribution using the existing archived survey:

```bash
SURVEY_DIR="/Volumes/USB_1TB/BitMagic-experiments/surveys/gov2-op68-survey-20261102-122955"
python3 experiments/opcode_playlist/playlist.py select \
  --survey "$SURVEY_DIR" --tool ./gov2_playlist \
  --opcode 68 --top 10000 --min-bytes 256 --ranking high-bytes
```

This reads the existing measurements and hashes selected raw input records;
it does not repeat compression. Ranking is descending target bytes with ties
resolved by ascending original index. New artifacts use `high-bytes` names;
the original high-bpi playlists and selection summary are preserved. The new
summary is `selection-summary-high-bytes.json`. Pass `--ranking high-bytes`
to run_replay.py to evaluate these playlists. For this survey their paired
raw replay inputs require approximately 4.32 GiB, plus comparison outputs.

Reuse SURVEY_DIR from the selection command:

```bash
python3 experiments/opcode_playlist/run_replay.py \
  --survey "$SURVEY_DIR" --bin ./bminv --tool ./gov2_playlist \
  --threads 10 --progress 1000 --opcode 68 --top 10000 --min-bytes 256 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-op68-playlist-eval-$(date +%Y%m%d-%H%M%S)"
```

This compares baseline and opcode-68 correction selection on each dataset's
own frozen playlist: four streams, raw-ID verification and independent threaded
decode. It reuses the existing correction harness. The optional --input-dir
relocates raw inputs but still requires matching identity and selected-record
hashes. Original inputs are never moved. Extracted inputs and index mappings
stay in OUT/inputs; compressed experiment streams are deleted after successful
verification/reporting by the existing harness. Logs, indices, hashes and
reports remain. Original indices live in inputs/gov2.DATASET.mapping.csv and
playlist metadata; internal compressed collection indices run from 0 in replay
order. Current search options evaluate every block, so replay order does not
change per-vector choices.

The replay baseline must exactly reproduce the surveyed selected target bytes,
IDs and record counts. This detects configuration drift or an incompatible
binary before interpreting candidate savings. Final report.txt compares complete
stream savings, correction decisions and opcode efficiencies with the frozen
full survey for both datasets. These deliberately difficult workloads are for
algorithm development; their gains must not be extrapolated to the collection.

## Validation

`python3 experiments/opcode_playlist/test_playlist.py` tests exact ratio ranking,
tie rules, zero/high index membership, duplicate rejection, truncation and extra
stream bytes. Focused profiler ASan/UBSan tests remain in opcode_profile.
A ten-vector distributed Gov2 fixture exercised full-survey CSV accounting,
paired top-three playlist creation, exact membership restore, raw extraction,
and all four replay encode/verify/decode phases with frozen-baseline matching.
This is focused validation, not execution of the full requested survey.

## Persistent survey library

Completed, reusable survey datasets belong under
`/Volumes/USB_1TB/BitMagic-experiments/surveys/`, outside dated disposable
experiment runs. Preserve survey CSVs, frozen playlists, ranked measurements,
metadata and reports here. Replay benchmark outputs remain in the dated runs
tree. Raw Gov2 inputs and code stay in their original local locations.

Do not move an active survey. Once run_survey.py finishes successfully:

```bash
python3 experiments/opcode_playlist/archive_survey.py --survey "$SURVEY_DIR"
SURVEY_DIR="/Volumes/USB_1TB/BitMagic-experiments/surveys/$(basename "$SURVEY_DIR")"
```

The helper requires completed paired survey and selection reports, validates
playlist/CSV hashes, refuses an existing destination, and uses a rename on the
same USB filesystem. It does not copy large data or change frozen artifacts.
relocation.json records provenance. Replay normalizes its report paths to the
actual --survey location, while retaining the original creation path. Relative
playlist/CSV references remain valid. Reusable surveys should be excluded from
routine disposable experiment-stream cleanup.
