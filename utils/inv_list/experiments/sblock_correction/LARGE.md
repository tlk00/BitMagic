> Historical findings: the large/polarity and quarter-trial runners were retired
> during release cleanup on 2026-11-10. Commands below are archival.

# Optional positive/inverse superblocks (opcode 98)

`set_sblock_bienc_large` is an experimental format, disabled by default through
`serializer::set_sblock_large_options(false)`. inv_list enables it explicitly
with `-sblock-large`; the highest-compression harness forwards `--sblock-large`.
Default serialization takes the existing route. Old opcode-68 wire encoding
and decoding remain unchanged.

## Format

The logical address span remains 256 bit-blocks (16,777,216 positions).
Opcode 98 now begins with one polarity bit: 0 lists set positions; 1 lists
zero positions. This revises the unreleased opcode-98 format without a decoder
for its previous positive-only layout. Released opcode formats are unchanged.
Opcode 98 reuses opcode-68-style endpoint, address and optional global minimum
correction fields, followed by `bic_encode_u32_cm` interiors. All fields share
one bit accumulator. The former length-width header bit is the new range bit;
the positive count field contains an unsigned 16-bit residual:

| Range bit | Count | Residual |
|---|---|---|
| 0 | 65,536..131,071 | count - 65,536 |
| 1 | 131,072..196,607 | count - 131,072 |

Inverse count has no base correction. Its range bit selects the width:
0 stores counts 0..65,535 in 16 bits; 1 stores larger counts in 24 bits.
This reuses an existing header bit, adding no extra selector field. Counts 0/1/2 omit the BIC interior;
the existing endpoint fields remain present. Inverse search considers only
groups with more than half their positions set, and at most the configured
number of zeros (default 196,607). `set_sblock_large_inverse_limit(N)` /
`-sblock-large-inverse-limit N` controls search effort; 0 is the positive-only
control. The wire decoder supports counts beyond this search cap.

Inverse restore sets runs of ones between the zero positions, intersected
with a requested range if present. It preserves OR deserialization semantics
when the destination already contains bits and never sets the reserved final
ID. Original-ID profiling counts restored ones, rather than listed zeros.
Decoder buffering uses the actual count in the fixed prefix. Heap scratch
also retains room for subsequent legacy superblock records.

The optional global correction is minimum adjacent gap minus one, clipped to
65535, as for ordinary 68. No adaptive correction or wave transforms were added
to this format. Decoder checks count, endpoint span, correction restoration and
address bounds before allocating the ID vector and invoking BIC. Stream record
bounds include the larger maximum. Scratch arrays and candidate buffers grow on
the heap; no count-sized stack arrays are introduced.

## Selection and scope

The experiment first constructs the exact ordinary stream with all enabled
selectors, then walks its deserialization index. It trials eligible groups and
accepts a new record only if the complete replacement is smaller. This compares
against actual constituent-block encoding, including BWT/chains, GAP-XOR,
opcode-65 selection, run framing and word padding. A zero/full run crossing the
group end is retained as a residual run with its bytes included in the trial.
Ties preserve the original bytes. Groups whose leading set bits are already
owned by a preceding cross-boundary full run are conservatively skipped and
counted as `boundary_skips`; this initial experiment does not split that prefix.

This post-selection implementation prioritizes an exact fallback over encoding
speed: it adds a baseline index pass and temporary BLOB storage on eligible
vectors. Reference XOR, bookmarks and destructive serialization are explicitly
excluded for this experiment. Input vectors below 65,536 IDs bypass the index
pass. All current ordinary superblock alternatives remain available in the
baseline, but opcode 98 only targets groups beyond their existing count limit.

The report reconciles accepted opcode-98 counts/bytes and predicted savings
against the final stream. Displacement counters show the baseline spans
replaced (a crossing run can also retain a residual run); the ordinary opcode
profile reports actual final-stream composition. Existing transform search
work still occurs before replacement, so candidate-search savings counters
must not be added to high-count savings. Final shared-BWT selection counts
are adjusted for displaced records and retain their reconciliation checks.

## Recursion and validation

CM BIC halves the array on its recursive branch and iterates over the other
branch. Maximum depth at 196,607 entries is approximately 18 levels, versus
16 at the old cap. This was tested, not assumed to eliminate all stack risks:
the focused suite passes on a 256 KiB pthread stack, in 32-/64-bit address
builds, and under ASan/UBSan on the normal test stack. Focused cases include
both range endpoints, just-outside counts, exact fallback, disabled byte
identity, end-of-address-space framing, mixed old/new records, crossing zero
runs, streaming output, bounded decode/index/profile, ranged decode and
truncated input. This is not a full stress-suite completion.

## Initial measurement

The first 1,000 vectors of the high-byte unsorted playlist (not the full
collection prefix) passed raw verification, independent decoding and exact
opcode/savings reconciliation. Baseline was 161,140,094 bytes; selection was
158,990,625 bytes: 2,149,469 bytes saved (1.333913%). Range 0 selected 919
records and range 1 selected 74, with no boundary skips. This targeted smoke
sample does not establish whole-collection savings or a stable speed result.

Evidence: `/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-large-smoke-20261102/report.txt`.

A full unsorted paired run was launched at
`/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-unsorted-large-full-20261102`.
Its final `report.txt` is produced only after both streams pass verification
and independent decoding; launch alone is not a validation result.

The full unsorted run completed successfully at 18:44 on 2026-11-02:
1,107,205 vectors and 5,979,715,441 IDs. Baseline 5,440,396,028 bytes;
selected 5,429,693,606 bytes; saved 10,702,422 bytes (0.196721%).
Range 0 selected all 4,715 eligible groups (8,056,972 bytes saved); range 1
selected 1,790 of 1,928 (2,645,450 bytes saved). No boundary skips occurred.
Opcode 98 represents 716,476,707 IDs in 700,139,965 bytes (7.817588 bits/ID).
Raw verification, independent decoding, opcode accounting and exact savings
reconciliation all passed. Encode/verify/decode seconds were
1356.940/53.371/48.623 baseline and 1378.629/67.588/60.184 selected.
This paired measurement shows a 1.60% encode-time increase and 23.78%
decode-time increase; a repeated isolated timing run is needed before treating
the latter as a stable cost. The size improvement is real but substantially
smaller than the targeted smoke sample suggested.

## Preserved positive-only baseline and polarity experiment

The complete verified pre-polarity measurements, opcode composition and SHA-256
are frozen locally in `OP98_POSITIVE_BASELINE.json`, sourced from
`gov2-all-op98-20261102-191605`. Sorted: 2,644,388,979 bytes, 3.537812 bits/ID,
6,150 opcode-98 records / 367,440,266 opcode bytes, 14,485,865 bytes saved
against all best options without 98. Unsorted: 5,429,693,606 bytes,
7.264150 bits/ID, 6,505 records / 700,139,965 bytes, 10,702,422 bytes saved.
Each dataset contains 1,107,205 vectors and 5,979,715,441 original IDs.

`run_polarity.py` runs two arms with the revised format: inverse cap 0
(positive control) and cap 196,607 (polarity selection). Both arms use every
retained highest-compression selector on full sorted and unsorted inputs.
This isolates inverse gains from the extra leading bit and word padding in
the revised positive format. Every stream is raw-verified, independently
decoded and opcode-profiled before deletion. The final table reports sizes,
weighted bits/ID, encode/verify/decode timings, inverse eligibility/wins,
over-cap groups and exact savings reconciliation. It also compares final
sizes with the frozen previous positive-only format.

Build separately; sources/scripts remain local, outputs go to the USB:

```bash
make neon
python3 experiments/sblock_correction/run_polarity.py \
  --bin ./bminv --threads 10 --progress 10000 --inverse-limit 196607 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-op98-polarity-$(date +%Y%m%d-%H%M%S)"
```

No `--max-vectors` means both full collections. The driver checks USB mount
and free space. Individual opcode reports are retained under
`MODE/DATASET/profile/report.txt`; the final paired table is `report.txt`.

Focused polarity validation passed under ASan/UBSan and in a 64-bit-address
build on a 256 KiB pthread stack. Cases cover inverse counts 0/1/2, both
positive count boundaries, 262,144-entry inverse wire decoding, mixed legacy
records after a short inverse list, final-address-space framing, existing
destination OR semantics, indexed ranges, truncation and cap fallback.
The complete paired harness passed with 10 threads on a synthetic dense
fixture (one inverse winner in each dataset arm; exact 4,335-byte saving)
and on the first real vector of each Gov2 dataset. The real-vector smoke had
no inverse winner: its dense group exceeded the initial zero-count cap.
These focused runs do not replace the requested full-collection experiment.
Evidence: `/Volumes/USB_1TB/BitMagic-experiments/tmp-op98-polarity/integration/report.txt`
and `/Volumes/USB_1TB/BitMagic-experiments/tmp-op98-polarity/real-smoke/report.txt`.

### Original design note (now implemented as the optional trial)

User design note (2026-11-02): measure whether the new
superblock approach also benefits very dense groups by encoding the positions
of zero bits, with an explicit polarity flag, and restoring the complement
within the exact superblock span. Compare complete records against all current
selectors. Collect original and complemented cardinalities, candidate/winner
counts, byte savings and decode costs before choosing a format.

The original opcode 98 represented only 65,536..196,607 listed IDs. A very dense
group can have fewer than 65,536 zeros, so its complemented list falls outside
that format. There is no symmetric inverted-superblock counterpart in this
experiment. A future format should cover these lower list counts as well as
the higher ranges; the current range-bit plus residual cannot express them.
Evaluate a 24-bit count (possibly count minus one), or a compact low-count
mode/escape that preserves short metadata for small complemented lists.
24 bits can represent counts through 16,777,215 directly; the complete
16,777,216-position span needs count-minus-one or an explicit full-group case.
Empty/full cases, endpoint omission for lists of size 0/1/2, partial final
groups and range/index decoding require explicit treatment. A polarity mode
alone does not require materializing a dense list: buffers should follow the
chosen sparse side. Preserve released opcode formats; stage any new layout
as an optional experiment with strict complete-size fallback.

## Run

Build separately, from utils/inv_list:

```bash
make neon
python3 experiments/sblock_correction/run_large.py \
  --bin ./bminv --dataset unsorted --threads 10 --progress 10000 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-unsorted-large-$(date +%Y%m%d-%H%M%S)"
```

The default is the full collection. Optional `--max-vectors N` limits a smoke
run; `--input-dir DIR` can replay extracted playlists. Baseline and selected
streams run sequentially, each followed by raw-ID verification, independent
threaded decode and opcode profiling. Verified `.bv` streams are removed by
the existing child harness; logs, indices, SHA-256 and reports remain. Sources
stay local, outputs go to a freshly created USB directory. The driver checks
USB availability/free space and prints a paired summary table at the end.

## Retain inverse mode for other datasets (2026-11-03)

The full Gov2 polarity run (`gov2-op98-polarity-20261102-202756`) had no
inverse candidates within the 196,607-zero cap and one over-cap group per
dataset. The first-vector trial with cap 8,388,607
(`gov2-op98-dense-probe-20261102-223531`) evaluated both groups and accepted
neither; all streams verified and matched their positive controls by SHA-256.
The user chose to retain inverse encoding for future datasets. Its existing
range bit now selects 16/24-bit zero counts; positive count ranges keep their
current meaning. The superseded fixed-24-bit inverse experiment is unreleased.

Encoding enumerates a temporary zero mask (`group mask -= original vector`).
Decoding sets the one-runs between decoded zero positions. This has the same
reconstructed value as `full group AND NOT zero mask`, while preserving the
normal OR-merge behavior for a pre-populated destination. Clearing zero
exceptions directly in the destination would incorrectly erase existing bits.

Count-width revision validation: ASan/UBSan boundary, mixed-record, indexed
range, OR-merge, truncation and streaming tests passed, including inverse
1,024- and 65,536-zero groups. Native NEON inv_list rebuild passed. The denser
fixture exposed an ordinary full-block run lookahead beyond the allocated top
block table; a bounds guard fixes that read without changing the wire format.
