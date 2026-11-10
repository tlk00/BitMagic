# GAP decomposition experiment

Build `bminv` separately from `utils/inv_list` with `make neon`. The experiment
is disabled unless `-gap-xor` is passed. Existing GAP encodings and opcodes are
unchanged. The new opcode 94 wraps a cleaned ordinary GAP component and one XOR
residual stored as bounded BIC positions. Each full encoded block competes
with its ordinary GAP encoding, and ties use the ordinary encoding. Experimental opcodes 90 and 93 and their decoders have been removed. Their
numbers are available for reuse. Opcode 94 is the sole GAP-XOR format; the
`-gap-xor` option enables it directly. Released BitMagic formats are unchanged.

The runner compresses, verifies against the raw IDs, decodes, and writes a
summary of complete stream bytes, timing, and candidate/winner counts.
`--reference-dir` reuses a previously verified control stream and
reads its first N BLOB offsets, which avoids writing a second 10K stream.
Generated streams are retained unless `--discard-streams` is specified. Raw
Gov2 inputs remain in the repo; output goes to USB_1TB. Run with `--help` for
all options. The focused test covers BIC residual decoding, ordinary/operation/stream
round trips, indexed skipping and ten workers:

```
clang++ -std=c++17 -O1 -g -fsanitize=address -DBMNEONOPT -I../../src \
  experiments/gap_xor/test_gap_xor_codec.cpp -o /private/tmp/test_gap_xor_codec
/private/tmp/test_gap_xor_codec
```

## Full Gov2 results that motivated BIC-only emission

All 1,107,205 vectors in each collection verified and decoded. The original
three-mode GAP-XOR search saved 3,502,266 bytes (0.130182%) on Gov2.sorted and
760,672 bytes (0.013965%) on Gov2.unsorted against the matching BIT-selector
control. Among emitted GAP-XOR blocks, the residual-mode counts were:

| Collection | Gamma | BIC | Nested GAP | BIC share |
| --- | ---: | ---: | ---: | ---: |
| sorted | 12,912 | 651,186 | 53,676 | 90.72% |
| unsorted | 0 | 305,578 | 186 | 99.94% |

These are selection frequencies, not an isolated BIC-versus-gamma byte
comparison. BIC-only emission simplifies search and may lose some of the
original saving; the new full Gov2.sorted run measures that loss directly.
The old complete-run reports remain on USB_1TB under
`gov2-sorted-gap-xor-full-20261027` and `gov2-unsorted-gap-xor-full-20261027`.

## Boundary encoding follow-up

Legacy GAP v3 already stores the upper endpoint in reverse coordinates as
`tail_delta = 65535 - max_v`. Its current payload is an 8- or 16-bit fixed
field, sometimes with three low bits in the header. A gamma trial should encode
`tail_delta + 1`, because gamma has no zero value. Compare the complete record,
including a mode bit and the final bit-writer flush, against the current field.
For a fixed 8-bit tail, gamma is shorter only for sufficiently small values;
for a 16-bit tail it has a wider opportunity. The exact break-even depends on
header reuse and word padding, so a dry bit count alone is insufficient.

The new residual BIC mode derives its upper bound from count and minimum delta:
`65535 - (n-1)*(min_delta-1)`. It stores no separate maximum. The legacy
cleaned GAP component and any nested GAP residual each use their own encoder
and correction coordinates, so a shared upper boundary for three lists is not
currently a free metadata saving. A later opcode could measure shared bounds
after mapping all three arrays into compatible coordinates; this experiment
keeps the old component formats intact for compatibility.

## Retired adaptive GAP-XOR experiment

The 10K sorted/unsorted comparison found compact framing savings of 4,504
and 2,772 bytes. Adaptive correction added only 503 and 6 bytes, respectively.
Both obsolete formats (90 and 93) were uncommitted experiments, not released
formats requiring compatibility. Only compact BIC framing (94) remains.
The superblock adaptive codec and its 64-entry minimum are unchanged.

## Exception coding hypothesis rejected (2026-10-29)

The evenly sampled 10K sorted/unsorted probe verified and decoded both
streams. On selected GAP-XOR records, raw gamma saved only 4/0 bytes,
corrected gamma 16/0, and explicit-bound BIC 1,380/120. These were diagnostic
projections including metadata and padding, not emitted alternate streams.
The hypothesis of a useful compression improvement was rejected for these
samples; retain the existing BIC format. There were 55 sorted candidate
splits potentially rescued by corrected gamma, not 55 distinct block wins.
No unsorted candidates were rescued. The diagnostic code and runner were
removed. Reports remain at
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-29/runs/inv_list/gov2-gap-exception-20261029-184350`.
