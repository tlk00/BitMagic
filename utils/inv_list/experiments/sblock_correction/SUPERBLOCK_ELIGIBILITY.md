# Future experiment: higher-cardinality superblock routing

User proposal, 2026-11-02: allow more IDs in a superblock at explicitly selected
highest compression levels. This is separate from the current opcode-68
playlists, which contain only groups accepted by the existing routing policy.
No routing or format change is implemented by this note.

Implementation follow-up: the separately optional opcode-98 experiment now
implements the agreed two higher count ranges; see [LARGE.md](LARGE.md).
The existing count-limited routing described below remains the baseline.

## Current code review

- `serializer::serialize_to()` in `src/bmserial.h` surveys each top-level group
  once using `blocks_manager::is_sparse_sblock()` in `src/bmblocks.h`.
- Address span is fixed: 256 blocks of 65,536 bits, or 16,777,216 bit positions.
  Raising the ID cap is different from expanding that span.
- `sparse_max_l5=48`, `sparse_max_l6=256`. The survey's count gate is
  `count <= sparse_cutoff * 256`: level 5 therefore allows at most 12,288 IDs.
  Level 6 permits 65,536 through that gate, but `serialize_to()` additionally
  requires `bv_count < 65536`, making the effective cap **65,535 IDs**.
  The current code does not impose a roughly 30K cap.
- `set_sparse_cutoff()` caps its parameter at 256, so it cannot raise the
  current highest limit. Grouping requires more than one nonempty ordinary
  block; groups below five IDs are rejected. Reference XOR disables
  sparse-group routing. Compression level must be at least 5.
- GAP tolerance defaults to the relaxed 100 policy: it does not reject eligible
  groups simply because their original GAP representation is efficient. Lower
  tolerances can retain the former GAP-length rejection rule.
- Once count/type routing succeeds, the serializer commits to a superblock and
  skips its constituent blocks. The enabled selector compares whole-group
  formats 68/92/90/97, but this routing point does not trial the complete
  alternative of separately encoding all constituent blocks.

## Bounds to audit before any expansion

Opcode 68 uses an 8-bit count or `delta16` count. Current 90/97 decoders enforce
`n <= 65535`; adaptive 92 decoding has the same upper bound. Encoder reserve
calculations, adaptive scratch allocation and stream-reader record-size bounds
also assume approximately 65,536 entries. Expanding the count gate alone is
unsafe. Preserve existing opcode-68 encoding/decoding compatibility; evaluate a
new optional experimental format or explicitly version an unreleased format
only after deciding which codecs actually warrant larger counts.

## Proposed staged measurement

Survey rejected groups above 65,535 IDs across full sorted and unsorted inputs,
using cardinality bands and complete constituent-block bytes. Stratify by
byte contribution, block types, number of nonempty blocks and largest block
contribution. The current opcode-68 playlists cannot cover this population.
Then compare full trial superblock records against the existing separate-block
encoding, including all framing, indices, endpoints and padding. Enable wider
grouping only on a strict win, initially in explicitly selected highest
compression modes. Keep address span fixed for this experiment; changing the
span is an additional independent design dimension.

## User refinement: opcode-implied count base

A dedicated higher-count ("2x") opcode can imply `count >= 65,536` and store
only `count - 65,536` in an unsigned 16-bit field. This represents counts
65,536 through 131,071 inclusive without widening that field or changing the
existing opcodes. Count 131,072 would require another range/escape or a different
endpoint convention; define that explicitly if exact doubling to that inclusive
endpoint is desired. This is count-range expansion within the existing address
span, not a doubling of the span itself.

Thus length-field widening is not inherently necessary. The new decoder must
restore the implicit count base before checking bounds, allocating scratch and
decoding. Grow encoder and decoder scratch storage dynamically from the decoded
or measured entry count, using checked size calculations and preserving input
validation. Audit fixed record bounds and helper restrictions as well as the
main array. Keep this option reserved for the future high-compression experiment;
it is not part of the current unsplit/quarter playlist control.
