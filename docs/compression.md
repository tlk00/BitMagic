# Compression levels

BitMagic separates its compact in-memory representations from the encodings
used to serialize them. This guide describes the current serializer presets;
a level chooses a policy for selecting representations, rather than a single
codec applied to every byte. All these encodings are lossless.

## Serialization of bvectors

### Choosing a level

`bm::serializer<BV>::set_compression_level()` selects a preset for a
`bm::bvector<>`. **Level 6 is the default.** Levels 1–10 are described below;
the API also accepts level 0, outside the compression presets covered here.

The table describes enabled techniques, not a guarantee that every input will
use them. Block contents, the in-memory block representation, eligibility
conditions, and candidate size determine the actual encoding. A higher level
does not mean that all data is converted to BWT or to any other one format.

| Level | Enabled techniques | Transform pipeline | Encoding selection |
|---|---|---|---|
| **1** | Basic framing, empty/full regions, fixed-width GAP boundaries, and BIT payloads. | Basic representation and framing. | Direct routing by block type. |
| **2** | Adds zero-word runs, hierarchical nonempty-region maps, and single-bit representation. | Lightweight zero-region representations. | Heuristic BIT-block selection. |
| **3** | Adds set-/clear-bit arrays for BIT blocks and BIC for existing GAP blocks. | Position/run-boundary lists → integer encoding. | Heuristics and block-type rules. |
| **4** | Delta-based Elias gamma coding for GAP boundaries and eligible position lists. | Differences → gamma. | Heuristics plus local size checks and fallback. |
| **5** | Broader BIC selection, sparse superblock grouping, and lighter BIC estimation settings. | Position lists → correction where applicable → BIC. | BIC-oriented cost models and fallback. |
| **6 — default** | Balanced BIC, broader sparse-group settings, and dynamic-range reduction; established formats. | Broader range-reduction → BIC policy. | Balanced cost models and fallback. |
| **7** | Adds BIT-array endpoint and superblock correction selectors; existing formats unchanged. | Same transform family as 6. | Adds complete-record endpoint and correction selection. |
| **8** | Adds GAP-XOR-BIC and adaptive, fixed/adaptive-wave, large-count, and inverse superblocks. | GAP cleanup/XOR exceptions; adaptive or wave corrections → BIC. | Adds size-based competition between eligible block/superblock formats. |
| **9** | Adds whole-block bitwise BWT and retained chain candidates; competing payload layouts share one BWT. | Adds BWT and chain pipelines, including chains after BWT. | Per-block selection: BWT must save **≥16 bytes** compared with the best non-BWT encoding of the **same 65,536-bit block**, including all record metadata and padding. |
| **10** | Same candidate family as 9. | Same pipelines as 9. | Same per-block selection, with a **≥1-byte** minimum BWT saving. |

Selection operates on individual blocks or eligible superblocks, depending on
the codec. The BWT saving threshold applies to one regular block—not to a
superblock or the entire `bvector<>`.

Levels are not ten strictly cumulative algorithms: in particular, level 4
changes the GAP policy, and levels 5–6 replace the lower-level BIT-block
selection rules. Presets 7–10 extend the higher-level selection policy.
The BWT thresholds compare complete records, including transform metadata,
framing, and alignment padding. Ties retain the non-BWT candidate. Level 10
is the most size-oriented preset, not a claim of globally optimal compression.

Levels 8–10 can emit newer formats that older BitMagic readers do not support.
Choose a level appropriate to the readers deployed by your application.
The default remains 6; levels 6–7 do not enable the newer opcode family.
The decoder follows the recorded format and does not need the encoder's level.

A **transform pipeline** combines reversible preprocessing with an encoder.
Simple pipelines already occur at lower levels (for example, differences
followed by gamma coding at level 4); levels 8–10 add richer compositions.
Transforms need not preserve the working representation's size: wave compaction
shrinks the coordinate domain but adds a restoration map. Both costs count.

**Candidate selection** is distinct from the pipeline itself. Heuristics
estimate costs; exact size scoring calculates encoded length without emitting
a candidate; trial encoding produces candidate bytes. Exact scoring and trial
encoding can make the same decision when both include all metadata and padding.
Level 7 expands complete-record selection, but lower levels already perform
some size checks and fallbacks. Selection is local to eligible candidates,
not an exhaustive search of every possible encoding.

Raw-storage fallback is also used by other compressors: [RFC 8878, section
3.1.1.2.2](https://www.rfc-editor.org/rfc/rfc8878.html#section-3.1.1.2.2)
recommends an uncompressed Zstandard block when compression would enlarge it.
This is an example of the selection principle, not a description of BitMagic's
wire format.

### API example

```cpp
#include "bm.h"
#include "bmserial.h"

void save_vector(const bm::bvector<>& bv,
                 bm::serializer<bm::bvector<>>::buffer& output)
{
    bm::serializer<bm::bvector<>> serializer;
    serializer.set_compression_level(10);

    // Optional application policy: require 32 bytes of BWT savings.
    // Omit this override to retain level 10's 1-byte threshold.
    serializer.set_block_selection_bwt_min_gain(32);
    serializer.serialize(bv, output);
}
```

Select the level **before** applying custom options. Calling
`set_compression_level()` again resets the preset's codec selectors and
parameters, including an earlier BWT threshold override. Invalid levels above
10 are ignored. The BWT gain setter requires block selection to be enabled
(as it is at levels 9–10); its accepted range is 0–65,535 bytes. A value of 0
still requires a strictly smaller record.

In-memory optimization and serialization level are separate controls.
`bvector::optimize()` can change a block's in-memory representation to GAP
when appropriate; selecting a serialization level is not a replacement for
that operation. The representation matters because GAP and BIT blocks follow
different encoder paths.

### Serialization samples

- [bvsample04](../samples/bvsample04/sample4.cpp): serialize to memory, manage buffers, and deserialize.
- [bvsample14](../samples/bvsample14/sample14.cpp): select a compression level and perform set operations on serialized BLOBs.
- [bvsample22](../samples/bvsample22/sample22.cpp): compress intervals and use bookmarks for selective range deserialization.

These examples illustrate the APIs; use the level table above for the current
presets rather than treating a sample's chosen level as a recommendation.

## Other containers

## Appendix: Methods and terminology

| Term | Brief explanation |
|---|---|
| **BIT block** | Bitmap covering 65,536 positions; stores individual bits in memory. |
| **[RLE — run-length encoding](https://en.wikipedia.org/wiki/Run-length_encoding)** | Represents repeated values as runs; for example, `111100011` has run lengths 4, 3, and 2. |
| **GAP block** | BitMagic's RLE representation: a starting bit value and alternating run-end positions, further encoded with fixed-width integers, gamma, or BIC. GAP is a name, not an acronym; GAP allocation levels differ from serializer levels. |
| **[Elias gamma coding](https://en.wikipedia.org/wiki/Elias_gamma_coding)** | Variable-length integer code favoring small positive values; zero requires a separate convention. |
| **[BIC — binary interpolative coding](https://pages.di.unipi.it/pibiri/papers/BIC.pdf)** | Compresses sorted integers by recursively encoding middle values within bounds implied by ordering. BitMagic uses an adapted implementation of BIC, combined with coordinate corrections and endpoint-layout selection, rather than an unchanged reference implementation from the literature. |
| **[Delta coding](https://en.wikipedia.org/wiki/Delta_encoding)** | Stores differences between successive values instead of the values themselves. |
| **Minimum-gap correction** | For minimum adjacent gap `d`, transforms `x[i]` to `x[i] - i * (d - 1)`, shrinking the range while preserving strict order. Unlike delta coding, the result remains a coordinate list. |
| **Adaptive correction** | Changes correction between windows of 16, 32, 64, or 128 list entries; records change flags and gamma-coded absolute values or signed differences for reconstruction. |
| **Endpoint layout** | Store endpoints explicitly to tighten BIC bounds, or derive bounds to avoid that metadata. Selection balances payload savings against endpoint storage costs. |
| **[XOR — exclusive OR](https://en.wikipedia.org/wiki/Exclusive_or)** | Toggles a bit wherever the corresponding mask bit is one; applying the same mask again reverses the operation. |
| **GAP-XOR-BIC** | BitMagic pipeline: remove isolated one-bit or zero-bit runs, encode the cleaned GAP base, and BIC-encode globally corrected exception positions. Original = base XOR exceptions; the ordinary GAP record remains the size fallback. |
| **Reference-XOR** | Encodes differences from a reference vector's block; distinct from GAP cleanup within one block. |
| **Superblock** | Group of 256 regular blocks (16,777,216 positions), encoded together to share headers and exploit spacing. Large-count variants can store set-bit or clear-bit lists. |
| **Wave / wave compaction** | BitMagic calls a fixed-size bit window a “wave.” Removing empty windows compacts the coordinate space before BIC encoding; a separately encoded window map restores the original positions. Savings must exceed map overhead. Superblock windows are fixed at 1,024 bits or selected from powers of two from 128 to 4,096. |
| **Polarity** | Choice of representing ones or zeros; for wave maps, occupied or empty regions. |
| **Hierarchical flags** | Group summaries represent uniform flag groups compactly; mixed groups retain individual flags. |
| **[BWT — Burrows–Wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform)** | Reversible rearrangement, not compression by itself. BitMagic transforms cyclic rotations of a 65,536-bit BIT block, one bit per symbol, preserving length and population count. Competing payload encodings share the transform; GAP blocks are not expanded for these trials. |
| **Primary index** | BWT metadata identifying the original sequence's row among sorted cyclic rotations. |
| **Periodic chain** | Predicted bits at `start + k * step`: step 1 is consecutive positions, step 2 skips every other position. Retained candidates include chains before or after BWT, with residuals for exceptions. |
| **Residual / exceptions** | Differences needed to reconstruct the original data from a simplified base or prediction. |
| **Interleaved start/end array** | Alternating chain starts and ends, corrected and BIC-encoded together with an explicit first start. Step-1 formulas choose between this layout and gamma coding. |
| **Opcode / record** | An opcode identifies an encoded format; the complete record includes its metadata, payload, and alignment padding. |
| **Complete-record gain** | Bytes saved after counting all record overhead, not just the compressed payload. |
| **Transform pipeline** | Reversible preprocessing followed by encoding; restoration metadata contributes to the total size. |
| **Candidate selection** | Chooses among eligible representations using heuristics, exact size scoring, or trial encoding. |
| **Exact size scoring / trial encoding** | Calculate candidate size without emitting bytes, or actually encode it; complete-record comparisons include metadata and padding. |

## References

These sources describe general algorithms and formats. BitMagic's adaptations
and selection policies are described above; they are not implied by the
standard algorithm names.

1. Michael Burrows and David Wheeler. [A Block-sorting Lossless Data Compression Algorithm](https://www.cs.jhu.edu/~langmea/resources/burrows_wheeler.pdf). Digital Systems Research Center, Research Report 124, 1994.
2. Giulio Ermanno Pibiri. [On Implementing the Binary Interpolative Coding Algorithm](https://pages.di.unipi.it/pibiri/papers/BIC.pdf).
3. [RFC 8878: Zstandard Compression and the application/zstd Media Type, section 3.1.1.2.2](https://www.rfc-editor.org/rfc/rfc8878.html#section-3.1.1.2.2). Raw-block fallback when compression would enlarge a block.
4. Wikipedia. [Run-length encoding](https://en.wikipedia.org/wiki/Run-length_encoding).
5. Wikipedia. [Elias gamma coding](https://en.wikipedia.org/wiki/Elias_gamma_coding).
6. Wikipedia. [Delta encoding](https://en.wikipedia.org/wiki/Delta_encoding).
7. Wikipedia. [Exclusive or](https://en.wikipedia.org/wiki/Exclusive_or).
8. Wikipedia. [Burrows–Wheeler transform](https://en.wikipedia.org/wiki/Burrows%E2%80%93Wheeler_transform).
