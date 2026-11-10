> Historical findings: the large/polarity and quarter-trial runners were retired
> during release cleanup on 2026-11-10. Commands below are archival.

# Offline four-quarter trial

`quarter_trial.cpp` reads raw vectors, optimizes them, serializes using the
retained highest selectors, verifies each baseline vector, and profiles final
opcode-68 records. It trials only these final records, so their baseline already
won competition against 90/92/97. Other opcodes are unchanged in the projected
collection size. No production opcode, writer, or decoder is added or modified.

The hypothetical record contains one opcode byte, gamma(superblock index+1),
four quarter-occupancy bits, then each occupied quarter's payload in address
order. Quarters divide the 16,777,216-bit address range into four equal ranges.
Each payload stores gamma(count), one correction-presence bit, optional
gamma(correction), gamma(first+1), and (for at least two IDs)
gamma(quarter-width-corrected-last). Interior corrected values use
`bic_encode_u32_cm`, bounded by corrected endpoints. Singleton quarters store
only the first endpoint. All fields share one accumulator; it flushes once
at record end. Local correction is min adjacent gap minus one, clipped to
65535 like opcode 68. It is selected only when its exact payload bit count
including correction metadata is smaller than the uncorrected payload.
Every candidate, including losers, is decoded and compared against all original
relative IDs, and full byte consumption is checked. Strict complete-record
savings are summed; ties fall back to the existing encoding.

The no-local-correction control uses the identical quarter framing but disables
all corrections. Its difference from adaptive quarters isolates local correction
benefit. Differences versus opcode 68 also include changed boundary metadata
and quarter partitioning; those effects must not be attributed solely to
min-delta correction. CSV length buckets are ceilings of log2(original count).

Build from the repository root:

```bash
c++ -std=c++17 -O2 -DBMNEONOPT -Isrc \
  utils/inv_list/experiments/sblock_correction/quarter_trial.cpp \
  -o /tmp/gov2_quarter_trial
/tmp/gov2_quarter_trial --selftest
```

Run separately for each extracted raw playlist input:

```bash
/tmp/gov2_quarter_trial PATH/inputs/gov2.sorted PATH/sorted.csv >PATH/sorted.log 2>&1
```

The scanner is single threaded. Its elapsed time includes original encoding,
profiling, baseline decoding, candidate trials and candidate validation; it is
not a production encoder benchmark. This is a format exploration, not a complete
new-opcode integration or whole-collection verification.

## Unsplit control

The next control uses the same payload fields on the entire unsplit superblock:
gamma(count), correction flag/value, gamma(first+1), gamma(span-last), and CM
BIC interiors. It stores one opcode byte and gamma(superblock index+1), with
no quarter mask. It measures unsplit-only wins and how much quarters add after
both ordinary opcode 68 and this control are available. Distinct hypothetical
opcodes identify these two layouts, so no additional layout bit is assumed.
Each layout and their joint best choice preserve the existing encoding on ties.
The scanner now disables byte-order/GAP-length headers exactly like bminv and
uses its sorted bulk-insertion path. Original quarter results used one additional
byte-order header byte per vector; opcode-68 record costs were unchanged.
