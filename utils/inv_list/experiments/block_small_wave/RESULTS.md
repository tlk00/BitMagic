# Retired: opcode-65 small-wave packing

Retired on 2026-11-02 after a full Gov2.unsorted trial. No production opcode
was introduced. The exploratory scanners, scratch codec, tests, reporting
scripts, execution harnesses and build targets have been removed, including
the preceding 8/16-bit occupancy scanner. Generic opcode size/cardinality
profiling is retained independently.

## Measurement

Full collection: 1,107,205 vectors, 5,979,715,441 original IDs,
5,440,396,028 serialized bytes. Final opcode 65 represents 1,322,869 records,
3,670,095,480 original IDs and 2,606,384,553 bytes (5.681344 bits/original ID).

The occupancy scan found 74.123% empty aligned 8-bit waves and 58.428% empty
16-bit waves. Every examined opcode-65 block had empty waves. Occupancy alone
did not predict a compression gain.

The trial removed empty waves, remapped retained positions, and encoded a BIC
wave map plus a BIC residual. It tried occupied/empty map polarity and both
residual endpoint layouts. Complete costs included marker, metadata and final
32-bit accumulator padding. Adaptive width additionally paid one selection bit.
No extra minimum-delta correction was applied.

| Variant | Eligible records | Strict wins | Candidate bytes before fallback | Overhead vs opcode 65 |
|---|---:|---:|---:|---:|
| Fixed 8-bit | 1,322,869 | 0 | 2,945,512,697 | 339,128,144 |
| Fixed 16-bit | 1,322,869 | 0 | 2,821,746,429 | 215,361,876 |
| Adaptive 8/16-bit | 1,322,869 | 0 | 2,821,911,617 | 215,527,064 |

Map overhead outweighed residual savings. The trial establishes failure of
these layouts on final opcode-65 Gov2.unsorted records. GAP-coded maps and
additional correction methods were not measured; they remain hypotheses,
without retained implementation. Occupied-wave adjacency/run statistics were
not collected, so occupancy histograms cannot establish GAP map size.

Raw-ID verification, independent decode and byte accounting passed. Scratch
cost/write/restore logic also passed synthetic positive fixtures and sanitizer
tests. Trial scoring took 34.475 seconds with 10 workers; no alternative
collection stream was emitted.

Evidence retained on USB:

- Occupancy: `/Volumes/USB_1TB/BitMagic-experiments/2026-11-01/runs/inv_list/gov2-unsorted-tiny-waves-20261101-211700`
- Trial: `/Volumes/USB_1TB/BitMagic-experiments/2026-11-01/runs/inv_list/gov2-unsorted-small-wave-trial-20261101-223344`

See [opcode efficiency results](../opcode_profile/RESULTS.md) for the retained
profiler and full collection measurements.
