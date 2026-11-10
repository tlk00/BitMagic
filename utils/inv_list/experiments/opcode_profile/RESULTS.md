# Full Gov2.unsorted opcode efficiency, 2026-11-02

The retained highest-compression configuration encoded and raw-verified
1,107,205 vectors containing 5,979,715,441 original IDs. Independent decoding
and opcode byte/cardinality reconciliation passed. The stream is
5,440,396,028 bytes, or 7.278468 weighted bits/original ID. SHA-256 matches
the previous verified highest-compression unsorted baseline:
`bd4d71607c3bcb9fa1ea5afd248f2d4c43d561f8ad1cff1391d90fb8bf067c9a`.

| Opcode | Name | Bytes | Stream share | Original IDs | Total bits/ID |
|---|---|---:|---:|---:|---:|
| 65 | set_block_arr_bienc_v3s | 2,606,384,553 | 47.908% | 3,670,095,480 | 5.681344 |
| 68 | set_sblock_bienc_v3 | 2,211,498,269 | 40.650% | 1,467,976,630 | 12.051954 |
| 61 | set_block_gap_bienc_v3 | 257,384,689 | 4.731% | 211,745,675 | 9.724295 |
| 94 | set_block_gap_xor_framed | 132,778,824 | 2.441% | 116,671,846 | 9.104429 |
| 72 | set_block_bwt | 130,990,166 | 2.408% | 250,759,757 | 4.178985 |
| 11 | set_block_bit | 64,569,033 | 1.187% | 214,924,969 | 2.403407 |
| 97 | set_sblock_bienc_wave_adaptive | 20,232,955 | 0.372% | 17,861,352 | 9.062228 |
| 96 | set_block_bwt_unpacked | 7,997,803 | 0.147% | 15,393,524 | 4.156451 |
| 90 | set_sblock_bienc_wave1024 | 144,159 | 0.003% | 108,826 | 10.597394 |

Opcodes 65 and 68 account for 88.5576% of stream bytes, making their payloads
the principal size targets. Top-level markers plus vector headers account for
only 0.128755%; this does not measure internal metadata or padding separately.

The denominator is original set-bit positions restored by the record, including
inverted/transformed formats. Total bits/ID includes the marker and all payload
metadata/padding. Payload bits/ID excludes only the marker. Zero-ID records show
n/a. These are weighted ratios of totals. Different opcode populations have
different densities and ranges; lower bits/ID alone is not evidence that one
codec would beat another on the same input. In particular, raw BIT blocks are
dense, whereas opcode 68 covers sparse superblocks.

Encode wall time including profiling was 1,379.215 seconds; raw verification
63.739 seconds; independent decode 57.544 seconds. Cardinality profiling adds
a full decode pass, so this is not an isolated encoder timing comparison.

Evidence: `/Volumes/USB_1TB/BitMagic-experiments/2026-11-02/runs/inv_list/gov2-unsorted-opcode-bpi-20261102-091640`

Decision: retain generic opcode bytes/IDs/bits-per-ID collection in BM
serialization code, presentation in bmdbg/harness, and reconciliation checks.
The unsuccessful small-wave experiment and preceding occupancy scanner are
retired; their [findings](../block_small_wave/RESULTS.md) remain available.
GAP wave-map coding was proposed but never measured or implemented.
