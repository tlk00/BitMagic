# Opcode-95 retirement validation

Opcode 95, the unreleased packed-wave BWT experiment, has been removed from
both encoding and decoding. Its map/scatter implementation, wave-size probes,
polarity/codec comparison modes, scratch statistics, CLI flags, and obsolete
experiment runners/tests were removed. Historical findings were retained.

The retained shared 72/96 selector still performs one original-polarity BWT
per eligible block. Opcode 96 (`bm::set_block_bwt_unpacked`) is implemented by
`bm::bwt_bic_encoder` / `bm::bwt_bic_decoder` in `src/bmbwt.h`. Its wire
layout, including reserved header bits, is unchanged. Its reader retains
existing polarity and formula support. Released formats, chain encodings,
GAP-XOR, and superblock formats were not retired or renumbered. Opcode 95 is
unassigned; streams containing that unreleased experiment are no longer supported.

## Focused checks

- Standalone `make neon` completed.
- Focused native-NEON ASan/UBSan test passed with dense/sparse patterned and
  random input: round trips, profile counts, one-BWT accounting, and exact
  equality to fixtures captured with the pre-cleanup implementation.
- Opcode-96 formula and inverted-polarity decoding passed explicit tests;
  materializing and skip-only decoding consumed exactly the encoded record.
- Script syntax, Python compilation and path-scoped whitespace checks passed.
- Source scan confirms removal of opcode-95 dispatch and retired options.

## Paired 10K verification

Run: `/Volumes/USB_1TB/BitMagic-experiments/2026-10-31/runs/inv_list/gov2-bwt-retired-10k-20261031-190434`.
Reference: `gov2-bwt-waves-10k-20261031-174949/baseline` under the same runs folder.

Both sorted and unsorted prefixes use 10,000 vectors and 10 workers, with
adaptive/wave superblocks, GAP-XOR, shared BWT/chain selection, the 16-byte BWT
cutoff and opcode-65 endpoint selection. The harness requires raw-ID
verification, independent decoding, opcode accounting, no opcode-95 records,
and SHA-256 identity with the previous verified streams.

Both datasets passed. New streams are byte-identical to the verified reference
by SHA-256, with zero opcode-95 records and matching transform/profile counts.
Each prefix contains 4,818,820,706 IDs.

| Dataset | Vectors | Serialized bytes | Bits/ID | Encode s | Verify s | Decode s |
|---|---:|---:|---:|---:|---:|---:|
| Sorted | 10,000 | 1,457,318,527 | 2.419378 | 691.048 | 10.683 | 7.314 |
| Unsorted | 10,000 | 3,561,956,968 | 5.913409 | 949.554 | 11.838 | 7.386 |

Sorted selected 128,906 opcode-72 records and 70,669 opcode-96 records, with
634,155 BWT transforms for 634,155 eligible blocks. Unsorted selected 27,423
opcode-72 records and 1,651 opcode-96 records, with 1,362,922 transforms for
1,362,922 eligible blocks.

SHA-256 hashes (identical before and after cleanup):

- Sorted: `0bb703a6ca07503ef77d99b3329eb9af785bfedcc54de7c5b870487caa7d81b2`
- Unsorted: `d016ffc7feb7b72463dd2794a3d370f87e3588e73c3139591593e8fec2001b04`

The run folder contains `retirement-report.txt`, per-dataset reports, logs,
indexes, and hashes. Verified generated `.bv` files were removed by the
harness after verification, independent decoding, accounting and hashing.
Raw inputs remain local. These are 10K checks, not a full stress regression.
Single-pass timings are reported for context, not as a speedup claim.
