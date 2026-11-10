# Gov2.sorted GAP level-2/3 forced-BWT experiment (completed)

This temporary experiment reconstructed every raw Gov2.sorted vector as an
optimized `bvector<>`, converted each GAP level-2/3 block to a BIT block in
memory, and forced a whole-block BWT record. Native BIT blocks used the
level-6 step-2 selector with a 16-byte BWT minimum gain and the adaptive
gamma-versus-fixed-interleaved step-1 formula. The control used the same
selector on native BIT blocks and retained ordinary GAP blocks.

The full run completed after recovering and verifying an interrupted
8,125-vector prefix. Both complete streams verified against the raw IDs and
decoded all 1,107,205 vectors (5,979,715,441 IDs). The forced stream
converted 1,053,101 GAP level-2 and 959,483 level-3 blocks.

| Stream | Bytes | Bits/original ID |
| --- | ---: | ---: |
| Control | 2,690,277,862 | 3.599205 |
| Forced GAP2/3 BWT | 2,695,162,898 | 3.605741 |

Forced BWT was 4,885,036 bytes (0.181581%) larger. Its measured full-stream
decode took 95.407 seconds versus 43.718 seconds for the control, using
eight read workers; verification took 101.309 versus 48.817 seconds. These
are single-run timings. The interrupted encode and resumed tail do not give
a controlled full-collection encode-time comparison.

Conclusion: forcing BWT for *all* high-level GAP blocks did not improve this
collection. This result does not rule out selective per-block wins. A future
probe should compare complete native-GAP and BWT records separately for
levels 2 and 3, including framing and a meaningful gain threshold, before
adding a selector.

The generated streams, indexes, phase logs, and report remain at
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/gov2-sorted-gap23-forced-bwt-20261027-093511`.
The experiment-only serializer flag, conversion path, CLI option, resume
support, and runnable harness have been removed from the checkout.
