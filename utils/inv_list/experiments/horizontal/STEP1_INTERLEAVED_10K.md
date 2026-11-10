# Gov2.sorted step-1 interleaving: 10K format comparison

The fixed-minimum interleaved format is the smallest of the three single-format
streams. The experiment compressed the first 10,000 Gov2.sorted vectors at
level 6 with ten NEON workers, 512-vector batches, independent residual step 2,
and a 16-byte BWT gain threshold. All four streams, including the reused
explicit-first gamma control, were verified against the source and decoded.
There were 643,520 eligible BIT blocks in each stream.

| Formula | Complete BV bytes | Saving vs gamma | Compression | Verify | Decode |
|---|---:|---:|---:|---:|---:|
| Explicit-first gamma control | 1,467,475,175 | — | — | — | — |
| Fixed interleaving | 1,466,934,194 | 540,981 (0.036865%) | 627.112 s | 9.991 s | 6.524 s |
| Shared measured minimum | 1,467,051,296 | 423,879 (0.028885%) | 614.508 s | 9.750 s | 6.489 s |
| Separate measured minima | 1,467,063,722 | 411,453 (0.028038%) | 614.535 s | 9.750 s | 6.594 s |

Fixed is 117,102 bytes smaller than shared and 129,528 bytes smaller than
separate minima. The complete-record comparison shows 54,801 BIT blocks whose
selected record improved by 700,443 bytes and 29,940 that regressed by 159,462
bytes relative to gamma. These are small net gains; the single-pass timing
differences are insufficient to claim a speed ranking. The previous four-mode
adaptive run was 1,466,732,118 bytes, 202,076 bytes smaller than fixed-only,
but required all three interleaved codecs and their decoders. This cleanup
accepts that measured size cost to keep only gamma and the best standalone
interleaved format. The new two-choice adaptive result remains to be measured
on the full collection.

The code now retains fixed interleaving (`-periodic-step1-interleaved 1`) and
per-formula gamma-versus-fixed selection (`... 4`). Modes 2 and 3, along with
experimental wire opcodes 86–89 and their decoders, were removed. Previously
written streams using those opcodes require the earlier binary. The remaining
experimental options are still disabled by default.

The detailed 10K report and phase logs are in
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/gov2-sorted-step1-variants-v10000-20261025`.
Its three verified `.bv` and `.idx` pairs were deleted after the report was
saved to free USB space. Raw Gov2 input remains in the repository checkout.

From `utils/inv_list`, build separately with `make neon`. Run the full
collection on the mounted USB with a fresh output directory and the existing
verified explicit-first control:

```bash
bash experiments/horizontal/run_gov2_step1_interleaved.sh \
  --bin ./bminv \
  --input /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.sorted \
  --control-dir /Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/gov2-sorted-step1-first-start-20261025-152910 \
  --out "/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/gov2-sorted-step1-gamma-fixed-full-$(date +%Y%m%d-%H%M%S)" \
  --mode auto --threads 10 --batch-vectors 512 \
  --level 6 --min-gain 16 --progress 1000
```

The runner checks the mount and free space, compresses, verifies, decodes, and
reports complete stream savings against the reused gamma control.
