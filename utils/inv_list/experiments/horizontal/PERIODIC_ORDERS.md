# Ordered periodic-chain experiment


USB storage (USB_1TB is an exFAT volume mounted at `/Volumes/USB_1TB`): the Gov2 input files stay in the repository. Gov2 experiment output directories are stored under `/Volumes/USB_1TB/BitMagic-experiments/2026-10-23/`, with links at their former paths. The USB wrapper checks that the volume is mounted before starting; the core harness also checks every path under `/Volumes` so an unplugged drive cannot silently redirect output to the system disk.

To continue the interrupted 10K sample-128 run on the USB, from `utils/inv_list`:

```bash
bash /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list/experiments/horizontal/run_gov2_periodic_orders_usb.sh --resume
```

This reuses the complete sorted BLOBs and `unsorted.baseline`, then discards and retries the incomplete `unsorted.independent` BLOB. Each remaining BLOB is written to the USB. Add `--prune-verified` if space becomes tight; that option removes each BLOB only after verification and decoding, retaining logs and recorded sizes. For a new USB run use the same wrapper without `--resume`; pass `--sample 1` if evaluating every eligible BIT block.

Build the inv_list app separately:

```bash
cd /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list
make
```

Run both Gov2 distributions, first 10,000 vectors, sampling about 1/128 eligible BIT blocks:

```bash
bash /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list/experiments/horizontal/run_gov2_periodic_orders_usb.sh \
  --bin ./bminv \
  --sorted /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.sorted \
  --unsorted /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.unsorted \
  --out /Volumes/USB_1TB/BitMagic-experiments/2026-10-23/runs/inv_list/gov2-periodic-orders-first10000-s128 \
  --max-vectors 10000 --progress 25 --threads 8 --batch-vectors 512 \
  --sample 128 --min-step 1 --max-step 128 --budget 2304 \
  --max-chains 256 --formula-chains 256 --min-length 4 --holes 1 --level 6
```

If a run fills the disk, resume the same directory with identical experiment settings and `--resume --prune-verified`. The harness verifies and decodes each complete BLOB before removing it, writes its byte count to a `.size` sidecar, and keeps all compression logs for the final report. A partial BLOB is discarded and recompressed; its verification and decode are repeated. Changing the sampling rate, step range, or other experiment settings requires a new output directory. For example, the first 10K run with sample 128 can be recovered by adding `--resume --prune-verified` to its original command. The pruning option is also useful from the beginning of a new run.

The script performs no build. It refuses an existing output directory and stops on any failed compression, verification, decoding, or report-accounting check. Without pruning it retains ten BLOBs: five modes for each dataset. `--sample 1` evaluates every eligible BIT block and can be much slower. GAP blocks remain outside the experiment.

Modes:

* Baseline: ordinary serialization, no periodic search.
* Independent: all step/hole searches use the original polarity-normalized block; choose the best single-step subset.
* Dense: ascending steps, updating the residual after each accepted addition.
* Sparse: descending steps, updating the residual.
* Middle: 64,65,63,66,...,1,128 for the default range, updating the residual.

Polarity is normalized once. At each step, both hole policies search independently on the current residual. Retain up to `--max-chains` chains per search using the existing pool policy. Rank retained chains by net popcount reduction (length minus twice the holes), then length and start for deterministic ties. Evaluate prefixes of size 1,2,4,... up to the remaining formula capacity; include the actual pool/capacity endpoint. Each subset is sorted by start for BIC.

The best strict improvement in complete residual-plus-accumulated-formula bytes is committed once per step. Ordinary encoding wins ties. This greedy policy cannot cross an intermediate size regression to reach a better later combination. The total formula is capped by `--formula-chains`, independently of pool size. Saturation can stop a schedule early; selected-count histograms expose this. The default 2,304-trial budget covers 128 steps times two hole policies times nine subset sizes. Smaller budgets can bias schedule coverage; exhaustion is reported. The old `-periodic-variants` option applies only to the legacy engine.

The new formula uses opcode 76 inside the existing residual frame. Its shared bitstream contains polarity (1 bit), group count minus one (7 bits), and each group's step minus one (7 bits), small-count flag plus count minus one (5 bits for up to 32, otherwise 8), BIC starts, and gamma(length minus three). There is one terminal 32-bit flush for all groups. Residual framing still costs four bytes, and residual codec padding remains separate. Existing opcode-75 formulas remain readable. No Rice, residue-reordering, or adaptive length-code selection is included in this first comparison.

Outputs:

* `report.txt`: measured sizes/times, selected chain/group counts, schedule wins and density/GAP strata.
* `metrics.tsv`: per-mode summaries.
* `blocks.tsv`: identical sampled coordinates compared across all four experimental modes, including original normalized popcount and GAP count.
* `chain_lengths.tsv`: discovered, retained, and accepted chain-length histograms. Lengths 1–32 are exact bins; larger bins double in width. Discovery counts are observations across searches, including both hole policies and evolving residuals, not distinct chains.
* Per-phase logs and timings, and `run.txt` with exact commands.

The calculated per-block best of four is an oracle comparison, not another emitted BLOB. Incremental step savings sum to actual stream savings. Formula bit totals include header, starts, lengths and padding; fixed framing is reported separately. All modes use the same compact formula, so the comparison isolates independent versus ordered subset searches; it does not isolate compact coding from subset selection relative to the older experiment.

Optional focused correctness tests (32-/64-bit identifiers, AddressSanitizer):

```bash
bash /Users/anatoliykuznetsov/dev/BitMagic/utils/inv_list/experiments/horizontal/build_periodic_order_tests.sh \
  /Users/anatoliykuznetsov/dev/BitMagic /tmp/bm-periodic-order-tests
/tmp/bm-periodic-order-tests/order-test32
/tmp/bm-periodic-order-tests/order-test64
```

Tests cover count escapes, step 128, 128 groups, overlapping XOR predictions, polarity, exact consumption, truncated formulas, all schedules, budget limits, multi-group acceptance, round trips and ordinary fallback. These are focused tests, not the library stress suite or a completed 10K benchmark.
