> Historical experiment notes. The exploratory measurement code and old
> runners described below were retired on 2026-11-01. Use
> [current highest-compression verification](../current/README.md).
> Retained codec correctness tests remain where applicable.

# Opcode-65 endpoint layout measurement

Build `bminv` separately with `make neon` in `utils/inv_list`. Then run the
harness in that same directory:

```bash
bash experiments/block_endpoint/run_gov2_block_endpoint.sh \
  --bin ./bminv --datasets unsorted --threads 10 --progress 1000 \
  --out /Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-unsorted-endpoints-$(date +%Y%m%d-%H%M%S)
```

The default is the full Gov2.unsorted collection. `--max-vectors 10000` runs
only the first 10K vectors. `--datasets sorted|both` is also supported. Raw
input remains in the repository; output and temporary streams stay on the USB.
The harness verifies every input ID, decodes the collection, checks the opcode
profile, and writes a summary table to `block-endpoint-report.txt`. Verified
`.bv` files are deleted unless `--keep-streams` is specified.

The `-block-endpoint-trial` option observes only FINAL emitted opcode-65
records. It decodes their position lists and independently writes two scratch
records with the exact existing wire layout: full-range BIC and explicit first
and last endpoints. For the emitted mode, its re-encoded bytes must match the
original record exactly. Complete candidate bytes include opcode, header, BIC,
and final 32-bit padding. Strictly smaller alternatives count as wins; ties
preserve the emitted mode. No production writer, decoder, or opcode is changed.

The report splits wins by the current heuristic's decision and gives aggregate
header/BIC/padding costs, per-record gain buckets, and paired-trial time. It
reconciles record count and baseline bytes with the independent top-level
opcode profile. Projected savings do not represent emitted compression gains.

A focused round-trip and byte-identity test is in `test_block_endpoint.cpp`.

## Optional encoder selector

The original `-block-endpoint-select` implementation wrote both existing opcode-65
endpoint layouts to bounded scratch storage, compared complete record bytes,
and retained the smaller one. Ties keep the original heuristic choice. The
option is disabled by default and leaves the opcode-65 decoder unchanged.

Run `run_gov2_block_endpoint_select.sh` for paired baseline and selected
streams with the established BWT, GAP-XOR, 68/92/90 superblock settings.
Both are independently verified and decoded. The paired endpoint probe should
find zero missed wins in the selected opcode-65 records; the final report
compares complete collection bytes and timings. Timings include identical
post-serialization endpoint diagnostic passes in both runs.

The focused `test_block_endpoint_select.cpp` checks round trips for both
streams and confirms that every emitted opcode-65 record already uses the
smaller of its two existing layouts.

## Compression-level design note

Keep this exact two-layout selection available as a size-first option when
defining future compression levels. On the full Gov2.unsorted collection,
the paired verified run saved 1,388,448 bytes (0.025514%) and increased
single-run encode wall time from 949.888 to 1,006.761 seconds (about 6%).
Both runs included the same diagnostic probe. The opcode-65 decoder and wire
format did not change; the selector remains disabled by default. These data
support a higher-effort size-oriented setting, but do not alone justify making
it a general default. Re-evaluate on other workloads and with a cheaper
choice rule before assigning it to a compression level. Full run:
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-unsorted-endpoint-select-20261030-210708`.

On the full Gov2.sorted collection, the same paired configuration saved
839,538 bytes (0.031495%), from 2,665,619,612 to 2,664,780,074 bytes.
Encode wall time rose from 664.094 to 724.104 seconds (about 9%) in one run;
both sides again included the same endpoint diagnostic pass. All 1,107,205
vectors and 5,979,715,441 IDs verified and decoded. The baseline's emitted
opcode-65 records had 500,424 bytes of endpoint-only opportunity, while the
selected stream had none. Actual collection savings were larger because the
smaller array candidates changed top-level codec choice: opcode 65 gained
14,453 records and BWT lost 9,985 records. Decode wall times were 51.509 and
51.868 seconds, too close in a single run to infer a speed effect. The sorted
result supports retaining the mode as a size-first candidate, but its encode
cost still argues for a cheaper choice rule before making it a balanced
compression-level default. Full run:
`/Volumes/USB_1TB/BitMagic-experiments/2026-10-30/runs/inv_list/gov2-sorted-endpoint-select-20261030-224741`.

## Exact bounded dry scoring (2026-11-01)

The retained selector now writes the heuristic layout, counts the alternative
using the same templated u16 BIC traversal, and writes the alternative only for
a strict complete-record byte win. Metadata and final 32-bit padding are charged.
Counting stops when a strict win becomes impossible. Ties retain the heuristic.
The old per-serializer 64 KiB trial buffer is removed. The decoder and wire
format are unchanged, and selection remains optional.

`bench_endpoint.cpp` is a single-thread serializer-only driver. It reads raw
Gov2 records, constructs and optimizes the vector outside the timer, times
serialization, separately times decoding, checks every round trip and writes
concatenated serialized records. `benchmark_dry.py` runs separately built old
and new implementations in ABBA order, reports serialization and wall times,
and requires identical stream SHA-256 hashes. It intentionally omits BWT and
other optional searches to isolate endpoint selection. Build old/new binaries
with the same compiler flags against the corresponding complete header states.
Do not rebuild or run other CPU-intensive tests during the timings.

Example from `utils/inv_list` after separately building both drivers:

```bash
python3 experiments/block_endpoint/benchmark_dry.py \
  --before experiments/block_endpoint/bench_before \
  --after experiments/block_endpoint/bench_after \
  --input /Volumes/USB_1TB/int_set_bench/Gov2Flat/gov2.unsorted \
  --max-vectors 10000 --out /Volumes/USB_1TB/BitMagic-experiments/op65-dry-10k
```

For full validation with all retained compression options, use
`experiments/current/run_gov2_highest.sh --datasets both --threads 10`.
