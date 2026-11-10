# Retained superblock wave encodings

Opcode **90**, `set_sblock_bienc_wave1024`, encodes fixed 1024-bit waves.
Opcode **97**, `set_sblock_bienc_wave_adaptive`, selects among
128/256/512/1024/2048/4096-bit waves using three width bits, ordinary BIC maps,
and exact dry scoring before writing the winning record once. The selector
builds 128-bit occupancy once and merges it for larger widths; corrected
positions still require a traversal for each width. Released formats and
68/92 fallback remain intact.

Both formats are retained. Exhaustive adaptive widths are intended for the
highest size-first compression levels, enabled explicitly, not a balanced
default. See [RESULTS.md](RESULTS.md) for full sorted/unsorted evidence.

The exploratory statistics and paired measurement harness are retired. Use
[../current/README.md](../current/README.md) and run_gov2_highest.sh for current
full verification, opcode reporting and reference byte-identity checks.
No opcode-97 width/savings instrumentation remains in the codec.

Focused tests (repository root; substitute each source):

```bash
clang++ -std=c++17 -O1 -g -fsanitize=address,undefined -DBMNEONOPT -Isrc \
  utils/inv_list/experiments/sblock_wave_adaptive/test_sblock_wave_adaptive.cpp \
  -o /tmp/test_sblock_wave_adaptive
/tmp/test_sblock_wave_adaptive
```

`test_sblock_wave_adaptive.cpp` checks exhaustive score/choice parity against
an independent test-only writer, all width/layout round trips, large map
indexes, boundary cases, invalid selectors, and actual opcode-97 profiling.
`test_sblock_wave.cpp` covers fixed opcode-90 layouts, fallback and integration.
The small independent model in test_model.h is test code only; it is not
included by bminv or installed serialization headers.

The fixed-wave API is `set_sblock_wave1024_options(true)`; its model is
`sblock_wave1024_model`. The CLI remains `-sblock-wave`. Renaming changes no
opcode values or payload bytes; fixed/adaptive sanitizer tests and synthetic
SHA-256 comparisons against pre-rename streams passed.
