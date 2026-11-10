> Historical experiment. The opcode-retirement comparison runner and report
> generator have been removed. Use [current preset validation](../current/README.md).
> Findings and the shared-BWT C++ test remain available below.

# Retained shared BWT formats

Opcode 95 (packed BWT waves) and its encoder/decoder are retired. Opcode 96
retains the current original-polarity direct-BIC layout. Decoder-only support
for its earlier inverse-polarity and chain-formula variants was retired on
November 3; those unreleased streams require the older binary. Its current
writer bytes are unchanged, including fixed zero fields and the empty-formula
marker. The `--bwt-shared` selector still compares 72/96 using one
original-polarity BWT per block. Released formats remain supported.

Use [current highest-compression validation](../current/README.md) for full
sorted/unsorted verification and byte comparison against the latest full run.
The retained focused test is `test_bwt_shared.cpp`; it checks round trips,
profiling, pre-cleanup byte fixtures when supplied, and rejection of the
retired variants during materializing and skip traversal.

Historical measurements remain in [WAVE_SIZE_FINDINGS.md](WAVE_SIZE_FINDINGS.md),
[HISTORY.md](HISTORY.md) and [RETIREMENT_VALIDATION.md](RETIREMENT_VALIDATION.md).
Commands, counters and source references there describe their original
experiments and may no longer apply to the current code.
