# Canonical-block SHA-256 results

Measured on 2026-11-11 on macOS 26.6.2 ARM64, Apple Clang 21, native BitMagic NEON and OpenSSL 3.6.4. Both full Gov2 datasets passed all phases. Each dataset contains 1,107,205 vectors and 5,979,715,441 IDs.

The runs used eight workers, buffered reads, a shared 1,024 MiB staging budget and a 2,048-vector batch cap. The coordinator overlaps preparation with worker decoding/checking. These are single sequential observations with filesystem caching effects. Some fixture/build activity overlapped the sorted run. The two builds used the same algorithm and optimization, with direct Clang compilation for sorted and CMake Release for unsorted. Binary checksums and commands are preserved in the run directories.

| Phase | Gov2.sorted seconds | Gov2.unsorted seconds |
|---|---:|---:|
| exact-6-vs-10 | 41.66 | 43.25 |
| sha-6-vs-10 | 97.12 | 185.08 |
| build-bv-cache | 43.65 | 90.73 |
| build-raw-cache | 53.45 | 102.04 |
| cached-level-7 | 43.86 | 91.15 |
| cached-level-8 | 51.92 | 91.36 |
| cached-level-9 | 53.15 | 91.55 |
| cached-level-10 | 54.18 | 95.60 |

All raw-derived cache records matched the corresponding level-6-derived cache records byte-for-byte. Each cache is 44,288,296 bytes, approximately 42.24 MiB. Source SHA-256 provenance in the cache headers matches the saved reference raw-input and level-6 stream hashes. Cached levels 7–10 read their candidate BLOBs/index and the small cache, with no golden BV or raw reads.

The canonical-block approach reduced candidate-validation input traffic but increased CPU checking time enough to lose elapsed time in these runs. On Gov2.sorted, exact pair checking spent 2.90 summed worker-seconds in equality, while cached level 10 spent 251.16 summed worker-seconds in hashing/checking. On Gov2.unsorted the corresponding values are 4.88 and 611.49. Worker durations include scheduling delays and differ from process CPU time.

For sorted level 10, the pair pass read 5,350,850,871 BLOB bytes; the cached pass read 2,644,146,026 BLOB bytes. Cache validation and record reading add approximately 84.5 MiB, plus index reads. The same canonical 8 KiB is hashed for each occupied block regardless of its population, so sparse occupied blocks can expand the hashing workload substantially. These measurements establish a CPU/I/O tradeoff on this machine. HDD throughput and Windows hardware behavior require separate measurement.

All reported peak charged staged bytes were within 1,073,741,824 bytes. That charge includes both slots and record storage, while indices, decoded bvectors, codec/hash scratch, allocator metadata and OS cache are additional.

Independent Python hashlib/bit-placement oracle tests passed on sparse, dense, FULL, high-ID, holed and empty vectors. Buffered/mapped checks, large raw vectors, exact/SHA/cache mismatch paths, removed-golden cached validation, cache corruption/count/truncation checks, raw truncation and oversized-length rejection passed. The fixtures also passed with ASan/UBSan. A warning-enabled syntax check passed. Whole-Gov2 timing runs used buffered reads; whole-Gov2 mmap timing and Windows execution were not performed.

## Artifact cleanup

The temporary run directories, caches and experimental builds were removed at
the user's request during integration into inv_list. The recorded measurements
above remain available. Integrated use is documented in
`../../reference/README-inv-test.md`.
