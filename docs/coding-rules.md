# BitMagic coding rules

This file is the source of truth for the project's coding conventions collected
from development discussions. Read the current version before changing code.
Update rules here instead of maintaining duplicate rules in agent memory.
Apply them to new and modified code; keep unrelated refactoring out of a change.

## Documentation and declarations

- Document headers with their purpose, active use, and relevant limitations.
  Distinguish production code from experimental or development-only utilities.
- Give each class and structure a Doxygen description. Mark internal utilities
  explicitly with `@internal` and use the appropriate Doxygen group (for example,
  `@ingroup bvserial` for serialization utilities).
- Document every function and method: intent, every parameter, return value when
  applicable, preconditions, state changes, and ownership/lifetime requirements.
- Explain implementation-specific algorithms. For example, BitMagic BWT works on
  bits; describe how it differs from classic symbol-oriented BWT.
- When documenting an opcode or any other constant, include both its numeric 
  value and symbolic name, for example opcode 96 (`bm::set_block_bwt_unpacked`). 
  Verify both against the code.
- Put each data member on its own line with a short Doxygen annotation:

  ```cpp
  unsigned count = 0;  ///< Number of events recorded by the serializer.
  ```

- Keep documented declarations in classes/structures and move implementations
  outside them. Header definitions that need `inline` must retain it.
- Put `inline` on a separate line, after any template declarations:

  ```cpp
  template<class Alloc>
  inline
  void codec_workspace<Alloc>::init()
  {
      // Reset state while retaining storage.
  }
  ```

  Function templates generally do not require `inline` solely to permit header
  definitions. Non-template functions and full specializations still need the
  appropriate linkage treatment. Follow the surrounding header's convention.

## Formatting and names

- Use readable spacing around operators and after control keywords.
- Put control-structure braces on separate lines. Expand even short branches;
  use one statement per line rather than compressed statements.

  ```cpp
  if (!ones && n < 128)
  {
      enc.put_8(static_cast<unsigned char>(128 + n));
      return;
  }
  ```

- Comment loop endings with the counter variable, and switch endings with the
  selector. For other loops, use a meaningful controlling variable or condition.

  ```cpp
  for (unsigned i = 0; i < count; ++i)
  {
      process(i);
  } // for i

  switch (sb_mode)
  {
  case 0:
      out.gamma(sb + 1);
      break;
  default:
      out.put_32_no(sb);
      break;
  } // switch sb_mode
  ```

- Explicitly qualify BM types, constants, and standalone function calls with
  `bm::`, including inside namespace `bm`. Examples: `bm::bwt_wave_bits`,
  `bm::set_sblock_bienc_large`, and `bm::bwt_sort_wave_rotations(...)`.
- Use the `Vect` suffix for vector-type template parameters, such as
  `CorrectionsVect` or `ValuesVect`.
- Use semantic names for buffers. When sharing storage across algorithms,
  choose a suitably general name that describes all its uses.

## Containers, allocation, and workspace lifetime

- Avoid STL dependencies in low-level (algorithmic and SIMD kernels) library headers 
  and preserve non-STL use.
  Prefer BM containers and existing BM algorithms over `std::vector` and
  `<algorithm>` utilities. Use BM sort utilities; keep general-purpose utilities
  such as heap sort next to the existing sorting functions.
- use of STL is fine in experimental code in diagnostics and rendering code like bmdbg.h 
- use of STL is ok hight level code (examples and utilities) and in code implementing 
  parallel algorithms with threading dependencies and mutexes
- Propagate the application's allocator. Derive it from `BV::allocator_type`
  through the serializer or pass an allocator template parameter. Do not hardcode
  `bm::standard_allocator` in internal storage that must honor a custom BV
  allocator, including debug allocators.
- Prefer retained serializer/deserializer workspaces over local heap containers
  constructed and destroyed for every block. Avoid large temporary blocks on
  the stack; use allocator-backed, correctly aligned storage.
- Keep scratch capacity across calls, but bound growth by the codec's needs.
  Reuse storage rather than accumulating memory for each vector processed.
- Check whether temporary vectors can reserve or initialize enough storage
  before repeated work. Lazy allocation is acceptable. Where justified, document
  a speculative initial capacity and its rationale; do not impose arbitrary
  eager allocations on every algorithm.
- Make state reset explicit with an `init()` method and ensure it is called
  before use. Reset all per-operation state, including early-return paths, while
  retaining capacity. Document any configuration or cached state intentionally
  preserved by a reset.
- Prefer existing optimized BM block operations over handwritten loops, such as
  block inversion and BIC-to-bitset decoding. Meet SIMD alignment requirements
  using the BV allocator before selecting such operations. 

## Structure and dependencies

- Separate producer and consumer responsibilities into encoder/decoder or
  serializer/deserializer classes. Use descriptive names rather than a generic
  model that mixes both interfaces where that separation is appropriate.
- Keep closely related implementations together when this removes an
  unnecessary header split. Internal function-grouping helpers may use templated
  structures directly in `bm` rather than extra case-specific namespaces.
- Reuse common validation utilities such as `bm::serialization_require()`;
  avoid codec-local copies with identical behavior. Use generic serialization
  errors for shared validation, not names tied to a retired transform.
- Put report rendering and presentation helpers in `bmdbg.h`. Keep models and
  compression logic independent of reporting. Template report arguments where
  needed so unrelated users do not need every serialization type declared.
- Keep development statistics and observation histories optional. Avoid their
  allocation and processing overhead in ordinary application use. Propagate
  the policy through all serializer types, not only bvector serialization.
- Include required headers explicitly; do not depend on incidental transitive
  includes. Update build dependencies when retiring or merging headers.

## Correctness and validation

- Add `BM_ASSERT` for meaningful internal invariants and bounds. Keep runtime
  validation for invalid input; assertions must not replace required error
  handling. Test release builds must evaluate BM assertions, with and without
  SIMD. Ordinary release applications, examples, and performance builds should
  not enable internal assertions solely for test purposes.
- Respect deserialization's OR semantics when decoding into an existing bvector.
  Cover zero, nonzero, GAP, and full target blocks where applicable.
- New opcodes must be classified correctly for indexing, bookmarks, range/gather
  decoding, and skipping. Test indexed gather with multiple blocks/superblocks
  and partial requests. Whole-superblock decoding is currently acceptable.
- Prefer direct mathematical invariants and internal codec tests alongside
  serialization round trips. Use synthetic fixtures and opcode profiles to
  prove that the intended codec was actually exercised.
- Test relevant compression levels explicitly instead of relying exclusively on
  the default. Keep exhaustive level sweeps configurable through pedantic mode;
  bound combination counts so the suite remains practical. Do not multiply
  bookmark cases unless the test already covers them or the change needs them.
- Keep bvector serialization tests under `-bvser`; use independent sparse-vector
  serialization sections. Verify dispatch and harness coverage for new tests.
- Preserve useful implementation variation between `t.cpp` and `t64.cpp` while
  covering equivalent behavior. These test programs remain single-threaded;
  the harness may schedule independent processes in parallel.
- Floating-point tests should retain value comparisons and also check exact
  bit patterns where serialization requires bitwise preservation.
- Verify allocator propagation, workspace reuse after warm-up, and state reset
  across alternating inputs. Use guard checks and sanitizers where appropriate.
- Run focused checks appropriate to the change and state their scope accurately.
  Do not describe a smoke test as a completed full regression or benchmark run.

## Working-tree discipline

- Preserve unrelated local changes. Do not reset, stash, or commit without
  authorization. Keep edits narrowly scoped and reviewable.
- Preserve opcode numbering when retiring codecs: leave vacant slots rather
  than compacting the opcode space.
- Preserve existing wire-format and default-policy behavior unless the change
  explicitly calls for modifying them. Compare complete encoded records when
  choosing compression candidates and retain ordinary encoding on ties.
