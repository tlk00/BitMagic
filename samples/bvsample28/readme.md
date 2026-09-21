# Bit-vector XOR complexity and reference selection

This sample demonstrates `bmcomplexity.h`: exact structural statistics for one
`bm::bvector<>` and for the XOR product of two vectors without constructing the
XOR result. It builds eight small, deterministic vectors, optimizes their
physical representation, measures every standalone vector, compares all 28
unordered pairs, selects a lead vector for each pair and reports the best peer
for every vector.

The example is intentionally chatty. Its output explains which phase is
running and prints the intermediate statistics behind every decision.

## What problem does the sample solve?

XOR reference compression stores one vector (the **lead**) and represents a
related vector by its XOR difference from that lead. If vectors are similar,
the XOR residual is often structurally simpler than either input. Complements
can also produce a simple residual: an all-one range is complex by population
count but simple by run structure.

Population count alone does not express this property. For example, a dense
contiguous interval and a fragmented set can have the same count but very
different compression behavior. `bm::calc_complexity_xor()` supplements XOR
population with run and block-structure information.

## Step 1: construct a small collection

The collection contains no more than ten vectors and covers intentionally
different relationships:

- `A-base` contains intervals, a cross-block interval, one full block and a
  sparse regular region.
- `A-copy` is identical to `A-base`; its XOR is empty.
- `A-near` changes only a few parts of `A-base`.
- `A-complement` is the complement of `A-base` inside the eight-block demo
  universe. It does not invert the complete BitMagic address space.
- `B-base` forms a second pattern family.
- `B-near` is a small variation of `B-base`.
- `sparse` is an unrelated low-population vector.
- `fragmented` sets every third bit and is intentionally difficult for
  run-length compression.

Every vector is passed through `optimize()` before measurement. Most fields in
`bvector_complexity_statistics` describe logical content, but
`bit_to_gap_blocks` also depends on whether either XOR input is physically a
BIT block. Normalizing all inputs makes that observation reproducible.

## Step 2: measure standalone vectors

`bm::calc_complexity(bv)` returns:

- `count`: population count;
- `runs`: maximal runs of one bits across block boundaries;
- `block_runs`: one-runs counted independently inside every block;
- `zero_blocks`, `full_blocks`, `gap_blocks`, `bit_blocks`: content-based block
  classification;
- `gap_words`: total used 16-bit words, including one header per predicted GAP
  block, excluding allocation capacity slack;
- `zero_spans`, `full_spans`: consecutive spans of uniform blocks.

The sample also calls `bvector::calc_stat()` and prints `memory_used`. This is
an important distinction: `calc_stat()` describes the current allocation,
whereas `calc_complexity()` describes content and predicts whether mixed
blocks fit the default GAP policy. Neither the sample score nor the complexity
statistics promise an exact serialized byte count.

## Step 3: define an explicit demonstration cost

To rank candidates, the sample converts a complexity structure to one scalar:

```text
cost(s) = 8192 * bit_blocks
        +    2 * gap_words
        +        zero_spans
        +        full_spans
```

A dense BIT block has an 8192-byte payload. A GAP word occupies two bytes
(`sizeof(bm::gap_word_t)` in the code), so predicted GAP blocks contribute their
complete used payload, including their headers. An empty vector has cost zero.
Uniform spans retain a heuristic one-unit overhead per span.

For each mixed block classified as GAP, the implementation accumulates:

```text
gap_words += number_of_physical_zero_and_one_runs + 1
```

Every alternating run contributes an endpoint word; the additional word is
the GAP header. GAP-to-GAP XOR exposes this length directly. BIT-to-BIT and
cross-type XOR kernels also provide the run count needed to calculate it,
without constructing a GAP result. The same content-based calculation applies
to standalone vectors, keeping their costs comparable with XOR costs.

Used length differs from allocated capacity: `calc_stat()` accounts for GAP
capacity in memory usage, while its serialization allowance uses used length
times two plus three bytes per GAP block. This sample uses the GAP payload
length only; it does not include that serialization allowance or other record
headers. The result is therefore a payload-based heuristic, not a complete
serialized-size bound.

The previous formula charged `4 * block_runs` even inside fixed-size BIT
blocks. That double-counted their representation: eight fragmented BIT blocks
already occupy 65,536 payload bytes regardless of their internal run count.
Now a residual that still contains the same eight BIT blocks receives no
artificial saving just because it has fewer runs. Similarly, full blocks are
charged through their spans, not through a per-block run penalty.

This is a **heuristic coefficient**, not a mathematical metric and not a wire
format model. The weights live in the sample's `complexity_cost()` function so
an application can replace them with workload-specific weights, measured
serialized sizes or a learned predictor. The exact statistics remain useful
independently of this illustrative formula.

Implicit zero blocks cover most of BitMagic's address universe. They are not
added individually to the cost, because doing so would swamp the meaningful
structure of a small vector. Only the number of zero spans receives a small
weight.

## Step 4: compute pairwise XOR statistics

For every unordered pair `(A, B)`, the sample calls:

```cpp
bm::bvector_complexity_statistics xor_stats =
    bm::calc_complexity_xor(A, B);
```

The function does not create a result `bvector`. It uses specialized kernels
where available and small reusable temporary block storage for remaining
cross-representation cases.

XOR is symmetric:

```text
A XOR B == B XOR A
```

Consequently, all XOR complexity fields, the residual cost and the compatibility
score are symmetric. The sample calculates only the 28 unordered pairs rather
than all 56 ordered pairs.

## Step 5: map compatibility to 0 through 100

The cheaper operand becomes the lead, so the XOR residual replaces the more
expensive standalone operand. The score measures that predicted replacement:

```text
baseline = max(cost(A), cost(B))

score = 100                                      if XOR cost is zero
        0                                        if XOR cost >= baseline
        max(1, floor(100 * (baseline-XOR cost)/baseline)) otherwise
```

Interpretation:

- `100`: identical vectors; the residual is empty;
- a high value: the residual is much simpler than the vector it replaces;
- `0`: the residual does not reduce the pair-storage cost under this policy.

The lower bound reserves zero for no predicted improvement, while `floor`
reserves 100 for an exactly empty residual. Thus every positive non-perfect
saving is reported in the range 1 through 99.

Zero means "not a useful XOR reference according to this heuristic," not that
the vectors are incompatible in an algebraic sense. A bounded complementary
pair can score well because its XOR is one long, inexpensive full range.

## Step 6: select the lead vector

Although the XOR product is symmetric, the storage plan still needs a lead.
For the pair `(A, B)` there are two equivalent reconstructions:

```text
store A and (A XOR B), then reconstruct B
store B and (A XOR B), then reconstruct A
```

The residual is identical in both plans. The standalone lead cost is not. This
sample chooses the lower-cost operand as the lead and compares:

```text
independent pair cost = cost(A) + cost(B)
referenced pair cost  = min(cost(A), cost(B)) + cost(A XOR B)
```

`XOR improves` is printed only when the referenced pair cost is lower. This is
why lead selection matters even though XOR itself is symmetric.

In a collection-wide serializer, local pair selection is only the first step.
A production planner may also limit reference-chain depth, account for decode
frequency, require random-access roots, or optimize a global reference graph.

## Step 7: find matching peers

After the complete pair table, the sample chooses the highest symmetric score
for every vector. If scores tie, it prefers the lower referenced pair cost.
The expected families (`A-*` and `B-*`) should find each other, while the
fragmented vector should generally be a poor peer.
If every candidate scores zero, the sample prints `no improving peer`.

## Step 8: verify selected candidates with actual serialization

Output section 4 takes the positive-score best-peer pairs from section 3,
deduplicating pairs selected by both endpoints. It serializes each participating
standalone vector once, then materializes each selected XOR, calls `optimize()`,
and serializes that residual. All measurements use the same serializer with
default settings and report the actual buffer length, including BLOB headers.
Even an empty residual has a serialized representation; its measured size need
not equal the heuristic cost of zero.

Each pair occupies one row emphasizing measured percentage savings, alongside
the heuristic score, before/after byte totals, and the measured lead:

```text
independent bytes = serialized_size(A) + serialized_size(B)
referenced bytes  = serialized_size(lead) + serialized_size(A XOR B)
gain bytes        = independent bytes - referenced bytes
gain percent      = 100 * gain bytes / independent bytes
```

The comparison uses the cheaper serialized operand as lead, which may differ
from the heuristic's predicted lead. Positive savings indicate gain; negative
savings indicate expansion, where independent storage should be retained.
Zero means no gain. This percentage describes the whole
pair, whereas the heuristic score describes replacement of the non-lead vector.

A separate `NOT recommended` table measures all zero-score pairs as negative
controls, using exactly the same procedure. These pairs are selected by their
score, not by measured outcome: any positive saving in this table reveals a
missed opportunity for the heuristic rather than being hidden. Even with the
best measured lead, a negative percentage demonstrates a loss from using XOR.

These are ordinary serialized BLOBs for a lead and a residual, not a persisted
collection with reference metadata. Reference IDs, framing and dependency
management would require additional bytes in a real collection format. The
sample does not write files or build a global reference graph. Measurements
for different pairs must not be summed as collection savings because pairs can
share vectors. Positive-score pairs not chosen as best peers are not measured;
this remains a demonstration rather than an exhaustive serialization planner.

## Build and run

From this directory:

```sh
make
./sample28
```

Or build the `bvsample28` target with the repository CMake configuration.

## Expected output and interpretation

The following output was captured from the supplied macOS run using default
serializer settings and 32-bit addressing. Build/link commands and the shell
prompt are omitted. Memory usage can vary by platform and allocation strategy;
serialized sizes can change with serializer settings or library revisions.
This is an illustrative reference run, not a byte-for-byte regression fixture.

```text
BitMagic XOR complexity and reference-selection demo
Demo universe: [0, 524287] (8 blocks)
Building related, complementary and unrelated vectors...
```

### Output section 1

Section 1 separates population, global one-runs, block-local one-runs, predicted
block types and used GAP words. For A-base, cost is
`8192 + 2 * 12 + 3 + 1 = 8220`. Its 11,088-byte memory figure includes actual
allocation overhead and capacity slack. The fragmented vector pays for eight
BIT payloads, without an additional run penalty.

The demo data occupies eight blocks, but complexity statistics cover the common
BitMagic universe, including implicit trailing zeros. Thus 65,531 zero blocks
for A-base is expected in this 32-bit-address run. The difference between its
704 global runs and 705 block-local runs demonstrates cross-block merging.

```text
1. Measuring optimized standalone vectors
   cost is the sample's structural heuristic; memory is the actual in-memory allocation.

vector              count    runs blk-runs   zero   full    GAP    BIT GAP-words   zspan   fspan        cost     memory
A-base              88409     704      705  65531      1      3      1        12       3       1        8220       11088
A-copy              88409     704      705  65531      1      3      1        12       3       1        8220       11088
A-near              88500     707      708  65529      1      5      1        22       2       1        8239       11600
A-complement       435879     705      708  65529      3      3      1        12       2       3        8221       11088
B-base              69603    1203     1204  65531      0      4      1        14       2       0        8222       11344
B-near              70203    1205     1206  65530      0      5      1        20       3       0        8235       11600
sparse                 80      80       80  65528      0      8      0       176       1       0         353        4176
fragmented         174763  174763   174763  65528      0      0      8         0       1       0       65537       67664
```

### Output section 2

Section 2 evaluates all 28 unordered pairs. Identical vectors give an empty
residual and score 100. A-base and its bounded complement differ in every demo
bit, yet their residual has only one global run across eight full blocks. Its
cost of 2 accounts for one full span and one trailing zero span.

Related pairs score 99, whereas unrelated pairs score zero. These are heuristic
predictions, not measured savings percentages. BIT->GAP is zero throughout this
particular collection: it does not demonstrate a nonempty BIT-to-GAP residual.

```text
2. Computing every unordered XOR pair
   XOR statistics and score are symmetric. The lead is the cheaper standalone vector.

pair                          xor-count     runs  blk-runs    GAP    BIT GAP-words BIT->GAP   xor-cost   score  lead / decision
A-base / A-copy                       0        0         0      0      0         0        0          0    100%  A-base / XOR improves
A-base / A-near                     253        3         3      3      0        12        0         27     99%  A-base / XOR improves
A-base / A-complement            524288        1         8      0      0         0        0          2     99%  A-base / XOR improves
A-base / B-base                  158012     1907      1909      4      2        20        0      16427      0%  A-base / keep independent
A-base / B-near                  158612     1909      1911      5      2        26        0      16438      0%  A-base / keep independent
A-base / sparse                   88461      782       783      7      1       158        0       8509      0%  sparse / keep independent
A-base / fragmented              204232   174528    174529      0      8         0        0      65537      0%  A-base / keep independent
A-copy / A-near                     253        3         3      3      0        12        0         27     99%  A-copy / XOR improves
A-copy / A-complement            524288        1         8      0      0         0        0          2     99%  A-copy / XOR improves
A-copy / B-base                  158012     1907      1909      4      2        20        0      16427      0%  A-copy / keep independent
A-copy / B-near                  158612     1909      1911      5      2        26        0      16438      0%  A-copy / keep independent
A-copy / sparse                   88461      782       783      7      1       158        0       8509      0%  sparse / keep independent
A-copy / fragmented              204232   174528    174529      0      8         0        0      65537      0%  A-copy / keep independent
A-near / A-complement            524035        4        11      3      0        12        0         27     99%  A-complement / XOR improves
A-near / B-base                  158103     1910      1912      4      2        24        0      16435      0%  B-base / keep independent
A-near / B-near                  158703     1912      1914      5      2        30        0      16446      0%  B-near / keep independent
A-near / sparse                   88552      785       786      7      1       164        0       8521      0%  sparse / keep independent
A-near / fragmented              204263   174526    174527      0      8         0        0      65537      0%  A-near / keep independent
A-complement / B-base            366276     1908      1909      4      2        20        0      16427      0%  A-complement / keep independent
A-complement / B-near            365676     1910      1911      5      2        26        0      16438      0%  A-complement / keep independent
A-complement / sparse            435827      783       786      7      1       158        0       8509      0%  sparse / keep independent
A-complement / fragmented        320056   174528    174530      0      8         0        0      65537      0%  A-complement / keep independent
B-base / B-near                    1002        2         2      2      0         8        0         18     99%  B-base / XOR improves
B-base / sparse                   69661     1283      1284      7      1       174        0       8541      0%  sparse / keep independent
B-base / fragmented              197964   174364    174366      0      8         0        0      65537      0%  B-base / keep independent
B-near / sparse                   70261     1285      1286      7      1       178        0       8549      0%  sparse / keep independent
B-near / fragmented              198164   174362    174364      0      8         0        0      65537      0%  B-near / keep independent
sparse / fragmented              174789   174736    174736      0      8         0        0      65537      0%  sparse / keep independent
```

### Output section 3

Section 3 chooses the best peer per vector and the cheaper heuristic lead.
Pair costs here are heuristic units, not serialized bytes. A-base/A-copy and
B-base/B-near appear from both endpoints; these are the same unordered pairs,
not additional independent savings. Sparse and fragmented have no improving
peer under the predictor. This list is candidate selection, not a global
reference graph or collection storage plan.

```text
3. Best peer for each vector
   This ranks the symmetric scores; ties prefer the lower referenced pair cost.

   A-base          -> A-copy          score=100%  lead=A-base          pair-cost 16440 -> 8220
   A-copy          -> A-base          score=100%  lead=A-base          pair-cost 16440 -> 8220
   A-near          -> A-base          score= 99%  lead=A-base          pair-cost 16459 -> 8247
   A-complement    -> A-base          score= 99%  lead=A-base          pair-cost 16441 -> 8222
   B-base          -> B-near          score= 99%  lead=B-base          pair-cost 16457 -> 8240
   B-near          -> B-base          score= 99%  lead=B-base          pair-cost 16457 -> 8240
   sparse -> no improving peer
   fragmented -> no improving peer
```

### Output section 4

Section 4 validates four unique recommended pairs and all 21 zero-score pairs.
The recommended pairs save 32.63% to 46.43% of total pair bytes in this run.
Every negative control loses space, even with the cheaper measured lead.
No missed opportunity appears among the zero-score pairs in this collection;
this does not establish the predictor's accuracy on other workloads.

The score and pair savings have different denominators. For A-base/A-near,
score 99 describes predicted replacement of the non-lead vector. Measured
pair savings are `(95 - 64) / 95 = 32.63%`, including the retained lead.
An identical pair saves less than 50% because the empty XOR still has a
three-byte BLOB: `84 -> 45` saves 46.43%.

Negative savings can exceed 100% in magnitude. B-base/fragmented grows from
152 to 2,840 bytes, approximately 18.68 times the independent size; the extra
2,688 bytes are a 1768.42% expansion. XOR can destroy regular patterns that the
serializer compresses efficiently. These measured losses justify retaining
independent storage.

The capture's `scorepair savings` heading denotes two adjacent columns:
heuristic score and measured pair savings. All values below are preserved
from the supplied run.

```text
4. Actual serialization of selected pairs (unique pairs only)
   Default serializer settings; sizes in bytes. Materializing and optimizing each XOR.
   Totals include BLOB headers, but exclude reference IDs and collection framing.

   --- Recommended best-peer pairs ---
pair                           scorepair savings   independent -> lead+XOR bytes / result
A-base / A-copy                  100      46.43%        84 ->       45  gain; lead=A-base
A-base / A-near                   99      32.63%        95 ->       64  gain; lead=A-base
A-base / A-complement             99      44.71%        85 ->       47  gain; lead=A-base
B-base / B-near                   99      35.58%       104 ->       67  gain; lead=B-base

   --- NOT recommended: all zero-score pairs ---
   Negative savings mean expansion. A gain here is a missed opportunity for the heuristic.
pair                           scorepair savings   independent -> lead+XOR bytes / result
A-base / B-base                    0     -29.89%        87 ->      113  loss; lead=A-base
A-base / B-near                    0     -23.76%       101 ->      125  loss; lead=A-base
A-base / sparse                    0    -774.14%        58 ->      507  loss; lead=sparse
A-base / fragmented                0   -1178.52%       149 ->     1905  loss; lead=A-base
A-copy / B-base                    0     -29.89%        87 ->      113  loss; lead=A-copy
A-copy / B-near                    0     -23.76%       101 ->      125  loss; lead=A-copy
A-copy / sparse                    0    -774.14%        58 ->      507  loss; lead=sparse
A-copy / fragmented                0   -1178.52%       149 ->     1905  loss; lead=A-copy
A-near / B-base                    0    -153.06%        98 ->      248  loss; lead=B-base
A-near / B-near                    0    -139.29%       112 ->      268  loss; lead=A-near
A-near / sparse                    0    -646.38%        69 ->      515  loss; lead=sparse
A-near / fragmented                0   -1142.50%       160 ->     1988  loss; lead=A-near
A-complement / B-base              0     -29.55%        88 ->      114  loss; lead=A-complement
A-complement / B-near              0     -23.53%       102 ->      126  loss; lead=A-complement
A-complement / sparse              0    -759.32%        59 ->      507  loss; lead=sparse
A-complement / fragmented          0   -1170.67%       150 ->     1906  loss; lead=A-complement
B-base / sparse                    0    -521.31%        61 ->      379  loss; lead=sparse
B-base / fragmented                0   -1768.42%       152 ->     2840  loss; lead=B-base
B-near / sparse                    0    -416.00%        75 ->      387  loss; lead=sparse
B-near / fragmented                0   -1667.47%       166 ->     2934  loss; lead=B-near
sparse / fragmented                0    -305.69%       123 ->      499  loss; lead=sparse
```

## Questions and answers

### Is the score a metric?

No. The weighted coefficient is a ranking heuristic. It is clipped to zero,
and complementary vectors can have a simple XOR despite being far apart by
Hamming distance. Do not assume triangle inequality or metric-space behavior.

### Is the score symmetric?

Yes. Its baseline uses `max(cost(A), cost(B))`, and the XOR residual is
symmetric. Lead selection uses the corresponding `min`: retain the cheaper
operand and replace the more expensive operand with the residual.

### Does a high score guarantee a smaller serialized BLOB?

No. Serialization has headers, encoding choices, bookmarks and potentially
application-specific XOR machinery. Use actual serialization when an exact
decision is required. The coefficient is intended as a cheap screening signal
for a large candidate set.

### Why keep population count if the score is run-oriented?

Population is still valuable for diagnostics and for policies which consider
sparse-list or BIC-style representations. This simple formula deliberately
focuses on BIT/GAP/run structure; applications can incorporate `count`.

### Why report both `runs` and `block_runs`?

`runs` is the exact number of global one-runs. `block_runs` treats every block
independently and therefore retains block-local encoding pressure. A run ending
with one at a block boundary and continuing with one in the next block is one
global run but contributes to both blocks locally.

### Why can `bit_to_gap_blocks` change after optimization?

It records XOR residual blocks predicted to fit GAP representation when at
least one input was physically a BIT block. It is therefore partly an input
storage diagnostic. Optimize vectors consistently before comparing it.
