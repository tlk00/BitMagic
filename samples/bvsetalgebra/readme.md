# BitMagic logical operations tutorial

This tutorial covers various scenarios of how operations can be programmed to achieve best performance and flexibility under different use-cases and scenarios.

- Logical operations with vectors (2-way and 3-way)
- Logical operations between a vector and a list of integers representing a set
- Logical operations between a vector and a serilaized/compressed vector
- interoperability with STL containers
- Fused operations (combinations of different logical operations or pipelines of operations). 


##Use cases covered:

###OR operation (set UNION)

| Operation      | Description           |
| -------------- | --------------------- |
| `A |= B`       | 2-way operation       |
| `A = A | B`    | 3-way operation       |
| `A = A or B`   | 2-way interpreter mode|
| `A = A or {serialized-BLOB)`  | 2-way interpreter mode with serialized vector|
| `A |= {list of integers}`     | C-style array operation       |
| `A |= {STL container}`      | with STL container via iterator       |
| `A = (B or C or D or ...)`   | variable list aggregation (UNION ALL) |


### AND operation (set intersect).

| Operation      | Description           |
| -------------- | --------------------- |
| `A &= B`       | 2-way operation       |
| `A = A & B`    | 3-way operation       |
| `A = A and B`  | 2-way interpreter mode|
| `A = A and {serialized-BLOB)`  | 2-way interpreter mode with serialized vector|
| `A &= {list of integers}`     | C-style array operation       |
| `A &= {STL container}`      | with STL container via iterator |
| `A = (B and C and D and ...)`   | variable list aggregation (N-way intersect) |


### XOR operation (exclusive OR)

| Operation      | Description           |
| -------------- | --------------------- |
| `A ^= B`       | 2-way operation       |
| `A = A ^ B`    | 3-way operation       |
| `A = A xor B`  | 2-way interpreter mode|
| `A = A xor {serialized-BLOB)`  | 2-way interpreter mode with serialized vector|
| `A ^= {STL container}`      | with STL container via iterator |
|

### SUB operation (set minus / AND NOT).

| Operation      | Description           |
| -------------- | --------------------- |
| `A -= B`       | 2-way operation       |
| `A = A - B`    | 3-way operation       |
| `A = A sub B`  | 2-way interpreter mode|
| `A = A sub {serialized-BLOB)`  | 2-way interpreter mode with serialized vector|
| `A -= {list of integers}`     | C-style array operation       |
| `A -= {STL container}`      | with STL container via iterator |
| `A = B sub (C or D or E or ...)`  | variable list aggregation (N-way MINUS) |



### Invert 

| Operation      | Description           |
| -------------- | --------------------- |
| `A = NOT A`    | in-place set inverter |


### AND-OR (fused and-or)

| Operation      | Description           |
| -------------- | --------------------- |
| `A |= B & C`   | 3-way operation       |
| `A = A or (A and B) or (B and C) or …`    | N-way AND-OR aggregation  |


### AND-SUB (fused AND-MINUS)

| Operation      | Description           |
| -------------- | --------------------- |
| `T = (B and C … ) sub (D or T or … )`  | N-way AND-SUB aggregation |


### AND-SUB-OR  (AND-SUB with UNION ALL)

`bm::aggregator<>::pipeline<>` to run multiple AND-SUB operations on a common group of vectors for better CPU cache optimization.


| Operation      | Description           |
| -------------- | --------------------- |
| `T1 = (B and C … ) sub (D or T or … )`   | step 1      |
| `T2 = (D and T) sub (C or Y or …)`    |   step 2|
|...| ....|
| `TN = (D and C) sub (Z or X or …)`    |   step N|

Results of all steps are aggregated:
`T = T1 or T2 or ... TN`


