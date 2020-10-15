## BitMagic C++ Library

BitMagic was created as a Algebra of Sets toolkit for Information Retrieval but currently evolved into a more general Data Science components library forn memory compact big data structures and algorithms on succinct data vectors. 
BitMagic implements compressed bit-vectors and containers (vectors) based on ideas of bit-plane transform
and Rank-Select compression. All containers are serializable (with compression) for efficient storage and
network transfer. 

BitMagic offers sets of method to architect your applications to use HPC techniques to save 
memory (thus be able to fit more data in one compute unit) and improve storage and traffic patterns 
when storing data vectors and models in files or object stores (SQL or noSQL).

BitMagic facilitates two big classes of scenarios:
- limited RAM applications (WebAssembly, IoT, Edge computing, desktop/workstation apps)
- Big Data HPC (petabyte problems) where amount of data is astronomical and computation 
and storage has to be distributed (need for efficient transfer)

### Applications and Use Cases 

BitMagic was used as a building blocks for:

- Algebra of Sets for IR, database inverted index construction
- Simulation of logical schemes and topologies (like FPGAs)
- Multidimentional binary distances on compressed bit-vectors 
  building blocks for binary clusterizations, Self-organizing maps
- Construction of memory compressed bioinformatics models 
  - sequence alignments 
  - collections of variations and SNPs 
  - compression of sequence reads 
  - k-mer classification systems
- Visualization of data-sets in the memory constrained edge configurations
  (edge comuting, IoT, WebAssembly)
- Graph and Tree analysis, construction of compressive matrixes of association 
- Web and application logs analytics and systems reliability automation 
- Low latency network scheduling (orchestration) system

Please visit our use-case notes:
[http://bitmagic.io/use-case.html](http://bitmagic.io/use-case.html)

### Optimizations and SIMD

BitMagic library is a high performance library, implementing custom optimizations
for variety of platforms and build targets:
- x86 (platform specific available bit-scan instructions)
- x86 SIMD: SSE2, SSE4.2(POPCNT, LZCNT), AVX2 (BMI1/BMI2), AVX-512(work in progress)
- Arm (use of specific available bit-scan instructions)
- WebAssembly (use of WebAsm built-ins and platform specific tricks)

BitMagic uses a data-parallel vectorized design with a goal not just provide a best single 
threaded performance but to facilitate highly parallel compute on many-core systems. 

### Compression algorithms

BitMagic uses a suite of compression algorithms, filters and transformations to reduce 
memory footprint, storage costs and network data transfer.
[http://bitmagic.io/design.html](http://bitmagic.io/design.html)

- Hierachical compression of bit-vectors
- D-GAP (RLE) of bit-blocks
- Binary Inetrpolative Coding (BIC)
- Elias-Gamma Coding 
- Bitwise-Tranposition for vectors (also known as Bit Planes coding or bit-slicing)
- XOR compression filters
- Frequency based dictionary remapping (similar to Huffman codes)  


Please visit our tech.notes:
[http://bitmagic.io/articles.html](http://bitmagic.io/articles.html)


### Main Features (bm::bvector<>)

- compressed bit-vector container 
- iterator (bm::bvector<>::enumerator to decode the bitset to integers
- set algebraic operations: AND, OR, XOR, MINUS, NOT on bit-vectors and integer sets
- fast bit-vector ierator (enumerator) for bit-vector traversal, algorithms for functor-based traversals (similar to std::for_each)
- fast import (compression) of integer lists (C++ style bulk_insert_iterator or using C array)
- aggregator: fast vectorized logical AND, OR, AND-MINUS operations on groups of bit-vectors
- serialization/hybernation of bit-vector containers into compressed BLOBs for persistence (or in-RAM compression)
- set algebraic operations on compressed BLOBs (on the fly deserialization with set-algebraic function)
- statistical algorithms to efficiently construct similarity and distance metrics, measure similarity between bit-vectors, 
integer sets and compressed BLOBs
- operations with rank: population count distances on bit-vector. 
Rank-Select operations are often used in succinct data structures, BitMagic implements a 
compact RS Index accelerated with SIMD and BMI (PDEP)

### Ranges and Intervals (bmintervals.h)

BitMagic supports re-interpretation of bit-vectors as collections of non-overlapping ranges 
of 1s flanked with 0s (for example: 011110110). Regular set functions provide set intersect / unions
intreval operations implement interval iterator and searches for interval boundaries.
Ranges and intervals has great utility in bioinformatics, because genomics data are often annotated as 
coordinate ranges. BitMagic offers building blocks for effciient oprations on intervals encoded 
as bit-vectors (find interval start/end, check if range is an inetrval, iterate intervals

### Serialization with compression

BitMagic uses contept of two-stage serialization-deserialization. 
The focus is on fast deserialization. BitMagic implements API for fast vector range deserialization 
and gather deserialization of compressed BLOBs. The ultimate feature of BitMagic is ability to work
with compressed data.

#### Stage One: succinct memory 
This is the main in-RAM operational state, where vectors are kept in memory compact form.
Succinct is NOT a compression. It is possible to access random elements in containers,
decode blocks, iterate vectors, make updates, run search algorithms. Stage One offers 
transparent use, it vectors look much like STL. Succinct is memory compact but not 
fully compressed. 

#### Stage Two: compression
BitMagic can serialize all containers and vectors with additional compression based on block of heuristics 
and codecs. The workhorse coding techniques are: Binary Interpolative Coding (BIC) and Elias Gamma.

BitMagic containers are called "sparse" vectors but in fact its compression schemes works well for both 
sparse and dense data.

BitMagic is tested on Gov2 benchmark set of inverted lists and number of proprietory data sets.
[http://bitmagic.io/bm5-cmpr.html](http://bitmagic.io/bm5-cmpr.html)

#### Decompression 
Deserialization always go back to Stage One, so data are not completely decoded but instead  
succinct in RAM. Our goal here is to both reduce application memory footprint and improve deserialization 
latency. Decompression algorithms support deserialization of arbitrary ranges or even gather 
deserializatin of elements.


### Succinct vectors 
BitMagic supports succinct (memory compact) vectors based on bit-transposition transform  
also known as bit-plane compression (BPC) or bit-slicing and Rank-Select compression. 
Succint vectors are labeled "sparse" but they also work for dense vectors.

Bit transposition solves two purposes: free unused bit plains and isolate regularity and entropy into
separate (sparse) bit-vectors. Compression on bit-planes offers both superior memory performance and fast search. One of the design golas was to be able to perform searches on succinct data using vectorized
logical operations.

Succinct bit-transposed implementation works well for both integer vectors and string vectors and 
it rivals other succinct schemes like prefix trees. Succinct vectors can be both sorted and unsorted.
The idea here is similar to Apache Arrow-Parquet, but it takes it farther with bit-plane compression 
and extensive use of accelerated Rank-Select compression.

- bit-transposed representation offers best memory footprint for numeric vectors when data 
uses numbers with limited or variable bit rate. If data needs just 27 bits succinct vector will use 
just that and not the nearest natural 32-bit type. It is adaptive and completely automatic.
- If string vector needs just a few bits to represent a char in particular position 
(DNA strings, chemical compounds as smiles, etc. ) BitMagic has an option to do transparent remapping,
so DNA string vector would use just 2-3 bits (for ATGCN alphabet), remapping analyses frequencies and
similar to Huffman. 
- Rank-Select approach allows to collapse all NULL values and save RAM 
- For deeper compression BitMagic also implements XOR filter which finds possible correlations 
between bit-planes to reduce enthropy 
[http://bitmagic.io/bm-xor.html](http://bitmagic.io/bm-xor.html)



#### Main features
- sparse vector(s) for native int types using bit transposition and separate compression of bit-plains, 
with support of NULL values (unassigned) for construction of in-memory columnar structures. Bit-transposed
sparse vectors can be used for on-the fly compression of astronomical, molecular biology or other data,
efficient store of associations for graphs, etc.
- search algorithms for sorted and unsorted succinct vectors (vectors of ints or strings)
- algorithms on sparse vectors: dynamic range clipping, search, group theory image (re-mapping).
- all containers are serializable with compression (binary interpolative coding, elias gamma coding)

#### Serialization and versioning
BitMagic supports serialization evolution - if serialization format changes, 
old saved data remains readable by the new code. Old code will NOT be able to read new BLOBs.
BitMagic changes major version number when serialization format changes.


### Memory profiling/monitoring
BitMagic implements memory profiling calls for all vectors. Any vector can be sampled for memory footprint so the top level system can adapt memory managemnet based on the runtime memory profiling.
Typical use case is memory cache of objects with compression to RAM and then to eviction to disk based
on resource consumption and costs (dynamic balance of demand and supply).


### 64-bit vs 32-bit

Yes!
BitMagic supports 64-bit, can be used with 32-bit address space (less overhead) or full 64-bit address space.
32-bit address space is the default mode 2^31-1 elements should be a good fit for short to medium range
IR and data science systems. 64-bit address mode is available using #define BM64ADDR or #include "bm64.h".
Current 64-bit implementation allows 2^48-1 vector elements for large scale systems.

### WebAssembly

BitMagic compiles and work with WebAssmbly (emscripten). Latest versions includes 
multiple tweaks, specific for the platform. Performance numbers are close to native code 
without SIMD (sometimes afster). Sample compile line would look like:

`emcc -std=c++11 -s ALLOW_MEMORY_GROWTH=1 -O2 -s WASM=1 ... `

### ARM

BitMagic fully supports ARM CPU. All releases are stress tested with Raspberry Pi 4.
BitMagic implements some algorithmic tweaks and improvements specific for ARM 
(like use of LZCNT instruction). BitMagic succinct containers can be very useful on embedded 
systems for edge computing with limited amount of available memory.

ARM NEON SIMD support is a work in progress.


### C-library interface:

- BitMagic library provides C-library wrapper, which builds as a "true C" library.
For redistribution it does NOT require C++ Runtime, because it compiles without use of
STL, C++ memory allocation (operator new) or exceptions. Our goal here is to eventually
provide a bridge to other languiages of data science (Python) and languages of enterprise 
scale development (Java, Scala) via JNI. 

### Features In Progress:

- compressed binary relational and adjacency matrixes and operations on matrixes for 
Entity-Relationship acceleration, graph operations, social analyticsm materialized RDBMS joins, etc 

- succinct data structures and containers based on bit-transposed data representation and 
rank-select compression

- memory efficient dictionaries of strings as an alternative to suffix trees

### How to start with BitMagic?
---

BitMagic C++ is a header only library (easy to build and use in your project) and it comes 
with a set of examples. It is NOT advised to use tests as a code example to study library usage.
Tests do not illustate the best usage patterns and models and often intentionally inefficient.

API documentation and examples:
[http://www.bitmagic.io/apis.html](http://www.bitmagic.io/apis.html)

Tutorial for Algebra of Sets:
[http://bitmagic.io/set-algebra.html](http://bitmagic.io/set-algebra.html)

Use Cases and Application notes:
[http://bitmagic.io/use-case.html](http://bitmagic.io/use-case.html)

Technical notes on performance optmization:
[http://bitmagic.io/articles.html](http://bitmagic.io/articles.html)

Doxygen:
[http://bitmagic.io/doxygen/html/modules.html](http://bitmagic.io/doxygen/html/modules.html)



### License: 

Apache 2.0. 

**Important!** We ask you to explicitly mention BitMagic project in any derived work or our published 
materials. Proper reference on your product/project page is a REQUIREMENT 
for using BitMagic Library.


### Quality Assurance:

BitMagic library pays serious attention to code quality and test coverage.  
As a building blocks library BitMagic needs to be stable and conformant to be useful.

We do not rely on unit tests alone, our tests often use _"chaos testing"_ (aka fuzzing)
where stress tests are based on randomized, generated sets and randomized operations. 
We regularly build and run test suits for Release and Debug mode for various combinations of SIMD 
optimizations. 

All variants of test builds take days to run, so the working master branch is 
NOT guaranteed to be perfect all the time. For production please use stable 
github release branches or distributions from SourceForge:
https://sourceforge.net/projects/bmagic/files/



### If you want to contribute or support BitMagic library:
---

1. GitHub master accepts patch requests
Our branching policy is that master cannot be considered fully stable between the releases.
(for production stability please use release versions)

2. Need help with mappings to Python and other languages (BitMagic has C bindings)


### How to build BitMagic C++ library:
---

BitMagic C++ is a header-only software package and you probably can just take the
sources and put it into your project directly. All C++ library sources/headers are in src
directory.

However if you want to use our makefiles you need to follow the next simple
instructions:


###### Unix:
---

1. Traditional (in-place build)

Apply a few environment variables by runing bmenv.sh in the project root directory:

	$ source ./bmenv.sh
	
(please use "." "./bmenv.sh" to apply root environment variable)

use GNU make (gmake) to build installation.


	$make rebuild


or (DEBUG version)
 
	$gmake DEBUG=YES rebuild

The default compiler on Unix and CygWin is g++.
If you want to change the default you can do that in makefile.in
(should be pretty easy to do)

2. CMake based build
Project now comes with a set of makefiles for cmake, you can just build it or generate project files for any 
cmake-supported environment.


###### Windows:
---

If you use cygwin installation please follow general Unix recommendations.
MSVC - solution and projects are available via CMAKE. 

###### MacOS
---

XCODE - project files are available via CMAKE.

---

BitMagic library for C and JNI mappings.

BitMagic library is available for C language (this is work in progress).
The main objective of C build is to bridge BitMagic into other programming languages.
C build is in the subdirectory "lang-maps".

C build creates versions of BitMagic build for SSE and AVX and adds CPU identification,
so the upper level system can support dynamic CPU identification and code dispatch.

C build uses C++ compiler, but does not use RTTI, exceptions (simulated with long jump)
and C++ memory management, so it is C++ language neutral, without runtime dependencies.
Algorithms and behavior are shared between C and C++.

Current state of development: 
   - bit-vector functionality is available via C interface

#### Python support

Python support is pending and we need help here.
If you are enthusiastic about Python and think you can help please contact:
anatoliy.kuznetsov @ yahoo dot com



Fine tuning and optimizations:
---

All BM fine tuning parameters are controlled by the preprocessor defines (and 
compiler keys). 

---

BM library supports CXX-11. Move semantics, noexept, initalizer lists.

---

To turn on SSE2 optimization #define BMSSE2OPT in your build environment.
To use SSE4.2  #define BMSSE42OPT
SSE42 optimization automatically assumes SSE2 as a subset of SSE4.2. 
(you don’t have to use both BMSSE2OPT and BMSSE42OPT at the same time). 

To turn on AVX2 - #define BMAVX2OPT
This will automatically enable AVX2 256-bit SIMD, popcount (SSE4.2) and other 
compatible hardware instructions.

BM library does NOT support multiple code paths and runtime CPU identification.
You have to build specifically for your target system or use default portable
build.

To correctly build for the target SIMD instruction set - please set correct 
code generation flags for the build environment.

BitMagic examples and tests can be build with GCC using cmd-line settings: 

	make BMOPTFLAGS=-DBMAVX2OPT rebuild
or

	make BMOPTFLAGS=-DBMSSE42OPT rebuild

It automatically applies the right set of compiler (GCC) flags for the target 
build.

	CMAKE
	
	cd build
	cmake -DBMOPTFLAGS:STRING=BMSSE42OPT ..
	make

OR

	cmake -DBMOPTFLAGS:STRING=BMAVX2OPT ..


---

BM library supports "restrict" keyword, some compilers 
(for example Intel C++) generate better
code (out of order load-stores) when restrict keyword is helping. This option is 
turned OFF by default since most of the C++ compilers does not support it. 
To turn it ON please #define BM_HASRESTRICT in your project. Some compilers
use "__restrict" keyword for this purpose. To correct it define BMRESTRICT macro 
to correct keyword. 

---


If you want to use BM library in a "no-STL project" (like embedded systems) 
define BM_NO_STL.

This rule only applies to the core bm::bvector<> methods. 
Auxiliary algorithms, examples, etc would still use STL.

---


Follow us on twitter: [https://twitter.com/bitmagicio](https://twitter.com/bitmagicio)


Thank you for using BitMagic library!

e-mail:   info@bitmagic.io

WEB site: [http://bitmagic.io](http://bitmagic.io)

GitHub:   [https://github.com/tlk00/BitMagic](https://github.com/tlk00/BitMagic)

