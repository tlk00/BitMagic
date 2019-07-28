## BitMagic C++ Library


Algorithms and tools for Algebra od Sets for information retrieval, 
indexing of databases, scientific algorithms, ranking, clustering, unsupervised machine learning 
and signal processing.

BitMagic library uses compressed bit-vectors as a main vehicle for implementing set algebraic operations, 
because of high efficiency and bit-level parallelism friendly for parallel processing with SIMD. 


### Main Features (bm::bvector<>):

- compressed bit-vector container 
- iterator (bm::bvector<>::enumerator to decode the bitset to integers
- set algebraic operations: AND, OR, XOR, MINUS, NOT on bit-vectors and integer sets
- fast bit-vector ierator (enumerator) for bit-vector traversal, algorithms for functor-based
  traversals (similar to for_each)
- fast import (compression) of integer lists (C++ style bulk_insert_iterator or using C array)
- aggregator: fast logical AND, OR, AND-MINUS operations on groups of bit-vectors
- serialization/hybernation of bit-vector containers into compressed BLOBs for persistence (or in-RAM compression)
- set algebraic operations on compressed BLOBs (on the fly deserialization with set-algebraic function)
- statistical algorithms to efficiently construct similarity and distance metrics, measure similarity between bit-vectors, 
integer sets and compressed BLOBs
- operations with rank: population count distances on bit-vector. Rank-Select operations are often 
used in succinct data structures

### Serialization with compression

BitMagic can serialize-compress bit-vector for efficient storage and transfer.
Naturally BitMagic uses integer compression schemes: GAP-RLE-Delta combined with Elias Gamma coding or 
Interpolated Binary Coding. 
 
BitMagic serialization-compression is a hybrid model, where bit-vectors are partitioned in order to adaptively set the best 
representation. BitMagic compression is adaptive and it remains efficient for both sparse and dense 
distributions of bits. Dense blocks are coded using complementary compressed integer lists.

BitMagic is tested on Gov2 benchmark set of inverted lists.
[http://bitmagic.io/bm5-cmpr.html](http://bitmagic.io/bm5-cmpr.html)


### Succinct vectors 

BitMagic supports succinct (memory compact) vectors based on bit-transposition transform.
Bit transposition solves two purposes: free unused bit plains and isolate regularity and enthrothy into
separately compressed bit-vectors. Compression on bit-planes offers both superior memory performance and fast search.

Succinct bit-transposed implementation works well for both integer vectors and string vectors and 
it rivals other succinct schemes like prefix trees. Succinct vectors can be both sorted and unsorted.   

#### Main features

- sparse vector(s) for native int types using bit transposition and separate compression of bit-plains, 
with support of NULL values (unassigned) for construction of in-memory columnar structures. Bit-transposed
sparse vectors can be used for on-the fly compression of astronomical, molecular biology or other data,
efficient store of associations for graphs, etc.
- search algorithms for sorted and unsorted succinct vectors (vectors of ints or strings)
- algorithms on sparse vectors: dynamic range clipping, search, group theory image (re-mapping).
- all containers are serializable with compression

### Memory profiling/monitoring

BitMagic implements memory profiling calls for all vectors. Any vector can be sampled for memory footprint so the 
top level system can adapt memory managemet based on exact footprint your vectors. 
BitMagic memory managemnt is industry grade and ready for integration into large scale systems 
(Now try to evaluate the footprint of an STL map<>, for example).


### 64-bit

Yes!
BitMagic supports 64-bit, can be used with 32-bit address space (less overhead) or 64-bit address space.
32-bit address space is the default mode 2^31-1 elements should be a good fit for short to medium range
search systems. 64-bit address mode is available using #define BM64ADDR or #include "bm64.h".
Current 64-bit implementation now allows 2^48-1 elements for large scale systems.

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
with a set of examples. 

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
As a low level library BitMagic needs to be stable and conformant to be useful.

We do not rely on unit tests alone, our tests often use _"chaos testing"_ where stress tests
are based on randomized, generated sets and randomized operations. We regularly build
and run test suits for Release and Debug mode for various combinations of SIMD 
optimizations. 

All variants of test builds take days to run, so the working master branch is 
not guaranteed to be perfect.

All of the above does not guarantee it is bug free project.


### HPC and SIMD:

BitMagic provides fast algorithms optimized for various SIMD sets, exploit super-scalar 
features of modern CPUs, cache-friendly algorithms, thread-safe and data-parallel structures. 

- SSE2    - Intel SSE2 but gets gradually phased out, no new algorithms and compute kernels
- SSE4.2  - Supported and adds x86 POPCNT  
- AVX2    - Supported adds POPCNT and BMI1/BMI2 (best performance!)
- AVX-512 - Work in progress. Some algorithms are ready(and stable), but we cannot get performance right (yet).
            AVX-512 projects using BitMagic for now should just use AVX2 mode (it is compatible).


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

	$ . ./bmenv.sh
	
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
MSVC - solution and projects are available via CMAKE generation

###### MacOS
---

XCODE - project files are available via CMAKE generation

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

Yes, we need it! 
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

