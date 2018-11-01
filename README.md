## BitMagic Library


Algorithms and tools for integer set algebra operations used for information retrieval, 
indexing of databases, scientific algorithms, ranking, clustering and signal processing.

BitMagic library uses compressed bit-vectors as a main vehicle for implementing set algebraic operations, 
because of high efficiency and bit-level parallelism of this representation. 


### Main Features:

- compressed bit-vector container with mechanisms to iterate integer set it represents
- set algebraic operations: AND, OR, XOR, MINUS on bit-vectors and integer sets
- fast bit-vctor ierator (enumerator) for bit-vector traversal
- fast import (compression) of integer lists (C++ style bulk_insert_iterator or using C array)
- aggregator: fast logical AND, OR, AND-MINUS operations on groups of bit-vectors
- serialization/hybernation of bit-vector containers into compressed BLOBs for persistence (or in-RAM compression)
- set algebraic operations on compressed BLOBs (on the fly deserialization with set-algebraic function)
- statistical algorithms to efficiently construct similarity and distance metrics, measure similarity between bit-vectors, 
integer sets and compressed BLOBs
- operations with rank: population count distances on bit-vector
- sparse vector(s) for native int types using bit transposition and separate compression of bit-plains, 
with support of NULL values (unassigned) for construction of in-memory columnar structures. Bit-transposed
sparse vectors can be used for on-the fly compression of astronomical, molecular biology or other data,
efficient store of associations for graphs, etc.
- algorithms on sparse vectors: dynamic range clipping, search, group theory image (re-mapping).
Collection of algorithms is increasing, please check our samples and the API lists.


### Features In Progress:

- compressed binary relational and adjacency matrixes and operations on matrixes for 
Entity-Relationship acceleration, graph operations, social analyticsm materialized RDBMS joins, etc 

- portable C-library layer working as a bridge to high level languages like Python, Java, Scala, .Net 

### License: 

Apache 2.0. 
We ask you to explicitly mention BitMagic project in any work using it or derived from our code 
or published materials.


### SIMD support:

BitMagic provides fast algorithms optimized for various SIMD sets, exploit super-scalar features of modern CPUs, 
cache-friendly algorithms and thread-safe and data-parallel structures.

- SSE2    - Intel SSE2 command set is supported, but gets gradually phased out, no new algorithms and compute kernels
- SSE4.2  - Supported
- AVX2    - Supported and recommended as fastest for high-performance systems
- AVX-512 - Work in progress. Some algorithms are ready(and stable), but we cannot get performance right (yet).
            AVX-512 projects using BitMagic for now should just use AVX2 mode (it is compatible).



### If you want to contribute or support BitMagic library:
---

1. GitHub master accepts patch requests
Our branching policy is that master cannot be considered fully stable between the releases.
(for production stability please use release versions)

2. We need help with real-life cases and benchmarks

3. We need help with mappings to Python and other languages (BitMagic has C bindings)


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

---

API documentation and examples:
[http://www.bitmagic.io/apis.html](http://www.bitmagic.io/apis.html)


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

