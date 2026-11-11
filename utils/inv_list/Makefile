# ============================================================================
#
# BitMagic Library makefile
# (c) 2017 Anatoliy Kuznetsov.
#
# ============================================================================
# Permission is hereby granted, free of charge, to any person 
# obtaining a copy of this software and associated documentation 
# files (the "Software"), to deal in the Software without restriction, 
# including without limitation the rights to use, copy, modify, merge, 
# publish, distribute, sublicense, and/or sell copies of the Software, 
# and to permit persons to whom the Software is furnished to do so, 
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included 
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
# OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
# OTHER DEALINGS IN THE SOFTWARE.
# ============================================================================

PROJECTNAME=BitMagic Library Sparse Vector Utility

CXXSOURCES=$(wildcard *.cpp)

BASICINCS=src

CSOURCES= 

DESTINATION = bminv

DESTTYPE = 
DESTINATIONDIR=utils/inv_list

# SHA-256 is optional. macOS and native Windows use system crypto providers.
INV_UNAME := $(shell uname -s)
ifeq ($(BM_NO_SHA256),YES)
    override BMOPTFLAGS += -DBM_NO_SHA256
else ifneq ($(findstring MINGW,$(INV_UNAME)),)
    EXTERNLIBS += -lbcrypt
else ifneq ($(INV_UNAME),Darwin)
    INV_CRYPTO_LIBS := $(shell pkg-config --libs libcrypto 2>/dev/null)
    ifneq ($(strip $(INV_CRYPTO_LIBS)),)
        EXTERNLIBS += $(INV_CRYPTO_LIBS)
        CXX_INCLUDE_FLAGS += $(shell pkg-config --cflags libcrypto 2>/dev/null)
    else
        override BMOPTFLAGS += -DBM_NO_SHA256
        $(warning libcrypto not found via pkg-config; building exact-only inv_list)
    endif
endif

include ../../makefile.in

# Shared serializer settings must rebuild the app even without generated dependencies.
$(CXXOBJECTS): $(wildcard inv_list*.h) ../../src/bmbvhash.h

# Build this utility with BitMagic's native AArch64 NEON backend. Force the
# object rebuild because changing BMOPTFLAGS alone does not update timestamps.
.PHONY: neon
neon:
	@case "$$(uname -m)" in arm64|aarch64) ;; \
	    *) echo "NEON build requires an AArch64 host" >&2; exit 1 ;; esac
	$(MAKE) --no-print-directory -B BMOPTFLAGS=-DBMNEONOPT $(DESTINATION)

# Standalone sparse-vector collection experiment; it is not part of bminv.
.PHONY: svexp-neon
svexp-neon:
	@case "$$(uname -m)" in arm64|aarch64) ;; \
	    *) echo "NEON build requires an AArch64 host" >&2; exit 1 ;; esac
	$(CXX) $(CXXFLAGS) -DBMNEONOPT $(DBGFLAGS) \
	    experiments/sparse_vector/svexp.cpp -o svexp $(LDFLAGS)

# Standalone block-composition scanner for raw Gov2 collections.
.PHONY: gov2-block-stats-neon
gov2-block-stats-neon:
	@case "$$(uname -m)" in arm64|aarch64) ;; \
	    *) echo "NEON build requires an AArch64 host" >&2; exit 1 ;; esac
	$(CXX) $(CXXFLAGS) -DBMNEONOPT $(DBGFLAGS) \
	    experiments/block_stats/gov2_block_stats.cpp -o gov2_block_stats $(LDFLAGS)

# Serialized bvector playlists for opcode survey/replay experiments.
.PHONY: gov2-playlist-neon
gov2-playlist-neon:
	$(CXX) $(CXXFLAGS) -DBMNEONOPT $(DBGFLAGS) \
	    experiments/opcode_playlist/playlist.cpp -o gov2_playlist $(LDFLAGS)
