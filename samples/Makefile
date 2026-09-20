# ============================================================================
#
# BitMagic Library makefile
# (c) 2002-2009 Anatoliy Kuznetsov.
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

# Build only the samples selected by the top-level Makefile.
# Usage: make [all|clean|rebuild]
#        make bvsample01 [bvsample02 ...]
#        make SAMPLES="bvsample01 bvsample02" [all|clean|rebuild]
# Compiler and platform settings are inherited by the sample Makefiles.

SAMPLES = \
    bvsample01 \
    bvsample02 \
    bvsample03 \
    bvsample04 \
    bvsample05 \
    bvsample06 \
    bvsample07 \
    bvsample08 \
    bvsample09 \
    bvsample10 \
    bvsample11 \
    bvsample12 \
    bvsample14 \
    bvsample15 \
    bvsample16 \
    bvsample17 \
    bvsample18 \
    bvsample19 \
    bvsample20 \
    bvsample21 \
    bvsample22 \
    bvsample23 \
    bvsample24 \
    bvsample25 \
    bvsample26 \
    bvsample27 \
    bvsample01_64 \
    bvsetalgebra \
    bv3vlogic \
    xsample01 \
    xsample02 \
    xsample03 \
    xsample04 \
    xsample04a \
    xsample05 \
    xsample06 \
    xsample07 \
    xsample07a \
    xsample08 \
    xsample09 \
    xsample10 \
    xsample11 \
    strsvsample01 \
    strsvsample02 \
    strsvsample03 \
    strsvsample04 \
    strsvsample05 \
    strsvsample06 \
    strsvsample07 \
    strsvsample10 \
    svfsample01 \
    svfsample02 \
    svfsample03 \
    svfsample04 \
    svfsample05 \
    svsample01 \
    svsample02 \
    svsample03 \
    svsample04 \
    svsample05 \
    svsample06 \
    svsample07 \
    svsample07a \
    svsample08 \
    svsample09 \
    svsample10 \
    rscsample01 \
    rscsample02 \
    rscsample03 \
    rscsample04 \
    rscsample05 \
    rscsample06 \
    rscsample07

.DEFAULT_GOAL := all
.PHONY: all clean rebuild $(SAMPLES)

all: $(SAMPLES)

$(SAMPLES):
	$(MAKE) -C $@ all

clean rebuild:
	@set -e; for dir in $(SAMPLES); do \
		$(MAKE) -C $$dir $@; \
	done
