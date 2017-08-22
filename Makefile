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

TARGETS_BLD =  \
              samples/sample1 \
              samples/sample2 \
              samples/sample3 \
              samples/sample4 \
              samples/sample5 \
              samples/sample6 \
              samples/sample7 \
              samples/sample8 \
              samples/sample9 \
              samples/sample10 \
              tests/stress tests/perf 
        
SHELL=/bin/sh
MAKE := ${MAKE} --no-print-directory

ifeq (${MAKECMDGOALS},)
	MAKECMDGOALS=all
endif

.DIRS_BLD: 
	@for dir in $(TARGETS_BLD); \
	do \
        echo $$dir; \
		${MAKE} -C $$dir $(MAKECMDGOALS) ; \
		if [ "$$?" != "0" ]; then exit 127; fi; \
	done; \
    

rebuild all clean: .DIRS_BLD

relprep: .DIRS_BLD
	$(RM) -r CVS *.plg *.ncb *.opt
	$(RM) -rf samples/CVS tests/CVS platforms/CVS src/CVS
	$(RM) -rf html
	$(RM) -rf src/bm__* src/*vcproj.*
	$(RM) -rf *.suo
	$(RM) -rf debug release cvsenv.sh
	cd platforms; dos2unix *.mk; cd -
	dos2unix readme *.txt
	cd src; dos2unix *.h
	chmod -x src/*
	chmod -x Makefile makefile.in Doxyfile license.txt *.sln readme

install : .DIRS_INST
