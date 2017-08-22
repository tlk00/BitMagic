COMMON_DFLAGS = $(DBGFLAGS)
LINKER_DFLAGS = $(DBGFLAGS)


ifeq ($(COMPILER),GNU_CC)
	PLATFORM_CXXFLAGS = -D_REENTRANT -DBM64OPT -DA_OSF -DSS_64BIT_SERVER -D__USE_STD_IOSTREAM
	PLATFORM_CFLAGS = -D_REENTRANT
	PLATFORM_LDFLAGS =
	EXTERN_LIBS = $(EXTERN_LIBS_BASE)/lib-gnu
	CXX = g++
	CC = g++
	LD = g++
	OPT_FLAGS = -O2
	CC_PIC_FLAGS = -fpic
	CXX_PIC_FLAGS = -fpic
	SO_FLAGS = -shared -Wl,-expect_unresolved,*
	SO_LIBS =
else
	PLATFORM_CXXFLAGS = -pthread -DDIGITAL_UNIX=0x510 -DBM64OPT -DA_OSF -ieee -noimplicit_include -D__USE_STD_IOSTREAM
	PLATFORM_CFLAGS = -pthread -ieee
	PLATFORM_LDFLAGS = -pthread -ieee -use_non_shared_libcxx
	PLATFORM_CLDFLAGS = -pthread -ieee
	EXTERN_LIBS = $(EXTERN_LIBS_BASE)/lib
	CXX = cxx
	CC = cc
#	LD = cc /usr/lib/cmplrs/cxx/_main.o
	LD = cxx
	OPT_FLAGS = -O2 -arch host
	CXX_PIC_FLAGS =
	CC_PIC_FLAGS =
	SO_FLAGS = -shared -Wl,-expect_unresolved,*
	SO_LIBS =
	CXX_REPOSITORY = cxx_repository
	CXX_CACHE = cxx_repository
	COMPILE.cc = $(CXX) $(CXXFLAGS) -MDf.deps/$(*F).pp -c
        DEPCOMMAND:=$(CXX) -MDfdepend.d~ -c
	MAKEDEPS = \
		-egrep -v '\/usr\/' .deps/$(*F).pp > .deps/$(*F).P; \
		tr ' ' '\012' < .deps/$(*F).pp \
		  | sed -e 's/^\\$$//' -e '/^$$/ d' -e '/:$$/ d' -e 's/$$/ :/' \
		    | egrep -v '\/usr\/' >> .deps/$(*F).P; \
		rm .deps/$(*F).pp
endif

#SYS_LIBS = -lpthread -lrt -no_so -L/usr/lib/cmplrs/cxx -lcxxstd -lcxx -so_archive -lexc -lm
SYS_LIBS = -lpthread -lrt -lexc -lm
AR = ar
ARFLAGS = cr
DAEMON_LIBS = -lmach


COMMON_CXXFLAGS = $(COMMON_DFLAGS) $(PLATFORM_CXXFLAGS)
COMMON_CFLAGS = $(COMMON_DFLAGS)  $(PLATFORM_CFLAGS)
COMMON_LDFLAGS = $(LINKER_DFLAGS) $(PLATFORM_LDFLAGS)
COMMON_CLDFLAGS = $(LINKER_DFLAGS) $(PLATFORM_CLDFLAGS)

INSTALL = /usr/local/bin/install
INSTALLDIR = /usr/local/bin/install -d
LDD = ldd 2>&1
AWK = awk
TEST = /bin/test

