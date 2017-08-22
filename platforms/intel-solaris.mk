COMMON_DFLAGS = $(DBGFLAGS)
LINKER_DFLAGS = $(DBGFLAGS)

ARCH=$(shell uname -m)

ifeq ($(COMPILER),GNU_CC)
	OS_VER = -D__$(shell uname -s)_$(shell uname -r | sed 's/\./_/g')
	PLATFORM_CXXFLAGS = -D_REENTRANT $(OS_VER)
	PLATFORM_CFLAGS = -D_REENTRANT $(OS_VER)
	COMMON_LDFLAGS = $(LINKER_DFLAGS)
	COMMON_CLDFLAGS = $(COMMON_LDFLAGS)
	EXTERN_LIBS = $(EXTERN_LIBS_BASE)/lib-gnu
	CXX = g++ -mcpu=v8 -Wall
	CC = g++ -mcpu=v8 -Wall
	LD = g++
	CC_PIC_FLAGS = -fPIC
	CXX_PIC_FLAGS = -fPIC
	OPT_FLAGS = -g0 -O2
	SO_FLAGS = -shared
	SO_LIBS =
	LDSO = gcc
else
	CC_VER = $(shell CC -V 2>&1 | sed -e 's/.\../|b&|e/' -e 's/.*|b//' -e 's/|e.*//')
	PLATFORM_CXXFLAGS = -mt 
	PLATFORM_CFLAGS =
	ifneq ($(CC_VER),5.0)
		PLATFORM_CXXFLAGS += -template=no%extdef -D_RWSTD_COMPILE_INSTANTIATE
	endif
	COMMON_LDFLAGS = $(LINKER_DFLAGS) -xildoff -mt -staticlib=Crun
	COMMON_CLDFLAGS = $(LINKER_DFLAGS) -xildoff -mt
	EXTERN_LIBS = $(EXTERN_LIBS_BASE)/lib
	ULTRASPARC=sun4u
	ifeq ($(ARCH),$(ULTRASPARC))
		ifeq ($(SUN_ARCH),)
			SUN_ARCH = v8plusa
		endif
		CXX = CC -xchip=ultra2 -xarch=$(SUN_ARCH)
		CC = cc -xchip=ultra2 -xarch=$(SUN_ARCH)
		ifeq ($(SUN_ARCH),v9a)
			PLATFORM_CXXFLAGS += -DSS_64BIT_SERVER 
			LD = CC -xchip=ultra2 -xarch=v9
		else
			LD = $(CXX)
		endif
	else
		CXX = CC
		CC = cc
		LD = $(CXX)
	endif
	CC_PIC_FLAGS = -KPIC
	CXX_PIC_FLAGS = -KPIC
	OPT_FLAGS = -g0 -xO2
	SO_FLAGS = -G
	AR = $(CXX) -xar
	ARFLAGS = -o
	CXX_CACHE = SunWS_cache
endif

	SYS_LIBS = -lnsl -lsocket -lpthread -lposix4 -ldl
	DAEMON_LIBS = -lkstat

	COMMON_CXXFLAGS = $(COMMON_DFLAGS) $(PLATFORM_CXXFLAGS)
	COMMON_CFLAGS = $(COMMON_DFLAGS) $(PLATFORM_CFLAGS)
	INSTALL = /usr/local/bin/install
	INSTALLDIR = /usr/sbin/install -d
	AWK = awk
	TEST = /bin/test
