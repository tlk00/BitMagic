COMMON_DFLAGS = $(DBGFLAGS)
LINKER_DFLAGS = $(DBGFLAGS)

ARCH=$(shell uname -m)

ifeq ($(COMPILER),GNU_CC)
   LINKER_DFLAGS =
   #-g
   OS_VER = -D__$(shell uname -s)_$(shell uname -r | sed -e 's/\./_/g' -e 's/-.*//')
   PLATFORM_CXXFLAGS = -D_REENTRANT $(OS_VER) -D_GNU_SOURCE
   PLATFORM_CFLAGS = -D_REENTRANT $(OS_VER)
   COMMON_LDFLAGS = $(LINKER_DFLAGS) -export-dynamic
   COMMON_CLDFLAGS = $(COMMON_LDFLAGS)
   EXTERN_LIBS = $(EXTERN_LIBS_BASE)/lib
   CXX = g++ -Wall
   CC = gcc -Wall
   LD = g++
   CC_PIC_FLAGS = -fPIC
   CXX_PIC_FLAGS = -fPIC
   OPT_FLAGS = -g0 -O2
   SO_FLAGS = -shared
   SO_LIBS =
   SYS_LIBS = -lpthread -ldl -lrt
   COMMON_CXXFLAGS = $(COMMON_DFLAGS) $(PLATFORM_CXXFLAGS)
   COMMON_CFLAGS = $(COMMON_DFLAGS) $(PLATFORM_CFLAGS)
   INSTALL = /usr/bin/install
   INSTALLDIR = /usr/bin/install -d
   AWK = awk
   TEST = /usr/bin/test
else
   # TO DO: add Xlc support here
endif