COMPILER = GNU_CC

COMMON_DFLAGS = 
#-g -D_DEBUG
LINKER_DFLAGS = 
#-g

#	OS_VER = -D__$(shell uname -s)_$(shell uname -r | sed -e 's/\./_/g' -e 's/-.*//')
OS_VER = -DCYGWIN_NT
PLATFORM_CXXFLAGS = -D_REENTRANT $(OS_VER) -D_GNU_SOURCE  -Wno-deprecated
PLATFORM_CFLAGS = -D_REENTRANT $(OS_VER)
COMMON_LDFLAGS = $(LINKER_DFLAGS) -export-dynamic -L/usr/local/lib
COMMON_CLDFLAGS = $(COMMON_LDFLAGS)
EXTERN_LIBS = $(EXTERN_LIBS_BASE)/lib
CXX = g++ -msse4 -march=core2 -Wall
CC = gcc  -Wall
LD = g++
CC_PIC_FLAGS = -fPIC
CXX_PIC_FLAGS = -fPIC
OPT_FLAGS = -g0 -O2 -march=core2 -fomit-frame-pointer -pipe
SO_FLAGS = -shared
SO_LIBS =

SYS_LIBS = -lpthread 

COMMON_CXXFLAGS = $(COMMON_DFLAGS) $(PLATFORM_CXXFLAGS)
COMMON_CFLAGS = $(COMMON_DFLAGS) $(PLATFORM_CFLAGS)
INSTALL = /usr/bin/install
INSTALLDIR = /usr/bin/install -d
AWK = awk
TEST = /usr/bin/test
