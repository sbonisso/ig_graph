include ../Makefile.inc

CC=g++
CFLAGS+= -I.. -DCPPPROG
DBGFLAGS+= -I.. -DCPPPROG
LDFLAGS=
SOURCES = $(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
INC	= -I. -I ${CURDIR}/../ext/ -I ${CURDIR}/../ext/cpptest/cpptest/src/

all: $(OBJECTS) 

subdirs : $(OBJECTS)

#subdirs_debug : CFLAGS	+= -DDEBUG -g
subdirs_debug : CFLAGS	= $(DBGFLAGS)
subdirs_debug : $(OBJECTS)

.cpp.o:
	$(ECHO) $(CC) $(INC) $(CFLAGS) $< -o $@
	$(CC) $(INC) $(CFLAGS) $< -o $@
