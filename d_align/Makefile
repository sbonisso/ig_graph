include ../Makefile.inc

CC=g++
CFLAGS+= -I.. -DCPPPROG
DBGFLAGS+= -I.. -DCPPPROG
LDFLAGS=
SOURCES = $(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
#INC	= -I. -I ${CURDIR}/../ext/seqan-library-2.0.0/include/
INC	= -I. -I ${CURDIR}/../ext/

all: $(OBJECTS) 

subdirs : $(OBJECTS)

#subdirs_debug : CFLAGS	+= -DDEBUG -g
subdirs_debug : CFLAGS	= $(DBGFLAGS)
subdirs_debug : $(OBJECTS)

.cpp.o:
	$(ECHO) $(CC) $(INC) $(CFLAGS) $< -o $@
	$(CC) $(INC) $(CFLAGS) $< -o $@
