include ../Makefile.inc

CC=g++
CFLAGS+= -I..
DBGFLAGS+= -I..
LDFLAGS=
SOURCES = $(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(OBJECTS) 

subdirs : $(OBJECTS)

#subdirs_debug : CFLAGS	+= -DDEBUG -g
subdirs_debug : CFLAGS	= $(DBGFLAGS)
subdirs_debug : $(OBJECTS)

.cpp.o:
	$(ECHO) $(CC) $(CFLAGS) $< -o $@
	$(CC) $(CFLAGS) $< -o $@
