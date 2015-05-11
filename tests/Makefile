include ../Makefile.inc

CC=g++
CFLAGS=-c -Wall -std=c++11 -I.. -I../ext/seqan-library-2.0.0/include/
LDFLAGS=
SOURCES = $(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(OBJECTS) 

test : $(OBJECTS)

.cpp.o:
	$(ECHO) $(CC) $(CFLAGS) $< -o $@
	$(CC) $(CFLAGS) $< -o $@
