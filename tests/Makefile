include ../Makefile.inc

CC=g++
CFLAGS=-c -g -Wall -std=c++11 -I.. -I../ext/seqan-library-2.0.0/include/ -I../ext/
LDFLAGS=
SOURCES = $(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(OBJECTS) 

test : $(OBJECTS)

.cpp.o:
	$(ECHO) $(CC) $(CFLAGS) $< -o $@
	$(CC) $(CFLAGS) $< -o $@
