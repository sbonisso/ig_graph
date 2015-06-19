include ../Makefile.inc

CC=g++
#CFLAGS=-c -g -Wall -std=c++11 -I.. -I../ext/seqan-library-2.0.0/include/ -I../ext/cereal/include/
CFLAGS=-c -g -Wall -std=c++11 
INC	= -I.. -I../ext/ -I../ext/cereal/include/ -I ${CURDIR}/../ext/cpptest/cpptest/src/
LDFLAGS=
SOURCES = $(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(OBJECTS) 

test : $(OBJECTS)

.cpp.o:
	$(ECHO) $(CC) $(INC) $(CFLAGS) $< -o $@
	$(CC) $(INC) $(CFLAGS) $< -o $@
