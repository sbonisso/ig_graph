include ../Makefile.inc

CC=g++
CFLAGS=-c -g -Wall -std=c++11
INC	= -I.. -I../ext/ -I../ext/cereal/include/ -I ${CURDIR}/../graphs/ -I ${CURDIR}/../file_io/
LDFLAGS=
SOURCES = $(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)

all: $(OBJECTS) 

test : $(OBJECTS)

.cpp.o:
	$(ECHO) $(CC) $(INC) $(CFLAGS) $< -o $@
	$(CC) $(INC) $(CFLAGS) $< -o $@
