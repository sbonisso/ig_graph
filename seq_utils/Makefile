include ../Makefile.inc

CC=g++
#CFLAGS=-c -Wall -std=c++11
CFLAGS += -I..
DBGFLAGS+= -I..
LDFLAGS=
SOURCES = $(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
INC	= -I.. -I ${CURDIR}/../ext/cereal/include/

all: $(OBJECTS) 

subdirs : $(OBJECTS)

#subdirs_debug : CFLAGS	+= -DDEBUG -g
subdirs_debug : CFLAGS	= $(DBGFLAGS)
subdirs_debug : $(OBJECTS)

.cpp.o:
	$(ECHO) $(CC) $(INC) $(CFLAGS) $< -o $@
	$(CC) $(INC) $(CFLAGS) $< -o $@
