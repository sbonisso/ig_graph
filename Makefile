include Makefile.inc

DIRS	= graphs seq_utils file_io d_align
EXE	= iggraph

#INC	= -I. -I ${CURDIR}/ext/seqan-library-2.0.0/include/ -I ${CURDIR}/ext/
INC	= -I. -I ${CURDIR}/ext/seqan-include/ -I ${CURDIR}/ext/
LIBS	= -lemon -lboost_regex -lpthread
SOURCES = $(wildcard graphs/*.cpp) $(wildcard seq_utils/*.cpp) $(wildcard file_io/*.cpp)  $(wildcard d_align/*.cpp)

OBJECTS=$(SOURCES:.cpp=.o)
TESTDIR = tests
TESTSRC	= $(wildcard $(TESTDIR)/*.cpp)
TESTOBJ	= $(TESTSRC:.cpp=.o)

SLIBS = -lpthread -static -static-libgcc -static-libstdc++

all : subdirs $(EXE)

subdirs : 
	-for d in $(DIRS); do (cd $$d; $(MAKE) $@ ); done

$(EXE) : RunIgGraph.o $(OBJECTS)
	$(ECHO) $(CC) $(INC) -o $(EXE) $(OBJECTS) RunIgGraph.o $(LIBS)
	$(CC) $(INC) -o $(EXE) $(OBJECTS) RunIgGraph.o $(LIBS)

static : subdirs RunIgGraph.o $(OBJECTS)
	$(ECHO) $(CC) $(INC) -o $(EXE) $(OBJECTS) RunIgGraph.o $(LIBS)
	$(CC) $(INC) -o $(EXE) $(OBJECTS) RunIgGraph.o $(LIBS) $(SLIBS)

RunIgGraph.o: RunIgGraph.cpp
	$(ECHO) $(CC) $(INC) $(CFLAGS) $< -o $@
	$(CC) $(INC) $(CFLAGS)  $< -o $@

debug : CFLAGS = $(DBGFLAGS)
debug :	subdirs_debug $(EXE)

subdirs_debug : 
	-for d in $(DIRS); do (cd $$d; $(MAKE) $@ ); done

test : $(EXE)
	cd $(TESTDIR); $(MAKE) $@
	$(ECHO) $(CC) -o unit_test $(TESTOBJ) $(OBJECTS) -lcpptest $(LIBS)
	$(CC) -o unit_test $(TESTOBJ) $(OBJECTS) -lcpptest $(LIBS)

clean : 
	$(ECHO) cleaning up
	$(RM) RunIgGraph.o
	$(RM) $(OBJECTS)	
	$(RM) $(EXE)
	$(RM) $(TESTOBJ)
	$(RM) unit_test
