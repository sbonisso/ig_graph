#ifndef TESTREFERENCEMAP_H_
#define TESTREFERENCEMAP_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>

#include "file_io/FastaParser.h"
#include "file_io/FastaMult.h"
#include "file_io/FastaRefID.h"
#include "file_io/FastaWriter.h"		// for writing FASTA files
#include "file_io/FastaFlowParser.h"

#include "graphs/ReferenceMap.h"

using namespace std;

class TestReferenceMap : public Test::Suite {
public:
    TestReferenceMap();
    virtual ~TestReferenceMap();

protected:
    virtual void setup();      		// setup resources...
//    virtual void tear_down() = 0; 		// remove resources...

    void test_ref_map_1();
};

#endif
