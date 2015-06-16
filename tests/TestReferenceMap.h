#ifndef TESTREFERENCEMAP_H_
#define TESTREFERENCEMAP_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>

#include "file_io/FastaParser.hpp"
#include "file_io/FastaMult.hpp"
#include "file_io/FastaRefID.hpp"
#include "file_io/FastaWriter.hpp"		// for writing FASTA files

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
