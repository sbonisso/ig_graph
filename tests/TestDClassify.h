#ifndef TESTDCLASSIFY_H_
#define TESTDCLASSIFY_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>

#include "../prettyprint.hpp"

#include "file_io/FastaParser.h"
#include "file_io/FastaMult.h"
#include "file_io/FastaRefID.h"
#include "file_io/FastaWriter.h"		// for writing FASTA files
#include "file_io/FastaFlowParser.h"

#include "d_align/DClassify.hpp"

using namespace std;

class TestDClassify : public Test::Suite {
public:
    TestDClassify();
    virtual ~TestDClassify();

protected:
    virtual void setup();      		// setup resources...
    virtual void tear_down();
    
    string d_fasta = "data/igh_refs_simple/human_IGHD.fa";
    
    void test_d_align_id();
    void test_d_align_score();
    
    void test_d_align_nohit();
};

#endif
