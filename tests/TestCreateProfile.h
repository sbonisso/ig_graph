#ifndef TESTCREATEPROFILE_H_
#define TESTCREATEPROFILE_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>
#include <sstream>

#include "file_io/FastaParser.h"
#include "file_io/FastaMult.h"
#include "file_io/FastaRefID.h"
#include "file_io/FastaWriter.h"		// for writing FASTA files
#include "file_io/FastaFlowParser.h"

#include "graphs/ReferenceMap.h"
#include "graphs/CanonicalAntibodyGraph.h"
#include "graphs/CreateProfile.h"

using namespace std;

class TestCreateProfile : public Test::Suite {
public:
    TestCreateProfile();
    virtual ~TestCreateProfile();

protected:
    virtual void setup();      		// setup resources...
    virtual void tear_down(); 		// remove resources...
    
    void test_cp_preds_1();
    void test_cp_scores_1();

    void test_cp_preds_2();
    void test_cp_scores_2();
    
    void test_max_n();
    void test_cdr3();

    CanonicalAntibodyGraph *cab;
};

#endif
