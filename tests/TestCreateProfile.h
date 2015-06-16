#ifndef TESTCREATEPROFILE_H_
#define TESTCREATEPROFILE_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>
#include <sstream>

#include "file_io/FastaParser.hpp"
#include "file_io/FastaMult.hpp"
#include "file_io/FastaRefID.hpp"
#include "file_io/FastaWriter.hpp"		// for writing FASTA files

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
    
    void test_cdr3_1();
    void test_cdr3_2();
    
    CanonicalAntibodyGraph *cab;
};

#endif
