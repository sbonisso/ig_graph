#ifndef TESTCOLORPROFILEMATRIX_H_
#define TESTCOLORPROFILEMATRIX_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>
#include <sstream>

#include "prettyprint.hpp"
#include "file_io/FastaParser.hpp"
#include "file_io/FastaRefID.hpp"

#include "graphs/ReferenceMap.hpp"
#include "graphs/CanonicalAntibodyGraph.hpp"
#include "graphs/CreateProfile.hpp"

using namespace std;

class TestColorProfileMatrix : public Test::Suite {
public:
    TestColorProfileMatrix();
    virtual ~TestColorProfileMatrix();

protected:
    virtual void setup();      		// setup resources...
    virtual void tear_down(); 		// remove resources...
    
    void test_cp_matrix_scores_1();
    void test_cp_matrix_scores_2();
    
    void test_cp_matrix_partition_1();
    
    CanonicalAntibodyGraph *cab;
};

#endif
