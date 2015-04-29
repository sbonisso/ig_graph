#ifndef TESTCANONICALANTIBODYGRAPH_H_
#define TESTCANONICALANTIBODYGRAPH_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>
#include <sstream>

#include "file_io/FastaParser.h"
#include "file_io/FastaMult.h"
#include "file_io/FastaRefID.h"
#include "file_io/FastaWriter.h"	// for writing FASTA files
#include "file_io/FastaFlowParser.h"

#include "graphs/ReferenceMap.h"
#include "graphs/CanonicalAntibodyGraph.h"

using namespace std;

class TestCanonicalAntibodyGraph : public Test::Suite {
public:
    TestCanonicalAntibodyGraph();
    virtual ~TestCanonicalAntibodyGraph();

protected:
    virtual void setup();      		// setup resources...
    
    void test_paint_1();
    void test_paint_2();
    void test_paint_2_color_propagate();
    void test_paint_2_simple();

    void test_paint_1_color_propagate_tip_j_start();
    void test_paint_1_color_propagate_tip_j_end();
    
    void test_paint_2_color_propagate_tip1();
    void test_paint_2_color_propagate_tip2();
};

#endif
