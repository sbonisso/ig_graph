#ifndef TESTCANONICALANTIBODYGRAPH_H_
#define TESTCANONICALANTIBODYGRAPH_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>
#include <sstream>

#include "prettyprint.hpp"
#include "file_io/FastaParser.hpp"
#include "file_io/FastaMult.hpp"
#include "file_io/FastaRefID.hpp"
#include "file_io/FastaWriter.hpp"	// for writing FASTA files

#include "graphs/ReferenceMap.hpp"
#include "graphs/CanonicalAntibodyGraph.hpp"

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
