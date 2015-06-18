#ifndef TESTMULTIKCANONABGRAPH_H_
#define TESTMULTIKCANONABGRAPH_H_

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

#include "graphs/MultiKCanonAbGraph.hpp"
#include "graphs/CanonicalAntibodyGraph.hpp"

using namespace std;

class TestMultiKCanonAbGraph : public Test::Suite {
public:
    TestMultiKCanonAbGraph();
    virtual ~TestMultiKCanonAbGraph();

protected:
    virtual void setup();      		// setup resources...
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    string d_fasta = "data/igh_refs_simple/human_IGHD.fa";
    string j_fasta = "data/igh_refs_simple/human_IGHJ.fa";
    
    void test_multi_k_ident();
    void test_multi_k_1();
};

#endif
