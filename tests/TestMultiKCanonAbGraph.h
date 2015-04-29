#ifndef TESTMULTIKCANONABGRAPH_H_
#define TESTMULTIKCANONABGRAPH_H_

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

#include "graphs/MultiKCanonAbGraph.h"
#include "graphs/CanonicalAntibodyGraph.h"

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
