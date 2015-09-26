#ifndef TESTCABGRAPH_H_
#define TESTCABGRAPH_H_

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
//#include "graphs/CreateProfile.hpp"

using namespace std;

class TestCABGraph : public CanonicalAntibodyGraph {
public:
    TestCABGraph(); 
    TestCABGraph(int kval);
    vector<int> run_test(int kval,
    			 string ref_fasta,
    			 string read_fasta_file,
    			 int index, 
    			 char to_ch,
    			 string ref_id);
};

#endif
