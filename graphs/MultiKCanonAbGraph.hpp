#ifndef MULTIKCANONABGRAPH_H_
#define MULTIKCANONABGRAPH_H_

#include <string>
#include <vector>
#include <ctime>
#include <time.h>
#include <unordered_map>
#include "prettyprint.hpp"

#include "Utils.h"

#include "graphs/CanonicalAntibodyGraph.hpp"

#ifdef DEBUG
#define MULTIKAB_DEBUG_PRINT(EXPR) DEBUG_PRINT("[MULTI-K]", EXPR)
#else
#define MULTIKAB_DEBUG_PRINT(x) do {} while (0)
#endif

using namespace std;

class MultiKCanonAbGraph : public CanonicalAntibodyGraph {
public:
    MultiKCanonAbGraph();
    MultiKCanonAbGraph(int v_k, int d_k, int j_k);
    //virtual ~MultiKCanonAbGraph() {}
    //
    virtual void addVReferences(string v_fasta);
    virtual void addDReferences(string d_fasta);
    virtual void addJReferences(string j_fasta);    
    //
    virtual vector<int> getVPainting(string ref_id, string seq);
    virtual vector<int> getDPainting(string ref_id, string seq);
    virtual vector<int> getJPainting(string ref_id, string seq);
    //
    int getVK() { return v_k_; }
    int getDK() { return d_k_; }
    int getJK() { return j_k_; }
    
protected:
    
    int v_k_; int d_k_; int j_k_;    
};


#endif
