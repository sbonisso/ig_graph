#ifndef CANONICALANTIBODYGRAPH_H_
#define CANONICALANTIBODYGRAPH_H_

#include <string>
#include <vector>
#include <ctime>
#include <time.h>
#include <unordered_map>
#include "prettyprint.hpp"

#include "file_io/FastaParser.hpp"  		// for parsing FASTA files
#include "file_io/FastaRefID.hpp"

#include "graphs/ReferenceMap.hpp"

#include "tests/TestCanonicalAntibodyGraph.hpp"

#ifdef DEBUG
#define CANONAB_DEBUG_PRINT(EXPR) DEBUG_PRINT("[CANON_AB]", EXPR)
#else
#define CANONAB_DEBUG_PRINT(x) do {} while (0)
#endif

using namespace std;

enum class Segment { V_GENE, D_GENE, J_GENE };

class CanonicalAntibodyGraph {
public:
    CanonicalAntibodyGraph();
    CanonicalAntibodyGraph(int kval);
    virtual ~CanonicalAntibodyGraph() {}

    virtual void addVReferences(string v_fasta);
    virtual void addDReferences(string d_fasta);
    virtual void addJReferences(string j_fasta);    
    //
    virtual vector<int> getVPainting(string ref_id, string seq);
    virtual vector<int> getDPainting(string ref_id, string seq);
    virtual vector<int> getJPainting(string ref_id, string seq);
    //
    string getVRefID(int index) { return v_ids_[index]; }
    string getDRefID(int index) { return d_ids_[index]; }
    string getJRefID(int index) { return j_ids_[index]; }
    //
    string getVRefSeq(string ref_id) { return v_kmers_.getSeq(ref_id); }
    string getDRefSeq(string ref_id) { return d_kmers_.getSeq(ref_id); }
    string getJRefSeq(string ref_id) { return j_kmers_.getSeq(ref_id); }
    //
    int getNumV() { return this->v_kmers_.size(); }
    int getNumD() { return this->d_kmers_.size(); }
    int getNumJ() { return this->j_kmers_.size(); }
    //
    int getNumReferences();
    int getK();

    vector<string> getVRefs() { return v_ids_; }
    vector<string> getDRefs() { return d_ids_; }
    vector<string> getJRefs() { return j_ids_; }
	
    friend class TestCanonicalAntibodyGraph; // for unit testing
    /* friend TestColorProfileMatrix;     // for unit testing */
    
protected:
    int k_;
    
    ReferenceMap v_kmers_;
    ReferenceMap d_kmers_;
    ReferenceMap j_kmers_;

    vector<string> v_ids_;
    vector<string> d_ids_;
    vector<string> j_ids_;
    
    bool read_into_hash(string fasta_file, ReferenceMap &h, int k);
    
    vector<int> getPainting(string ref_id, string seq, ReferenceMap &h, int k, 
			    bool tip_trim);
    vector<int> initPainting(string ref_id, string seq, ReferenceMap &h, int k);
    
    void propagateColor(string ref_id, string seq, ReferenceMap &h, 
			vector<int> &row);
    int propagateColorBulge(string ref_id, string seq, ReferenceMap &h,
			    vector<int> &row, int start_ind, int end_ind);
    int propagateColorTip(string ref_id, string seq, 
			  ReferenceMap &h, vector<int> &row, bool tip_trim);
    void pushColor(int start_ind, int end_ind,
		   string &read_substr, string &ref_substr,
		   vector<int> &row, int ref_index, bool is_tip);

    bool hasAtLeastOneMatch(vector<int> &row);
};


#endif
