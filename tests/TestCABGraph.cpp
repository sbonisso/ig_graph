#include "TestCABGraph.hpp"

TestCABGraph::TestCABGraph() : CanonicalAntibodyGraph() {}

TestCABGraph::TestCABGraph(int kval) : CanonicalAntibodyGraph(kval) {}

vector<int> TestCABGraph::run_test(int kval,
				   string ref_fasta,
				   string read_fasta_file,
				   int index, 
				   char to_ch,
				   string ref_id) {
    this->addVReferences(ref_fasta);
    //
    ifstream in(read_fasta_file);
    if(!in.is_open()) { return vector<int>(); }
    //
    string read_id;  string seq = ""; 
    for(int i = 0; i < 2; i++) {
	getline(in, read_id);
	getline(in, seq);
    }
    //
    seq[index] = to_ch;
    //
    vector<int> v = this->initPainting(ref_id, seq, this->v_kmers_, this->k_);
    this->propagateColorTip(ref_id, seq, this->v_kmers_, v, true);
    return v;
}
