#include "MultiKCanonAbGraph.h"

MultiKCanonAbGraph::MultiKCanonAbGraph() {
    v_k_ = 21;
    d_k_ = 21;
    j_k_ = 21;
}

MultiKCanonAbGraph::MultiKCanonAbGraph(int v_k, int d_k, int j_k) {
    v_k_ = v_k;
    d_k_ = d_k;
    j_k_ = j_k;
}
/**
 *
 */
void MultiKCanonAbGraph::addVReferences(string v_fasta) {
    k_ = v_k_;
    CanonicalAntibodyGraph::addVReferences(v_fasta);
}
void MultiKCanonAbGraph::addDReferences(string d_fasta) {
    k_ = d_k_;
    CanonicalAntibodyGraph::addDReferences(d_fasta);
}
void MultiKCanonAbGraph::addJReferences(string j_fasta) {
    k_ = j_k_;
    CanonicalAntibodyGraph::addJReferences(j_fasta);
}
/**
 * return 'painting' of a read to a particular reference
 */
vector<int> MultiKCanonAbGraph::getVPainting(string ref_id, string seq) {
    return getPainting(ref_id, seq, v_kmers_, v_k_, false);
}
vector<int> MultiKCanonAbGraph::getDPainting(string ref_id, string seq) {  
    return getPainting(ref_id, seq, d_kmers_, d_k_, true);
}
vector<int> MultiKCanonAbGraph::getJPainting(string ref_id, string seq) {
    return getPainting(ref_id, seq, j_kmers_, j_k_, true);
}

