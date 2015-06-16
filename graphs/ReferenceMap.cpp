#include "ReferenceMap.hpp"


ReferenceMap:: ReferenceMap() {}
/**
 *
 */    
bool ReferenceMap::addReference(string ref_id, string seq) {
    if( h_.find(ref_id) != h_.end() ) { return false; }
    h_[ref_id] = unordered_map<string,int>();
    ref_h_[ref_id] = seq;
    return true;
}
/**
 *
 */
bool ReferenceMap::addKmerToReference(string ref_id, string kmer, int ind) {
    if( h_.find(ref_id) == h_.end() ) { return false; }
    h_[ref_id][kmer] = ind;
    return true;
}
/**
 *
 */
int ReferenceMap::getIndex(string ref_id, string kmer) {
    if( h_[ref_id].find(kmer) == h_[ref_id].end() ) { return -1; }
    return h_[ref_id][kmer];
}
/**
 *
 */
int ReferenceMap::size() { return (int)h_.size(); }
/**
 *
 */
string ReferenceMap::getSeq(string ref_id) { 
    if( ref_h_.find(ref_id) == ref_h_.end() ) { return ""; }
    return  ref_h_[ref_id];
}
/**
 * return a vector of reference IDs, i.e., the keys of the hash
 */
vector<string> ReferenceMap::getIDs() {
    vector<string> ref_ids;
    ref_ids.reserve(this->h_.size());
    for( const auto &pair : this->h_ ) {
	ref_ids.push_back(pair.first);
    }
    // now sort so consistent across object creations
    std::sort(ref_ids.begin(), ref_ids.end());
    return ref_ids;
}
