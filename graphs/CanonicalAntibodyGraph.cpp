#include "CanonicalAntibodyGraph.hpp"
/**
 *
 */
CanonicalAntibodyGraph::CanonicalAntibodyGraph() {
    k_ = 21;
}
CanonicalAntibodyGraph::CanonicalAntibodyGraph(int kval) {
    k_ = kval;
}
/**
 * add references from specified FASTA file
 */
void CanonicalAntibodyGraph::addVReferences(string v_fasta) {
    this->read_into_hash(v_fasta, this->v_kmers_, k_);
    vector<string> ref_ids = this->v_kmers_.getIDs();
    for(int i = 0; i < (int)ref_ids.size(); i++) {
	v_ids_.push_back(ref_ids[i]);
    }
}
/**
 *
 */
void CanonicalAntibodyGraph::addDReferences(string d_fasta) {
    this->read_into_hash(d_fasta, this->d_kmers_, k_);
    vector<string> ref_ids = this->d_kmers_.getIDs();
    for(int i = 0; i < (int)ref_ids.size(); i++) {
	d_ids_.push_back(ref_ids[i]);
    }
}
/**
 *
 */
void CanonicalAntibodyGraph::addJReferences(string j_fasta) {
    this->read_into_hash(j_fasta, this->j_kmers_, k_);
    vector<string> ref_ids = this->j_kmers_.getIDs();
    for(int i = 0; i < (int)ref_ids.size(); i++) {
	j_ids_.push_back(ref_ids[i]);
    }
}
/**
 * go through FASTA, populating the given hash with k-mers => index
 */
bool CanonicalAntibodyGraph::read_into_hash(string fasta_file, 
					    ReferenceMap &h,
					    int k) {
    FastaRefID<FastaParser> fp(fasta_file);
    fp.openFile();
	
    std::locale loc;
    while(fp.hasNextSequence()) {
	pair<string,string> entry = fp.getNextEntry();
	
	string readID = fp.getCurrIDLine().substr(1);  // remove the '>'
	string seq = entry.second;
	
	h.addReference(readID, seq);
	
	int len = ((int)seq.size() - k);
	for(int i = 0; i < len; i++) {
	    string curr_kmer = seq.substr(i, k);
	    h.addKmerToReference(readID, curr_kmer, i);
	}
    }
    fp.closeFile();
    return true;
}
/**
 * return 'painting' of a read to a particular reference
 */
vector<int> CanonicalAntibodyGraph::getVPainting(string ref_id, string seq) {
    return getPainting(ref_id, seq, v_kmers_, k_, false);
}
vector<int> CanonicalAntibodyGraph::getDPainting(string ref_id, string seq) {  
    return getPainting(ref_id, seq, d_kmers_, k_, true);
}
vector<int> CanonicalAntibodyGraph::getJPainting(string ref_id, string seq) {
    return getPainting(ref_id, seq, j_kmers_, k_, true);
}
/**
 * generic method for painting a read and propagating the color
 */
vector<int> CanonicalAntibodyGraph::getPainting(string ref_id, 
						string seq,
						ReferenceMap &h,
						int k,
						bool tip_trim) {    
    vector<int> v = initPainting(ref_id, seq, h, k);
    //CANONAB_DEBUG_PRINT("PROPAGATE:\t"<<v);
    if(!hasAtLeastOneMatch(v)) { 
	//CANONAB_DEBUG_PRINT("NO MATCH!\n");
	return v; 
    }    
    propagateColor(ref_id, seq, h, v);
    //CANONAB_DEBUG_PRINT("PROPAGATE2:\t"<<v);
    propagateColorTip(ref_id, seq, h, v, tip_trim);
    //CANONAB_DEBUG_PRINT("PROPAGATE3:\t"<<v);
    return v;
}
/**
 *
 */
bool CanonicalAntibodyGraph::hasAtLeastOneMatch(vector<int> &row) {
    int len = (int)row.size();
    for(int i = 0; i < len; i++) { 
	if(row[i] >= 0) return true;
    }
    return false;
}
/**
 * initial painting of the read based on shared edges in Ab graph
 */
vector<int> CanonicalAntibodyGraph::initPainting(string ref_id, 
						 string seq,
						 ReferenceMap &h,
						 int k) {
    //int len = ((int)seq.size()-k);
    int len = (int)seq.size();
    vector<int> v(len,-1);
    for(int i = 0; i <  len; i++) {
	string curr_kmer = seq.substr(i, k);
	v[i] = h.getIndex(ref_id, curr_kmer);	
    }    
    return v;
}
/**
 * propagates single color (reference) for a row of the color matrix
 * should take as input the output from getPainting()
 */
void CanonicalAntibodyGraph::propagateColor(string ref_id, 
					    string seq,
					    ReferenceMap &h,
					    vector<int> &row) {
    // scan painting for gaps, i.e., -1 with begining and ending (a bulge)
    int start_ind = -1;
    int end_ind = -1;
    int len = (int)row.size();
    for(int i = 1; i < len-1; i++) {
	if(row[i] < 0 && row[i-1] > 0 && start_ind < 0) { start_ind = i; }
	else if(row[i] < 0 && row[i+1] > 0 && 
		start_ind > 0 && end_ind < 0) { end_ind = i; }
	else if(row[i] > 0 && start_ind > 0 && end_ind > 0) { 
	    // prop color
	    propagateColorBulge(ref_id, seq, h, row, 
				start_ind-1, end_ind+1);
	    // reinstate flags
	    start_ind = -1;
	    end_ind = -1;
	}
	else {}
    }
}
/**
 * given the start/end index of a bulge, propagate the color in the row
 */
int CanonicalAntibodyGraph::propagateColorBulge(string ref_id,
						string seq,
						ReferenceMap &h,
						vector<int> &row, 
						int start_ind, 
						int end_ind) {
    string ref_seq = h.getSeq(ref_id);
    int bulge_len_ref = (row[end_ind] - row[start_ind]);
    //CANONAB_DEBUG_PRINT("BULGE LEN:\t"<<bulge_len_ref);
    if(bulge_len_ref < 0) { return -1; }
    //
    string ref_substr( ref_seq.substr(row[start_ind], bulge_len_ref) );
    string read_substr( seq.substr(start_ind, (end_ind-start_ind)) );
    //CANONAB_DEBUG_PRINT("SUBSTRS:\t"<<ref_substr<<"\n\t"<<read_substr);
    if( read_substr.size() != ref_substr.size() ) { return -2; }
    //
    int init_ref_pos = row[start_ind-1];
    int ref_index = init_ref_pos;
    //
    pushColor(start_ind, end_ind, read_substr, ref_substr, row, ref_index, false);
    return 0;
}
/**
 * helper method that 'pushes' color along a bulge/tip as specified
 */
void CanonicalAntibodyGraph::pushColor(int start_ind, int end_ind,
				       string &read_substr, string &ref_substr,
				       vector<int> &row, int ref_index, 
				       bool is_tip) {
    int index = 0;
    // CANONAB_DEBUG_PRINT("READ IND:\t"<<start_ind<<"\t"<<end_ind);
    // CANONAB_DEBUG_PRINT("REF IND:\t"<<ref_index);
    // CANONAB_DEBUG_PRINT("REF:\t"<<ref_substr);
    // CANONAB_DEBUG_PRINT("READ:\t"<<read_substr);
    for(int i = start_ind; i < end_ind; i++) {
	if(is_tip && ref_substr[index] != read_substr[index]) { break; }
	if(ref_substr[index] == read_substr[index]) {
	    row[i] = ref_index;
	}
	ref_index++;
	index++;
    }
}
/**
 * propagate color at tips at either end of the reference (found by the 
 * first/last occurence within the row of the read)
 */
int CanonicalAntibodyGraph::propagateColorTip(string ref_id, string seq,
					      ReferenceMap &h, vector<int> &row,
					      bool tip_trim) {
    int len = (int)row.size();
    // no tips on read
    if(row[0] >= 0 && row[len-1] >= 0) { return -1; }
    string ref_seq = h.getSeq(ref_id);
    int ref_len = (int)ref_seq.size();
    int num_tips = 0;
    if(row[0] < 0) {
	int end_ind = 0;
	while(end_ind < len-1 && row[end_ind] < 0) { end_ind++; }
	//CANONAB_DEBUG_PRINT("END_IND:\t"<<end_ind<<"\tLEN:\t"<<row.size());
	int ref_pos = row[end_ind];
	int curr_len = ref_pos;
	int start_ind = end_ind - ref_pos;  // project ref start to read
	if(start_ind >= 0) {
	    //string read_tip( seq.substr(0, end_ind) );
	    string read_tip( seq.substr(start_ind, curr_len) );
	    string ref_tip( ref_seq.substr(0, row[end_ind]) );
	    int ref_index = row[end_ind]-ref_pos;	
	    // if doesn't start at beginning, propagate
	    if(start_ind < end_ind) {
		pushColor(start_ind, end_ind, read_tip, ref_tip, row, ref_index, tip_trim); //true);
		num_tips++;
	    }
	}
    }
    if(row[len-1] < 0) {
    	int start_ind = len-1;
    	while(start_ind > 0 && row[start_ind] < 0) { start_ind--; }
	//int curr_len = ref_len-row[start_ind];
	int curr_len = std::min(ref_len-row[start_ind], len-start_ind);
	int end_ind = start_ind+curr_len;    	
    	string read_tip( seq.substr(start_ind, curr_len) );
    	string ref_tip( ref_seq.substr(row[start_ind], curr_len) );
    	int ref_index = row[start_ind];
	// CANONAB_DEBUG_PRINT("READ IND:\t"<<start_ind<<"\t"<<end_ind);
	// CANONAB_DEBUG_PRINT("READ TIP:\t"<<read_tip<<"\nREF TIP:\t"<<ref_tip);
    	pushColor(start_ind, end_ind, read_tip, ref_tip, row, ref_index, tip_trim); //true);
    	num_tips++;
    }
    
    return num_tips;
}
/**
 * get the total number of references, |V|+|D|+|J|
 */
int CanonicalAntibodyGraph::getNumReferences() {
    return (v_kmers_.size() + 
	    d_kmers_.size() + 
	    j_kmers_.size() );
}
/**
 * return value of k used
 */
int CanonicalAntibodyGraph::getK() { return k_; }
