#include "DClassify.hpp"

DClassify::DClassify() {
    match_s_ = 5; 
    mismatch_s_ = -4;
    indel_s_ = -4;
    num_report_ = 2;
}
/**
 *
 */
DClassify::DClassify(std::string dref_fasta) : DClassify() {
    read_in_refs(dref_fasta, drefs_h_);   
}

DClassify::~DClassify() {}
/**
 * set the number of D gene-segments to report
 */
void DClassify::set_num_d_report(int num_d) { 
    num_report_ = num_d; 
}
/**
 * reads in contents of FASTA file into a map
 */
void DClassify::read_in_refs(std::string ref_fasta, 
			     std::map<std::string,std::string> &ref_h) {
    FastaRefID<FastaParser> fp(ref_fasta);
    fp.openFile();
    while(fp.hasNextSequence()) {
	std::pair<std::string,std::string> entry = fp.getNextEntry();
	std::string read_id = fp.getCurrIDLine().substr(1);
	std::string seq = entry.second;
	ref_h[read_id] = seq;
    }
}
/**
 * given a sequence and two intervals for V and J regions, returns vector 
 * of predicted D IDs
 */
std::vector<DLabel> DClassify::classify_d(std::string seq,
					       std::pair<int,int> v_part,
					       std::pair<int,int> j_part) {
    std::vector<DLabel> lst;
    // score the segment
    std::vector<DLabel> ret_v = score_d(seq, v_part, j_part);
    // if range overlapped, return empty list
    if(ret_v.empty()) { return lst; }
    // otherwise sort it based on score
    std::sort(ret_v.begin(), ret_v.end(), std::greater<DLabel>());
    
    DCLASS_DEBUG_PRINT(ret_v);
    // report the ids
    lst.resize(num_report_);
    for(int i = 0; i < num_report_; i++) {
	lst[i] = ret_v[i];
    }
    return lst;
}
/**
 * given a sequence and the intervals of V and J regions, returns a vector 
 * of pairs: <ID,score>
 */
std::vector<DLabel> 
DClassify::score_d(std::string seq, std::pair<int,int> v_part, std::pair<int,int> j_part) {
    std::vector<DLabel> ret_v;
    if(v_part.second >= j_part.first) { return ret_v; }
    
    int num_d = (int)drefs_h_.size();
    ret_v.reserve(num_d);
    
    std::string sub_seq = seq.substr(v_part.second, j_part.first-v_part.second+1);
    
    for(auto const &entry : drefs_h_) {
	std::string ref_id = entry.first;
	std::string ref_seq = entry.second;
        pair<int,int> score_p = get_local_score(ref_seq, sub_seq);
	DCLASS_DEBUG_PRINT(ref_id<<"\t"<<score_p);
	//
	DLabel dl(ref_id, score_p.first, score_p.second);
	ret_v.push_back( dl );
    }
    
    return ret_v;
}
/**
 * given a two sequences (one reference seq and one read substring), return the 
 * local alignment score
 */
pair<int,int> DClassify::get_local_score(std::string ref_seq, std::string read_substr) {
    TAlign ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), ref_seq);
    assignSource(row(ali, 1), read_substr);
    
    int score = localAlignment(ali, 
			       seqan::Score<int>(match_s_,mismatch_s_,indel_s_), 
			       seqan::SmithWaterman());
    //return score;
    // get num identity
    typedef seqan::Iterator<TRow>::Type TRowIterator;
    TRowIterator it = begin(row(ali, 0));
    TRowIterator itEnd = end(row(ali, 0));
    std::string s_row1 = "";
    for(; it != itEnd; ++it) {
	if(isGap(it)) s_row1 += "-";
	else s_row1 += (*it);
    }
    TRowIterator it2 = begin(row(ali, 1));
    TRowIterator itEnd2 = end(row(ali, 1));
    std::string s_row2 = "";
    for(; it2 != itEnd2; ++it2) {
	if(isGap(it2)) s_row2 += "-";
	else s_row2 += (*it2);
    }
    DCLASS_DEBUG_PRINT("\n"<<s_row1<<"\n"<<s_row2<<"\n");
    int a_len = (int)s_row1.size();
    int count = 0;
    for(int i = 0; i < a_len; i++) { 
	if(s_row1[i] == s_row2[i]) { count++; }
    }
    
    return pair<int,int>(score, count);
}
