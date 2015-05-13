#include "CreateProfile.h"

CreateProfile::CreateProfile() {
    n_ = 0;
    comp_cdr3_ = true;
}

CreateProfile::CreateProfile(CanonicalAntibodyGraph *cab) : 
    CreateProfile(cab, true)
{}
/**
 *
 */
CreateProfile::CreateProfile(CanonicalAntibodyGraph *cab, bool cmp_cdr3=true) {
    cab_ = cab;
    n_ = (*cab_).getNumReferences();
    v_scores_.resize( (*cab_).getNumV() );
    d_scores_.resize( (*cab_).getNumD() );
    j_scores_.resize( (*cab_).getNumJ() );
    max_report_ = 2;
    comp_cdr3_ = cmp_cdr3;
}
CreateProfile::CreateProfile(CanonicalAntibodyGraph *cab, 
			     int max_n, 
			     bool cmp_cdr3=true) : CreateProfile(cab, cmp_cdr3){
    max_report_ = max_n;
    comp_cdr3_ = cmp_cdr3;
}
/**
 *
 */
int CreateProfile::getProfileSize() { return n_; }
/**
 *
 */
ColorProfileMatrix CreateProfile::getColorProfile(string seq) {
    vector<int> vect(n_,0);
    int len = (int)seq.size();
    ColorProfileMatrix cp(n_, len, cab_);
    
    CREATEPROFILE_DEBUG_PRINT("NUM REF:\t"<<(*cab_).getNumReferences());
    // fill in initial 
    int index = 0;
    index = fillProfile((*cab_).getVRefs(), seq, index, Segment::V_GENE, &cp);
    index = fillProfile((*cab_).getDRefs(), seq, index, Segment::D_GENE, &cp);
    index = fillProfile((*cab_).getJRefs(), seq, index, Segment::J_GENE, &cp);
    //
    //v_scores_ = cp.getNormalizedLMerProbScoring(cp.cp_mat_v_);
    vector<int> vscores = cp.getSimpleScoring(cp.cp_mat_v_);
    vector<int> dscores = cp.getSimpleScoring(cp.cp_mat_d_);
    vector<int> jscores = cp.getSimpleScoring(cp.cp_mat_j_);
    //
    for(int i = 0; i < (int)vscores.size(); i++) { 
	v_scores_[i] = (double)vscores[i]; 
    }
    for(int i = 0; i < (int)dscores.size(); i++) { 
	d_scores_[i] = (double)dscores[i]; 
    }
    for(int i = 0; i < (int)jscores.size(); i++) { 
	j_scores_[i] = (double)jscores[i]; 
    }
    //    
    return cp;
}
/**
 *
 */
string CreateProfile::getCDR3() {
    return cdr3_str_;
}
/**
 *
 */
int CreateProfile::fillProfile(vector<string> ref_ids,
			       string seq,			       
			       int index,
			       Segment seg,
			       ColorProfileMatrix *cp) {
    CREATEPROFILE_DEBUG_PRINT("REF IDS:\t"<<ref_ids.size());
    for(int i = 0; i < (int)ref_ids.size(); i++) {
	string curr_ref_id = ref_ids[i];
	//
	if( seg == Segment::V_GENE ) {
	    cp->cp_mat_v_[i] = (*cab_).getVPainting(curr_ref_id, seq);
	}
	else if( seg == Segment::D_GENE ) {

	    cp->cp_mat_d_[i] = (*cab_).getDPainting(curr_ref_id, seq);
	}
	else if( seg == Segment::J_GENE ) {	    
	    cp->cp_mat_j_[i] = (*cab_).getJPainting(curr_ref_id, seq);
	}
	else {}
	//
	index++;
    }
    return index;
}
/**
 * performs the decision, selecting the top m (num_ret) scoring 
 * references from the given vector
 */
vector<string> CreateProfile::getTopPredicted_old(int num_ret, 
					      vector<double> v,
					      Segment seg) {
    vector<double> vals(num_ret, -1);
    vector<int> index(num_ret, -1);
    vector<string> refs(num_ret, "");
    for(int i = 0; i < num_ret; i++) {
	vector<double>::iterator itr = std::max_element(v.begin(), v.end());
	int ind = std::distance(v.begin(), itr);
	vals[i] = v[ind];
	index[i] = ind;
	v[ind] = -1;
	if(seg == Segment::V_GENE) { refs[i] = (*cab_).getVRefID(ind); }
	else if(seg == Segment::D_GENE) { refs[i] = (*cab_).getDRefID(ind); }
	else if(seg == Segment::J_GENE) { refs[i] = (*cab_).getJRefID(ind); }
	else {}	
    }
    return refs;
}
/**
 * return vector of labels given the vector of scores and the segment type
 */
vector<string> CreateProfile::getTopPredicted(int num_ret, 
					      vector<double> v,
					      Segment seg) {
    vector<string> refs(num_ret, "");
    vector<int> max_inds = getNMaxIndex(num_ret, v);    
    //cout<<"MAX INDS:\t"<<max_inds<<endl;
    //
    for(int i = 0; i < num_ret; i++) {
	int ind = max_inds[i];
	if(seg == Segment::V_GENE) { 
	    refs[i] = (*cab_).getVRefID(ind);
	    v_max_ind_ = max_inds;
	    v_pred_scores_ = getPredScores(num_ret, max_inds, v);
	}
	else if(seg == Segment::D_GENE) { 
	    if(v[ind] == 0) refs[i] = "?";
	    else refs[i] = (*cab_).getDRefID(ind);
	    d_max_ind_ = max_inds;
	    d_pred_scores_ = getPredScores(num_ret, max_inds, v);
	}
	else if(seg == Segment::J_GENE) { 
	    refs[i] = (*cab_).getJRefID(ind);
	    j_max_ind_ = max_inds;
	    j_pred_scores_ = getPredScores(num_ret, max_inds, v);
	}
	else {}	
    }
    return refs;
}
/**
 * return indeces of the top N (num_ret) that have maximal scores in vector v
 */
vector<int> CreateProfile::getNMaxIndex(int num_ret, vector<double> &v) {
    vector<double> vals(num_ret, -1);
    vector<int> index(num_ret, -1);
    // find first maximal element
    vector<double>::iterator itr = std::max_element(v.begin(), v.end());
    int ind = std::distance(v.begin(), itr);
    index[0] = ind;
    // now find remainder in naive fashon
    for(int i = 1; i < num_ret; i++) {
	int max_val = -1;
	int max_index = -1;
	for(int j = 0; j < (int)v.size(); j++) {
	    if(std::find(index.begin(), index.end(), j) == index.end()) {
		if(max_val < v[j]) {
		    max_val = v[j];
		    max_index = j;
		}
	    }
	}
	index[i] = max_index;
    }
    return index;
}
/**
 * return predicted V, D, or J labels
 */
vector<string> CreateProfile::getPredictedV() {    
    //return getTopPredicted(max_report_, v_scores_, Segment::V_GENE);
    return pred_v_;
}
vector<string> CreateProfile::getPredictedD() {
    //return getTopPredicted(max_report_, d_scores_, Segment::D_GENE);
    return pred_d_;
}
vector<string> CreateProfile::getPredictedJ() {
    //return getTopPredicted(max_report_, j_scores_, Segment::J_GENE);
    return pred_j_;
}
/**
 * return scores of predicted V, D, or J labels
 */
vector<double> CreateProfile::getPredictedVScores() {
    // int num_ret = max_report_;
    // vector<double> n_scores(num_ret, -1);
    // vector<int> max_inds = getNMaxIndex(num_ret, v_scores_);
    // v_max_ind_ = max_inds;
    // for(int i = 0; i < num_ret; i++) { 
    // 	n_scores[i] = v_scores_[max_inds[i]];
    // }
    // return n_scores;
    return v_pred_scores_;
}
vector<double> CreateProfile::getPredictedDScores() {
    // int num_ret = max_report_;
    // vector<double> n_scores(num_ret, -1);
    // vector<int> max_inds = getNMaxIndex(num_ret, d_scores_);
    // d_max_ind_ = max_inds;
    // for(int i = 0; i < num_ret; i++) { 
    // 	n_scores[i] = d_scores_[max_inds[i]];
    // }
    // return n_scores;
    return d_pred_scores_;
}
vector<double> CreateProfile::getPredictedJScores() {
    // int num_ret = max_report_;
    // vector<double> n_scores(num_ret, -1);
    // vector<int> max_inds = getNMaxIndex(num_ret, j_scores_);    
    // j_max_ind_ = max_inds;
    // for(int i = 0; i < num_ret; i++) { 
    // 	n_scores[i] = j_scores_[max_inds[i]];
    // }
    // return n_scores;
    return j_pred_scores_;
}
/**
 * given a vector of scores (score_v) and a ranking of indeces (max_n_ind)
 * return a vector of the scores of the top predicted
 */
vector<double> CreateProfile::getPredScores(int num_ret,
					    vector<int> max_inds, 
					    vector<double> score_v) {
    vector<double> n_scores(num_ret, -1);
    for(int i = 0; i < num_ret; i++) { 
	n_scores[i] = j_scores_[max_inds[i]];
    }
    return n_scores;
}

void CreateProfile::compute(string seq) {
    ColorProfileMatrix cp_mat = this->getColorProfile(seq);
    
    this->pred_v_ = getTopPredicted(max_report_, v_scores_, Segment::V_GENE);
    this->pred_d_ = getTopPredicted(max_report_, d_scores_, Segment::D_GENE);
    this->pred_j_ = getTopPredicted(max_report_, j_scores_, Segment::J_GENE);
    
    // if couldn't find cdr3 (possibly from truncated sequence) return nothing
    if(comp_cdr3_ && !this->computeCDR3(seq, cp_mat)) {
	cdr3_str_ = "?";
    }

}

/**
 * 
 */
//string CreateProfile::computeCDR3(string seq, ColorProfileMatrix &cp_mat) {
bool CreateProfile::computeCDR3(string seq, ColorProfileMatrix &cp_mat) {
    //
    seqan::String<seqan::Dna> gSeq;
    seqan::assign(gSeq, seq);
    seqan::StringSet<seqan::String<seqan::AminoAcid> > tSeqs;
    seqan::translate(tSeqs, gSeq, seqan::WITH_FRAME_SHIFTS);
    //
    int best_frame = 0; 
    seqan::String<seqan::AminoAcid> frame_seq;
    int num_stop_bf = 100;
    for(int i = 0; i < (int)seqan::length(tSeqs); i++) {
	int num_stop = 0;
	for(int j = 0; j < (int)seqan::length(tSeqs[i]); j++) {
	    if(tSeqs[i][j] == '*') { num_stop++; }
	}
	if(num_stop < num_stop_bf) { 
	    best_frame = i;
	    num_stop_bf = num_stop;
	    frame_seq = tSeqs[i];
	}
    }
    //   
    vector<pair<int,int> > v_part = cp_mat.getPartitions(cp_mat.getVColorProfile());
    vector<pair<int,int> > d_part = cp_mat.getPartitions(cp_mat.getDColorProfile());
    vector<pair<int,int> > j_part = cp_mat.getPartitions(cp_mat.getJColorProfile());
    //v_range = cp_mat
    // cout<<v_part[v_max_ind_[0]]<<endl;
    // cout<<d_part[d_max_ind_[0]]<<endl;
    // cout<<j_part[j_max_ind_[0]]<<endl;

    //
    pair<int,int> v_range = v_part[v_max_ind_[0]];
    pair<int,int> j_range = j_part[j_max_ind_[0]];    
    if(v_range.first > v_range.second) { return false; }
    if(j_range.first > j_range.second) { return false; }
    
    // find 2nd cys - start codon after
    int cdr3_start = 0;
    for(int i = (v_range.second+best_frame)/3; i > 0; i--) {
	if(frame_seq[i] == 'C') { 
	    cdr3_start = ((i+1)*3 - best_frame);
	    break;
	}
    }
    int cdr3_end = 0;
    // take 5bp from start of J
    cdr3_end = j_range.first+5;
    
    // cout<<seq<<endl;
    // cout<<best_frame<<endl;
    // cout<<frame_seq<<endl;
    // cout<<v_range<<endl<<j_range<<endl;
    // cout<<(cdr3_start+best_frame)/3<<"\t"<<(cdr3_end+best_frame)/3<<endl;
    // cout<<cdr3_start<<"\t"<<cdr3_end<<endl;
    
    //return seq.substr(cdr3_start, cdr3_end-cdr3_start);
    cdr3_str_ = seq.substr(cdr3_start, cdr3_end-cdr3_start);
    return true;
}
/**
 * overload operator<< 
 */
ostream& operator<< (ostream &out, CreateProfile &cp) {
    //
    vector<string> vpred = cp.getPredictedV();
    vector<string> dpred = cp.getPredictedD();
    vector<string> jpred = cp.getPredictedJ();
    //
    vector<double> vscores = cp.getPredictedVScores();
    vector<double> dscores = cp.getPredictedDScores();
    vector<double> jscores = cp.getPredictedJScores();
    // 
    // out<<vpred<<"\t"<<dpred<<"\t"<<jpred<<"\t"
    //    <<vscores<<"\t"<<dscores<<"\t"<<jscores<<endl;
    int m = cp.max_report_;
    string delim = ",";
    for(int i = 0; i < m; i++) { out<<vpred[i]<<(i == m-1 ? "\t" : delim); }
    for(int i = 0; i < m; i++) { out<<dpred[i]<<(i == m-1 ? "\t" : delim); }

    // for(int i = 0; i < m; i++) { out<<jpred[i]<<(i == m-1 ? "\t" : delim); }
    // out<<cp.getCDR3(); 

    // if compute CDR3, leave space for col
    if(cp.comp_cdr3_) { 
    	for(int i = 0; i < m; i++) { out<<jpred[i]<<(i == m-1 ? "\t" : delim); }
    	out<<cp.getCDR3(); 
    }
    else {
    	for(int i = 0; i < m; i++) { out<<jpred[i]<<(i == m-1 ? "" : delim); }
    }
    
    // for(int i = 0; i < m; i++) { out<<vscores[i]<<(i == m-1 ? "\t" : delim); }
    // for(int i = 0; i < m; i++) { out<<dscores[i]<<(i == m-1 ? "\t" : delim); }
    // for(int i = 0; i < m; i++) { out<<jscores[i]<<(i == m-1 ? "\t" : delim); }
    out<<endl;    
    //
    return out;
}
