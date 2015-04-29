#include "CreateProfile.h"

CreateProfile::CreateProfile() {
    n_ = 0;
}
/**
 *
 */
CreateProfile::CreateProfile(CanonicalAntibodyGraph *cab) {
    cab_ = cab;
    n_ = (*cab_).getNumReferences();
    v_scores_.resize( (*cab_).getNumV() );
    d_scores_.resize( (*cab_).getNumD() );
    j_scores_.resize( (*cab_).getNumJ() );
    max_report_ = 2;
}
CreateProfile::CreateProfile(CanonicalAntibodyGraph *cab, int max_n) : 
    CreateProfile(cab) {
    max_report_ = max_n;
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
    
    return cp;
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
vector<string> CreateProfile::getTopPredicted(int num_ret, 
					      vector<double> v,
					      Segment seg) {
    vector<string> refs(num_ret, "");
    vector<int> max_inds = getNMaxIndex(num_ret, v);
    //cout<<"MAX INDS:\t"<<max_inds<<endl;
    //
    for(int i = 0; i < num_ret; i++) {
	int ind = max_inds[i];
	if(seg == Segment::V_GENE) { refs[i] = (*cab_).getVRefID(ind); }
	else if(seg == Segment::D_GENE) { refs[i] = (*cab_).getDRefID(ind); }
	else if(seg == Segment::J_GENE) { refs[i] = (*cab_).getJRefID(ind); }
	else {}	
    }
    return refs;
}
/**
 * 
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
//vector<string> CreateProfile::getPredictedV(int num_ret) {
vector<string> CreateProfile::getPredictedV() {    
    return getTopPredicted(max_report_, v_scores_, Segment::V_GENE);
}
//vector<string> CreateProfile::getPredictedD(int num_ret) {
vector<string> CreateProfile::getPredictedD() {
    return getTopPredicted(max_report_, d_scores_, Segment::D_GENE);
}
//vector<string> CreateProfile::getPredictedJ(int num_ret) {
vector<string> CreateProfile::getPredictedJ() {
    return getTopPredicted(max_report_, j_scores_, Segment::J_GENE);
}
/**
 * return scores of predicted V, D, or J labels
 */
//vector<double> CreateProfile::getPredictedVScores(int num_ret) {    
vector<double> CreateProfile::getPredictedVScores() {
    int num_ret = max_report_;
    vector<double> n_scores(num_ret, -1);
    vector<int> max_inds = getNMaxIndex(num_ret, v_scores_);
    for(int i = 0; i < num_ret; i++) { 
	n_scores[i] = v_scores_[max_inds[i]];
    }
    return n_scores;
}
//vector<double> CreateProfile::getPredictedDScores(int num_ret) {
vector<double> CreateProfile::getPredictedDScores() {
    int num_ret = max_report_;
    vector<double> n_scores(num_ret, -1);
    vector<int> max_inds = getNMaxIndex(num_ret, d_scores_);
    for(int i = 0; i < num_ret; i++) { 
	n_scores[i] = d_scores_[max_inds[i]];
    }
    return n_scores;
}
//vector<double> CreateProfile::getPredictedJScores(int num_ret) {
vector<double> CreateProfile::getPredictedJScores() {
    int num_ret = max_report_;
    vector<double> n_scores(num_ret, -1);
    vector<int> max_inds = getNMaxIndex(num_ret, j_scores_);
    for(int i = 0; i < num_ret; i++) { 
	n_scores[i] = j_scores_[max_inds[i]];
    }
    return n_scores;
}
/**
 * overload operator<< 
 */
ostream& operator<< (ostream &out, CreateProfile &cp) {
    //
    vector<double> vscores = cp.getPredictedVScores();
    vector<double> dscores = cp.getPredictedDScores();
    vector<double> jscores = cp.getPredictedJScores();
    //
    vector<string> vpred = cp.getPredictedV();
    vector<string> dpred = cp.getPredictedD();
    vector<string> jpred = cp.getPredictedJ();	
    //	
    // out<<vpred<<"\t"<<dpred<<"\t"<<jpred<<"\t"
    //    <<vscores<<"\t"<<dscores<<"\t"<<jscores<<endl;
    int m = cp.max_report_;
    string delim = ",";
    for(int i = 0; i < m; i++) { out<<vpred[i]<<(i == m-1 ? "\t" : delim); }
    for(int i = 0; i < m; i++) { out<<dpred[i]<<(i == m-1 ? "\t" : delim); }
    //for(int i = 0; i < m; i++) { out<<jpred[i]<<(i == m-1 ? "\t" : delim); }
    for(int i = 0; i < m; i++) { out<<jpred[i]<<(i == m-1 ? "" : delim); }
    
    // for(int i = 0; i < m; i++) { out<<vscores[i]<<(i == m-1 ? "\t" : delim); }
    // for(int i = 0; i < m; i++) { out<<dscores[i]<<(i == m-1 ? "\t" : delim); }
    // for(int i = 0; i < m; i++) { out<<jscores[i]<<(i == m-1 ? "\t" : delim); }
    out<<endl;    
    //
    return out;
}
