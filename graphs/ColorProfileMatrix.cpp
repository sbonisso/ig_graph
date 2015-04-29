#include "graphs/ColorProfileMatrix.h"

const double ColorProfileMatrix::THRESH_P = 0.8;  //0.999; //0.80; // 0.95

ColorProfileMatrix::ColorProfileMatrix() {
    this->n_ = 0;
    this->len_ = 0;
    this->cab_ = NULL;
    //
    this->lmer_len_ = (MutationNBProbabilities::getInstance()).getLMerLen();
    this->lmer_shift_ = (lmer_len_ == 6) ? 2 : lmer_len_-3;
}

ColorProfileMatrix::ColorProfileMatrix(int num_seq, 
				       int read_len, 
				       CanonicalAntibodyGraph *cab) {
    this->n_ = num_seq;
    this->len_ = read_len;
    this->cab_ = cab;
    //
    resize_mat(cp_mat_v_, (*cab_).getNumV(), this->len_);
    resize_mat(cp_mat_d_, (*cab_).getNumD(), this->len_);
    resize_mat(cp_mat_j_, (*cab_).getNumJ(), this->len_);
    //
    this->lmer_len_ = (MutationNBProbabilities::getInstance()).getLMerLen();
    this->lmer_shift_ = (lmer_len_ == 6) ? 2 : lmer_len_-3;
}
/**
 *
 */
void ColorProfileMatrix::resize_mat(vector<vector<int> > &mat, int nrow, int len) {
    mat.resize(nrow);
    for(int i = 0; i < nrow; i++) {
	mat[i].resize(len,-1);
    }
}
/**
 * generalized mutation prob lookup, given an lmer length and a shifting size
 * @param lmerLen length of l-mer used
 * @param leftShift length of the shift to be used in selecting the 
 *        proper l-mer from the sequence
 */
vector<double> ColorProfileMatrix::getLMerProbScore(//int lmerLen, 
						    //int leftShift,
						    vector<vector<int> > cp_mat){
    int lmerLen = lmer_len_;
    int leftShift = lmer_shift_;
    string nullLmerStr(lmerLen,'.');
    //int start_ref_index = 0;
    //int end_ref_index = (*cab_).getNumV()-1;
    //int num_ref = end_ref_index - start_ref_index + 1;    
    
    int n = (int)cp_mat.size();
    prob_scores_.resize(n, 0);
    vector<double> vscores(n, 0);
    if(n_ == 0) return vscores;
    // for each profile...
    for(int i = 0; i < n; i++) {
    //for(int i = 0; i < end_ref_index; i++) {
	// get start pos of read on reference
	//int posPad = refPosStarts[i];	
	string ref_id = (*cab_).getVRefID(i);
	string curr_ref_seq = (*cab_).getVRefSeq(ref_id);
	// compute by summing matches
	for(int j = 10; j < len_ && j < 300; j++) {
	    if(j+lmerLen >= (int)curr_ref_seq.size()) { 
		vscores[i] += log2(1.0);
		continue;
	    }
	    
	    string lmer( curr_ref_seq.substr(j-leftShift, lmerLen) );
	    
	    double prob = 0.0;
	    // if find mismatch, use mismatch l-mer
	    if(std::find(lmer.begin(), lmer.end(), '.') != lmer.end()) {
		prob = (MutationNBProbabilities::getInstance()).
		    getProb(nullLmerStr, j);
		    //getProb(nullLmerStr, posPad+j);
	    }
	    // otherwise, use l-mer of reference
	    else {
		double p = (MutationNBProbabilities::getInstance()).
		    getProb(lmer, j);
		    //getProb(lmer, posPad+j);
		if(cp_mat[i][j] == -1) { prob = p; }
		else { prob = (1.0-p); }
	    }
	    vscores[i] += log2(prob);
	    
	    COLORPROFMAT_DEBUG_PRINT(i<<"\t"<<j<<"\t"<<j<<"\t"//posPad+j<<"\t"
				     <<lmer<<"\t"<<cp_mat[i][j]<<"\t"
				     <<prob<<"\t"<<log2(prob)<<"\t"
				     <<vscores[i]);
	}
    }
    return vscores;
}
/**
 *
 */
vector<double> ColorProfileMatrix::getNormalizedLMerProbScoring(vector<vector<int> > cp_mat) {
    //vector<double> vscores = this->getLMerProbScore(lmer_len_, lmer_shift_);
    vector<double> vscores = this->getLMerProbScore(cp_mat);
    return this->transformVScores(vscores);
}
/**
 *
 */
vector<double> ColorProfileMatrix::transformVScores(vector<double> &vscores) {
    COLORPROFMAT_DEBUG_PRINT("VSCORES:\t"<<vscores);
    double probMut = -numeric_limits<double>::infinity();
    if(vscores.size() == 0) { return vscores; }      // no matching, pass empty

    COLORPROFMAT_DEBUG_PRINT("PROB MUT:\t"<<probMut);
    for(int i = 0; i < (int)vscores.size(); i++) { 
	if(vscores[i] == 0) { vscores[i] = probMut; } 
    }

    // scale/transform back to probabilities
    double maxProb = *std::max_element(vscores.begin(), vscores.end());
    for(int i = 0; i < (int)vscores.size(); i++) { 
	vscores[i] -= maxProb;          // subtract max from all
    }    
    double total = 0.0;
    for(int i = 0; i < (int)vscores.size(); i++) { 
	total += pow(2.0, vscores[i]);  // compute denom
    }
    for(int i = 0; i < (int)vscores.size(); i++) { 
	vscores[i] = pow(2.0,vscores[i]) / total;  // compute prob using diff
    }
    for(int i = 0; i < (int)vscores.size(); i++) { 
	prob_scores_[i] = vscores[i]; 
    }
    
    // now rank-score (in reverse, i.e., highest prob -> highest rank) 
    // such that can use standard score picking scheme
    vector<pair<int,double>> vscoreIndex(vscores.size());
    for(int i = 0; i < (int)vscores.size(); i++) { 
	vscoreIndex[i] = pair<int,double>(i, vscores[i]); 
    }
    std::sort(vscoreIndex.begin(), 
	      vscoreIndex.end(), 
	      [](const pair<int,double> a, const pair<int,double> b)
	      {return a.second > b.second; }
	);
    // the top percent will tie for highest rank, all others rank in order, ignoring ties
    double cutOff = ColorProfileMatrix::THRESH_P;
    double cumProb = 0.0;
    for(int i = 0; i < (int)vscoreIndex.size(); i++) {
	double currProb = vscoreIndex[i].second;
	int origIndex = vscoreIndex[i].first;
	if(cumProb < cutOff) { vscores[origIndex] = (int)vscores.size(); }
	else { vscores[origIndex] = (int)vscores.size() - i; }
	// for development purposes ...
	if(cumProb < cutOff) { 
	    top_prob_ids_.push_back( pair<int,double>(origIndex, currProb) ); 
	}
	cumProb += currProb;
    }
    
    COLORPROFMAT_DEBUG_PRINT("SORT STUFF");
    for(int i = 0; i < (int)vscoreIndex.size(); i++) { 
	COLORPROFMAT_DEBUG_PRINT(i<<"\t"<<vscoreIndex[i].first
				       <<"\t"<<vscoreIndex[i].second); 
    }
    COLORPROFMAT_DEBUG_PRINT("PROB SCORE:\t"<<vscores);
    
    return vscores;
}
/**
 * simplest scoring scheme, +1 for a match, 0 otherwise
 * @returns vector of integer scores for each profile in class
 */
vector<int> ColorProfileMatrix::getSimpleScoring(vector<vector<int> > cp_mat) {
    int n = (int)cp_mat.size();
    //int curr_len = len_ - (*cab_).getK();
    int curr_len = len_;
    vector<int> vscores(n, 0);
    // for each profile...
    for(int i = 0; i < n; i++) {
	vscores[i] = 0;
	// compute by summing matches
	for(int j = 0; j < curr_len; j++) {
	    if(cp_mat[i][j] >= 0) vscores[i] += 1;
	}
    }
    return vscores;
}
