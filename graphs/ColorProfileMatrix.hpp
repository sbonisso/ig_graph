#ifndef COLORPROFILEMATRIX_H_
#define COLORPROFILEMATRIX_H_

#include <vector>
#include <string>
#include <algorithm>

#include "prettyprint.hpp"
#include "Utilities.hpp"

#include "seq_utils/MutationNBModel.hpp"
#include "seq_utils/Encoding.hpp"

#include "graphs/CanonicalAntibodyGraph.hpp"

#ifdef DEBUG
#define COLORPROFMAT_DEBUG_PRINT(EXPR) DEBUG_PRINT("[COLOR_PROF_MAT]", EXPR)
#else
#define COLORPROFMAT_DEBUG_PRINT(x) do {} while (0)
#endif

using namespace std;

class ColorProfileMatrix {
public:
    // constructors
    ColorProfileMatrix();
    ColorProfileMatrix(int num_seq, int read_len, CanonicalAntibodyGraph *cab);
    
    // get methods
    int getReadLength() { return this->len_; }
    int getNumSeq() { return this->n_; }
    //
    vector<vector<int> > getVColorProfile() { return this->cp_mat_v_; }
    vector<vector<int> > getDColorProfile() { return this->cp_mat_d_; }
    vector<vector<int> > getJColorProfile() { return this->cp_mat_j_; }
    
    //vector<double> getLMerProbScore(int lmerLen, int leftShift, vector<vector<int> > cp_mat);
    vector<double> getLMerProbScore(vector<vector<int> > cp_mat);
    vector<double> getNormalizedLMerProbScoring(vector<vector<int> > cp_mat);
    //vector<double> getMemoizedProbScore(int lmerLen, int leftShift, int vStartIndex);
    
    // scoring methods
    //vector<double> getSimpleScoring(int start_ref_index, int end_ref_index);
    vector<int> getSimpleScoring(vector<vector<int> > cp_mat);
    /* vector<double> getSimpleMotifScoring(); */
    /* vector<double> getProbScoring(); */
    /* vector<double> getMutationMatrixScoring(); */
    
    // partition method
    vector<pair<int,int> > getPartitions(vector<vector<int> > cp_mat);
    vector<pair<int,int> > getPartitions(vector<vector<int> > cp_mat, bool fwd);
    
    //vector<vector<int> > cp_mat_;
    vector<vector<int> > cp_mat_v_;
    vector<vector<int> > cp_mat_d_;
    vector<vector<int> > cp_mat_j_;
    
private:
    int len_;
    int n_;
    CanonicalAntibodyGraph *cab_;
    
    void resize_mat(vector<vector<int> > &mat, int nrow, int ncol);
    
    vector<double> prob_scores_;
    vector<pair<int,double> > top_prob_ids_;
    
    vector<double> transformVScores(vector<double> &vscores);
    
    static const double THRESH_P;
    int lmer_len_;
    int lmer_shift_;
};


#endif
