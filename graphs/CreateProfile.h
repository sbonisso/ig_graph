#ifndef CREATEPROFILE_H_
#define CREATEPROFILE_H_

#include <string>
#include <vector>
#include <ctime>
#include <time.h>
#include <unordered_map>

#include "Utils.h"

#include "ReferenceMap.h"
#include "CanonicalAntibodyGraph.h"
#include "ColorProfileMatrix.h"

#include "unit_tests/TestCreateProfile.h"

//#ifdef DEBUG_CREATEPROFILE
#ifdef DEBUG
#define CREATEPROFILE_DEBUG_PRINT(EXPR) DEBUG_PRINT("[CREATE_PROFILE]", EXPR)
#else
#define CREATEPROFILE_DEBUG_PRINT(x) do {} while (0)
#endif


using namespace std;

class CreateProfile {
public:
    CreateProfile();
    CreateProfile(CanonicalAntibodyGraph *cab);
    CreateProfile(CanonicalAntibodyGraph *cab, int max_n);
    
    ColorProfileMatrix getColorProfile(string seq);

    int getProfileSize();
    //
    vector<double> getVScores() { return v_scores_; }
    vector<double> getDScores() { return d_scores_; }
    vector<double> getJScores() { return j_scores_; }
    //
    /* vector<string> getPredictedV(int num_ret); */
    /* vector<string> getPredictedD(int num_ret); */
    /* vector<string> getPredictedJ(int num_ret); */
    vector<string> getPredictedV();
    vector<string> getPredictedD();
    vector<string> getPredictedJ();
    //
    /* vector<double> getPredictedVScores(int num_ret); */
    /* vector<double> getPredictedDScores(int num_ret); */
    /* vector<double> getPredictedJScores(int num_ret); */
    vector<double> getPredictedVScores();
    vector<double> getPredictedDScores();
    vector<double> getPredictedJScores();
    //
    vector<string> getVRefVect() { return (*cab_).getVRefs(); }
    vector<string> getDRefVect() { return (*cab_).getDRefs(); }
    vector<string> getJRefVect() { return (*cab_).getJRefs(); }
    
    friend ostream& operator<< (ostream &out, CreateProfile &cp);
    
    friend class TestCreateProfile; // for unit testing

private:
    void init();

    int fillProfile(vector<string> ref_ids,
		    string seq,			       
		    int index,
		    Segment seg,
		    ColorProfileMatrix *cp);
    
    vector<string> getTopPredicted_old(int num_ret, vector<double> v, Segment seg);
    vector<string> getTopPredicted(int num_ret, vector<double> v, Segment seg);

    vector<int> getNMaxIndex(int num_red, vector<double> &v);

    int n_;
    int max_report_;

    CanonicalAntibodyGraph *cab_;
    vector<double> v_scores_;
    vector<double> d_scores_;
    vector<double> j_scores_;
};

#endif
