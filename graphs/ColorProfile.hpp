/*
 * ColorProfile.h
 *
 *  Created on: Oct 31, 2013
 *      Author: sbonisso
 */

#ifndef COLORPROFILE_H_
#define COLORPROFILE_H_

#include <vector>
#include <string>
#include "Utils.h"
#include "prettyprint.hpp"

#include "seq_utils/MutationModel.h"
#include "seq_utils/MutationProbabilities.h"
#include "seq_utils/MutationNBProbabilities.h"
#include "seq_utils/MutationMatrix.h"
#include "seq_utils/MutationMemoizedProbabilities.h"

#ifdef DEBUG_COLORPROFILE
#define COLORPROFILE_DEBUG_PRINT(EXPR) DEBUG_PRINT("[COLORPROFILE]", EXPR)
#else
#define COLORPROFILE_DEBUG_PRINT(x) do {} while (0)
#endif

using namespace std;

class ColorProfile {
public:
    // constructors
    ColorProfile();
    ColorProfile(vector<string> colorP, vector<string> refS);
    ColorProfile(vector<string> colorP, vector<string> refS, vector<int> readIns);
	
    vector<pair<int,double> > topProbIDs;

    // get methods
    vector<string> getColorProfile() { return colorProfile; }
    vector<int> getReadIndex() { return readIndex; }
    vector<string> getRefSequences() { return refSeqs; }

    // set readIndex vector
    void setReadIndex(vector<int> rI) { readIndex = rI; }
    // set the read sequence, needed for some scoring methods
    void setReadSequence(string seq) { readSeq = seq; }

    // access to specific elements
    int getNthReadID(int i) { return readIndex[i]; }
    string getNthColorProfile(int i) { return colorProfile[i]; }
    string getNthRefProfile(int i) { return refSeqs[i]; }
    double getProbOfIndex(int i) { return (i < (int)probScores.size()) ? probScores[i] : -1.0; }

    // get total number of profiles
    int getNumProfiles() { return (int)colorProfile.size(); }

    vector<double> getLMerProbScore(int lmerLen, int leftShift);
    vector<double> getMemoizedProbScore(int lmerLen, int leftShift, int vStartIndex);
    vector<double> getNormalizedLMerProbScoring(int vStartIndex);

    // scoring methods
    vector<double> getSimpleScoring();
    vector<double> getSimpleMotifScoring();
    vector<double> getProbScoring();
    vector<double> getMutationMatrixScoring();
//	vector<double> getRFMutationScoring();

    // partition method
    vector<pair<int,int>> getPartitions();

private:
    vector<string> colorProfile;
    vector<int> readIndex;
    vector<string> refSeqs;

    vector<int> refPosStarts;
    vector<double> probScores;

    string readSeq;

    bool matchesWRCYMotif(string lmer);
    bool matchesRGYWMotif(string lmer);
    bool matchesWANMotif(string lmer);

    vector<double> transformVScores(vector<double> &vscores);

    static const double THRESH_P;
//	static const int LMER_LEN =   4; //5; //4; //8; //6; //4;
//	static const int LMER_SHIFT = 1; //2; //1; //3; //2; //1;
    int lmer_len;
    int lmer_shift;
};


#endif /* COLORPROFILE_H_ */
