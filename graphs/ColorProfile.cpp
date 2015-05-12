/*
 * ColorProfile.cpp
 *
 *  Created on: Oct 31, 2013
 *      Author: sbonisso
 */

#include "ColorProfile.h"

// init class static constant
const double ColorProfile::THRESH_P = 0.8;  //0.999; //0.80;	// 0.95

ColorProfile::ColorProfile() {
    lmer_len = (MutationNBProbabilities::getInstance()).getLMerLen();
    lmer_shift = (lmer_len == 6) ? 2 : lmer_len-3;
}

ColorProfile::ColorProfile(vector<string> colorP, vector<string> refS) {
    colorProfile = colorP;
    refSeqs = refS;
    refPosStarts.resize(refSeqs.size(), 125);

    lmer_len = (MutationNBProbabilities::getInstance()).getLMerLen();
    lmer_shift = (lmer_len == 6) ? 2 : lmer_len-3;
}

//ColorProfile::ColorProfile(vector<string> colorP, vector<string> refS, vector<int> readIns) {
//	colorProfile = colorP;
//	refSeqs = refS;
//	readIndex = readIns;
//}
ColorProfile::ColorProfile(vector<string> colorP, vector<string> refS, vector<int> refPos) {
    colorProfile = colorP;
    refSeqs = refS;
    refPosStarts = refPos;

    lmer_len = (MutationNBProbabilities::getInstance()).getLMerLen();
    lmer_shift = (lmer_len == 6) ? 2 : lmer_len-3;
}
/**
 * scoring scheme that uses mutation matrix (log-odds scores of l-mer mutations) to asses quality
 */
vector<double> ColorProfile::getMutationMatrixScoring() {
    COLORPROFILE_DEBUG_PRINT("IN MUTATION MATRIX SCORING");
    COLORPROFILE_DEBUG_PRINT("num = "<<colorProfile.size());
    COLORPROFILE_DEBUG_PRINT("SEQ:\t"<<readSeq);

    probScores.resize(colorProfile.size(), 0); 	      // hack for unit test only
    vector<double> vscores(colorProfile.size(), 0);
    // for each profile...
    for(int i = 0; i < (int)colorProfile.size(); i++) {
	// compute by summing matches
	COLORPROFILE_DEBUG_PRINT("REF SEQ "<<i<<"\t"<<refSeqs[i]);
	for(int j = 1; j < (int)colorProfile[i].size()-4; j++) {
//			if(colorProfile[i][j] != '.') vscores[i] += 1;

	    string readLmer = readSeq.substr(j-1,4);	// should never find '.' in read sequence
	    string refLmer = (refSeqs[i]).substr(j-1,4);

	    // found a mismatch in the reference, means large deviation, more than a few mutations
	    if(std::find(refLmer.begin(), refLmer.end(), '.') != refLmer.end()) {
		vscores[i] += 0;
		COLORPROFILE_DEBUG_PRINT("MUTATION MATRIX\t"<<j<<"\t"<<refLmer<<"\t"<<readLmer);
	    }
	    else {

		double val = (MutationMatrix::getInstance()).getScore(refLmer, readLmer);
		vscores[i] += val;
		COLORPROFILE_DEBUG_PRINT("MUTATION MATRIX\t"<<j<<"\t"<<refLmer<<"\t"<<readLmer<<"\t"<<val<<"\t"<<vscores[i]);
	    }

	}
	probScores[i] = vscores[i];		// hack for unit test only
    }
    return vscores;
}
/**
 * scoring scheme revolving around probabilistic interpretation of matches/mismatches
 */
vector<double> ColorProfile::getProbScoring() {
//	int posPad = 125;

    probScores.resize(colorProfile.size(), 0);
    vector<double> vscores(colorProfile.size(), 0);
    if(colorProfile.size() == 0) return vscores;

    // for each profile...
    for(int i = 0; i < (int)colorProfile.size(); i++) {
	// get start pos of read on reference
	int posPad = refPosStarts[i];
//		VERBOSE_PRINT("PAD START:\t"<<posPad<<"\t"<<i);
	// compute by summing matches
//		for(int j = 0; j < (int)colorProfile[i].size(); j++) {
//		for(int j = 1; j < (int)colorProfile[i].size() && j < 174; j++) {
//		for(int j = 10; j < (int)colorProfile[i].size() && ((posPad+j) < 300); j++) {
	for(int j = 10; j < (int)colorProfile[i].size() && ((posPad+j) < 300); j++) {
	    string lmer(refSeqs[i].substr(j-1, 4));
//			string lmer(refSeqs[i].substr(j-3,8));

	    double prob = 0.0;
	    // if find mismatch, use mismatch l-mer
	    if(std::find(lmer.begin(), lmer.end(), '.') != lmer.end()) {
		prob = (MutationProbabilities::getInstance()).getProb("....", posPad+j);
//				prob = 0.05;
	    }
	    else {		// otherwise, use l-mer of reference
		double p = (MutationProbabilities::getInstance()).getProb(lmer, posPad+j);
//				double p = (MutationModel::getInstance()).getProb(lmer, posPad+j);
		if(colorProfile[i][j] == '.') { prob = p; }
		else { prob = (1.0-p); }
//				double decThresh = 0.75;
//				if(colorProfile[i][j] == '.' && p < decThresh) { prob = 0.1; }
//				else if(colorProfile[i][j] != '.' && p < decThresh) { prob = 0.9; }
//				else if(colorProfile[i][j] == '.' && p >= decThresh) { prob = 0.99; }
//				else if(colorProfile[i][j] != '.' && p >= decThresh) { prob = 0.01; }
//				else { prob = 1.0; }
	    }
	    vscores[i] += log2(prob);

//			if(std::find(lmer.begin(), lmer.end(), '.') != lmer.end()) continue; 		// skip if found a '.' in l-mer
//
//			double prob = (MutationProbabilities::getInstance()).getProb(lmer, posPad+j);
//
//			if(colorProfile[i][j] == '.') { vscores[i] += prob; }
//			else { vscores[i] += (1.0-prob); }

	    COLORPROFILE_DEBUG_PRINT(i<<"\t"<<j<<"\t"<<posPad+j<<"\t"<<lmer<<"\t"<<colorProfile[i][j]<<"\t"<<prob<<"\t"<<log2(prob)<<"\t"<<vscores[i]);

//			if(colorProfile[i][j] != '.') { vscores[i] += 1; }
	}
    }
    COLORPROFILE_DEBUG_PRINT("VSCORES:\t"<<vscores);
    double probMut = -numeric_limits<double>::infinity();
//	for(int i = 0; i < (int)colorProfile[0].size(); i++) {
//		probMut += log2( ((MutationProbabilities::getInstance()).getProb("....", refPosStarts[0]+i)) ); }
    COLORPROFILE_DEBUG_PRINT("PROB MUT:\t"<<probMut);
    for(int i = 0; i < (int)vscores.size(); i++) { if(vscores[i] == 0) { vscores[i] = probMut; } }

    // scale/transform back to probabilities
    double maxProb = *std::max_element(vscores.begin(), vscores.end());							// find max
    for(int i = 0; i < (int)vscores.size(); i++) { vscores[i] -= maxProb; }						// subtract max from all
    double total = 0.0;
    for(int i = 0; i < (int)vscores.size(); i++) { total += pow(2.0, vscores[i]); }				// compute denom
    for(int i = 0; i < (int)vscores.size(); i++) { vscores[i] = pow(2.0,vscores[i]) / total; }	// compute prob using difference
    for(int i = 0; i < (int)vscores.size(); i++) { probScores[i] = vscores[i]; }

    // now rank-score (in reverse, i.e., highest prob -> highest rank) such that can use standard score picking scheme
    vector<pair<int,double>> vscoreIndex(vscores.size());
    for(int i = 0; i < (int)vscores.size(); i++) { vscoreIndex[i] = pair<int,double>(i, vscores[i]); }
    std::sort(vscoreIndex.begin(), vscoreIndex.end(), [](const pair<int,double> a, const pair<int,double> b){ return a.second > b.second; });
    // the top percent will tie for highest rank, all others rank in order, ignoring ties
    double cutOff = ColorProfile::THRESH_P; //0.801; //0.999; //0.80;	// 0.95
    double cumProb = 0.0;
    for(int i = 0; i < (int)vscoreIndex.size(); i++) {
	double currProb = vscoreIndex[i].second;
	int origIndex = vscoreIndex[i].first;
	if(cumProb < cutOff) { vscores[origIndex] = (int)vscores.size(); }
	else { vscores[origIndex] = (int)vscores.size() - i; }
	// for development purposes ...
	if(cumProb < cutOff) { topProbIDs.push_back( pair<int,double>(origIndex, currProb) ); }

	cumProb += currProb;
    }


    COLOR_DEBUG_PRINT("SORT STUFF");
    for(int i = 0; i < (int)vscoreIndex.size(); i++) { COLOR_DEBUG_PRINT(i<<"\t"<<vscoreIndex[i].first<<"\t"<<vscoreIndex[i].second); }

    COLOR_DEBUG_PRINT("PROB SCORE:\t"<<vscores);
    return vscores;
}
/**
 * generalized mutation probability lookup, given an lmer length and a shifting size
 * @param lmerLen length of l-mer used
 * @param leftShift length of the shift to be used in selecting the proper l-mer from the sequence
 */
vector<double> ColorProfile::getLMerProbScore(int lmerLen, int leftShift) {
    probScores.resize(colorProfile.size(), 0);
    vector<double> vscores(colorProfile.size(), 0);
    if(colorProfile.size() == 0) return vscores;
    string nullLmerStr(lmerLen,'.');

    // for each profile...
    for(int i = 0; i < (int)colorProfile.size(); i++) {
	// get start pos of read on reference
	int posPad = refPosStarts[i];
	// compute by summing matches
	for(int j = 10; j < (int)colorProfile[i].size() && ((posPad+j) < 300); j++) {
	    //string lmer(refSeqs[i].substr(j-1, 4));
	    string lmer(refSeqs[i].substr(j-leftShift, lmerLen));

	    double prob = 0.0;
	    // if find mismatch, use mismatch l-mer
	    if(std::find(lmer.begin(), lmer.end(), '.') != lmer.end()) {
//				prob = (MutationProbabilities::getInstance()).getProb(nullLmerStr, posPad+j);
		prob = (MutationNBProbabilities::getInstance()).getProb(nullLmerStr, posPad+j);
	    }
	    else {		// otherwise, use l-mer of reference
//				double p = (MutationProbabilities::getInstance()).getProb(lmer, posPad+j);
		double p = (MutationNBProbabilities::getInstance()).getProb(lmer, posPad+j);
		if(colorProfile[i][j] == '.') { prob = p; }
		else { prob = (1.0-p); }
	    }
	    vscores[i] += log2(prob);

	    COLORPROFILE_DEBUG_PRINT(i<<"\t"<<j<<"\t"<<posPad+j<<"\t"<<lmer<<"\t"<<colorProfile[i][j]<<"\t"<<prob<<"\t"<<log2(prob)<<"\t"<<vscores[i]);
	}
    }
    return vscores;
}

vector<double> ColorProfile::getNormalizedLMerProbScoring(int vStartIndex) {
//	vector<double> vscores = this->getLMerProbScore(LMER_LEN, LMER_SHIFT);
    vector<double> vscores = this->getLMerProbScore(lmer_len, lmer_shift);
//	vector<double> vscores = this->getMemoizedProbScore(LMER_LEN, LMER_SHIFT, vStartIndex);
    return this->transformVScores(vscores);
}

vector<double> ColorProfile::getMemoizedProbScore(int lmerLen, int leftShift, int vStartIndex) {
    probScores.resize(colorProfile.size(), 0);
    vector<double> vscores(colorProfile.size(), 0);
    if(colorProfile.size() == 0) return vscores;
    string nullLmerStr(lmerLen,'.');

    // for each profile...
    for(int i = 0; i < (int)colorProfile.size(); i++) {
	int currColorIndex = this->readIndex[i]-vStartIndex;
	// get start pos of read on reference
	int posPad = refPosStarts[i];
	// compute by summing matches
	for(int j = 10; j < (int)colorProfile[i].size() && ((posPad+j) < 300); j++) {
	    string lmer(refSeqs[i].substr(j-leftShift, lmerLen));

	    double prob = 0.0;
//			double cmpProb = 0.0;
	    // if find mismatch, use mismatch l-mer
	    if(std::find(lmer.begin(), lmer.end(), '.') != lmer.end()) {
		prob = (MutationMemoizedProbabilities::getInstance()).getProbNull(posPad+j);
	    }
	    else {		// otherwise, use l-mer of reference
//				double pOld = (MutationNBProbabilities::getInstance()).getProb(lmer, posPad+j);

		double p = (MutationMemoizedProbabilities::getInstance()).getProb(currColorIndex, posPad+j);
		if(colorProfile[i][j] == '.') { prob = p; } //cmpProb = pOld;}
		else { prob = (1.0-p); }  //cmpProb = (1.0-pOld); }
	    }
	    vscores[i] += log2(prob);

	    COLORPROFILE_DEBUG_PRINT(i<<"\t"<<j<<"\t"<<posPad+j<<"\t"<<lmer<<"\t"<<colorProfile[i][j]<<"\t"<<prob<<"\t"<<log2(prob)<<"\t"<<vscores[i]);
	    //COLORPROFILE_DEBUG_PRINT(readIndex[i]<<":\t"<<i<<"\t"<<j<<"\t"<<posPad+j<<"\t"<<lmer<<"\t"<<colorProfile[i][j]<<"\t"<<prob<<"\t"<<cmpProb<<"\t"<<log2(prob)<<"\t"<<vscores[i]);
	    //COLORPROFILE_DEBUG_PRINT(readIndex[i]<<":\t"<<i<<"\t"<<j<<"\t"<<posPad+j<<"\t"<<lmer<<"\t"<<colorProfile[i][j]<<"\t"<<prob<<"\t"<<cmpProb<<"\t"<<abs(prob-cmpProb)<<"\t"<<log2(prob)<<"\t"<<vscores[i]);
	}
    }
    return vscores;
}
/**
 *
 */
vector<double> ColorProfile::transformVScores(vector<double> &vscores) {
    COLORPROFILE_DEBUG_PRINT("VSCORES:\t"<<vscores);
    double probMut = -numeric_limits<double>::infinity();
    if(vscores.size() == 0) { return vscores; } 	// no matching, pass empty

    COLORPROFILE_DEBUG_PRINT("PROB MUT:\t"<<probMut);
    for(int i = 0; i < (int)vscores.size(); i++) { if(vscores[i] == 0) { vscores[i] = probMut; } }

    // scale/transform back to probabilities
    double maxProb = *std::max_element(vscores.begin(), vscores.end());							// find max
    for(int i = 0; i < (int)vscores.size(); i++) { vscores[i] -= maxProb; }						// subtract max from all
    double total = 0.0;
    for(int i = 0; i < (int)vscores.size(); i++) { total += pow(2.0, vscores[i]); }				// compute denom
    for(int i = 0; i < (int)vscores.size(); i++) { vscores[i] = pow(2.0,vscores[i]) / total; }	// compute prob using difference
    for(int i = 0; i < (int)vscores.size(); i++) { probScores[i] = vscores[i]; }

    // now rank-score (in reverse, i.e., highest prob -> highest rank) such that can use standard score picking scheme
    vector<pair<int,double>> vscoreIndex(vscores.size());
    for(int i = 0; i < (int)vscores.size(); i++) { vscoreIndex[i] = pair<int,double>(i, vscores[i]); }
    std::sort(vscoreIndex.begin(), vscoreIndex.end(), [](const pair<int,double> a, const pair<int,double> b){ return a.second > b.second; });
    // the top percent will tie for highest rank, all others rank in order, ignoring ties
    double cutOff = ColorProfile::THRESH_P;  //0.80; //0.999; //0.80;	// 0.95
    double cumProb = 0.0;
    for(int i = 0; i < (int)vscoreIndex.size(); i++) {
	double currProb = vscoreIndex[i].second;
	int origIndex = vscoreIndex[i].first;
	if(cumProb < cutOff) { vscores[origIndex] = (int)vscores.size(); }
	else { vscores[origIndex] = (int)vscores.size() - i; }
	// for development purposes ...
	if(cumProb < cutOff) { topProbIDs.push_back( pair<int,double>(origIndex, currProb) ); }

	cumProb += currProb;
    }


    COLOR_DEBUG_PRINT("SORT STUFF");
    for(int i = 0; i < (int)vscoreIndex.size(); i++) { COLOR_DEBUG_PRINT(i<<"\t"<<vscoreIndex[i].first<<"\t"<<vscoreIndex[i].second); }

    COLOR_DEBUG_PRINT("PROB SCORE:\t"<<vscores);
    return vscores;
}
/**
 * use RF for mutation scoring
 */
//vector<double> ColorProfile::getRFMutationScoring() {
//
//	probScores.resize(colorProfile.size(), 0);
//	vector<double> vscores(colorProfile.size(), 0);
//	if(colorProfile.size() == 0) return vscores;
//
//	// for each profile...
//	for(int i = 0; i < (int)colorProfile.size(); i++) {
//		// get start pos of read on reference
//		int posPad = refPosStarts[i];
//		// compute by summing matches
//		for(int j = 10; j < (int)colorProfile[i].size() && ((posPad+j) < 300); j++) {
//			string lmer(refSeqs[i].substr(j-3,8));
//
//			double prob = 0.0;
//			// if find mismatch, use mismatch l-mer
//			if(std::find(lmer.begin(), lmer.end(), '.') != lmer.end()) {
//				prob = 0.035;
//			}
//			else {		// otherwise, use l-mer of reference
//				double p = (MutationModel::getInstance()).getProb(lmer, posPad+j);
//				if(colorProfile[i][j] == '.') { prob = p; }
//				else { prob = (1.0-p); }
//			}
//			vscores[i] += log2(prob);
//
//			COLORPROFILE_DEBUG_PRINT(i<<"\t"<<j<<"\t"<<posPad+j<<"\t"<<lmer<<"\t"<<colorProfile[i][j]<<"\t"<<prob<<"\t"<<log2(prob)<<"\t"<<vscores[i]);
//		}
//	}
//	COLORPROFILE_DEBUG_PRINT("VSCORES:\t"<<vscores);
//	double probMut = -numeric_limits<double>::infinity();
//	COLORPROFILE_DEBUG_PRINT("PROB MUT:\t"<<probMut);
//	for(int i = 0; i < (int)vscores.size(); i++) { if(vscores[i] == 0) { vscores[i] = probMut; } }
//
//	// scale/transform back to probabilities
//	double maxProb = *std::max_element(vscores.begin(), vscores.end());							// find max
//	for(int i = 0; i < (int)vscores.size(); i++) { vscores[i] -= maxProb; }						// subtract max from all
//	double total = 0.0;
//	for(int i = 0; i < (int)vscores.size(); i++) { total += pow(2.0, vscores[i]); }				// compute denom
//	for(int i = 0; i < (int)vscores.size(); i++) { vscores[i] = pow(2.0,vscores[i]) / total; }	// compute prob using difference
//	for(int i = 0; i < (int)vscores.size(); i++) { probScores[i] = vscores[i]; }
//
//	// now rank-score (in reverse, i.e., highest prob -> highest rank) such that can use standard score picking scheme
//	vector<pair<int,double>> vscoreIndex(vscores.size());
//	for(int i = 0; i < (int)vscores.size(); i++) { vscoreIndex[i] = pair<int,double>(i, vscores[i]); }
//	std::sort(vscoreIndex.begin(), vscoreIndex.end(), [](const pair<int,double> a, const pair<int,double> b){ return a.second > b.second; });
//	// the top percent will tie for highest rank, all others rank in order, ignoring ties
//	double cutOff = ColorProfile::THRESH_P;   //0.8; //0.999; //0.80;	// 0.95
//	double cumProb = 0.0;
//	for(int i = 0; i < (int)vscoreIndex.size(); i++) {
//		double currProb = vscoreIndex[i].second;
//		int origIndex = vscoreIndex[i].first;
//		if(cumProb < cutOff) { vscores[origIndex] = (int)vscores.size(); }
//		else { vscores[origIndex] = (int)vscores.size() - i; }
//		// for development purposes ...
//		if(cumProb < cutOff) { topProbIDs.push_back( pair<int,double>(origIndex, currProb) ); }
//
//		cumProb += currProb;
//	}
//
//
//	COLOR_DEBUG_PRINT("SORT STUFF");
//	for(int i = 0; i < (int)vscoreIndex.size(); i++) { COLOR_DEBUG_PRINT(i<<"\t"<<vscoreIndex[i].first<<"\t"<<vscoreIndex[i].second); }
//
//	COLOR_DEBUG_PRINT("PROB SCORE:\t"<<vscores);
//	return vscores;
//}
/**
 * simplest scoring scheme, +1 for a match, 0 otherwise
 * @returns vector of integer scores for each profile in class
 */
//vector<int> ColorProfile::getSimpleScoring() {
vector<double> ColorProfile::getSimpleScoring() {
//	vector<int> vscores(colorProfile.size(), 0);
    vector<double> vscores(colorProfile.size(), 0);
    // for each profile...
    for(int i = 0; i < (int)colorProfile.size(); i++) {
	// compute by summing matches
	for(int j = 0; j < (int)colorProfile[i].size(); j++) if(colorProfile[i][j] != '.') vscores[i] += 1;
    }
    return vscores;
}
/**
 * return partition of each reference as <start,end> pair
 * @returns a vector of pairs, each pair representing the start_pos and end_pos of the reference assignment to the read
 */
vector<pair<int,int>> ColorProfile::getPartitions() {
    vector<pair<int,int>> partitionPairs(colorProfile.size(), pair<int,int>(-1,-1));
    for(int i = 0; i < (int)colorProfile.size(); i++) {
	int len = (int)colorProfile[i].size();
	// find index of first match
	int strt_pos = 0;
	for( ; colorProfile[i][strt_pos] == '.' && strt_pos < len; strt_pos++) { }
	// find index of last match
	int end_pos = len-1;
	for( ; colorProfile[i][end_pos] == '.' && end_pos > 0; end_pos--) { }
	partitionPairs[i].first = strt_pos;
	partitionPairs[i].second = end_pos;
    }

    return partitionPairs;
}
//vector<int> ColorProfile::getSimpleMotifScoring() {
vector<double> ColorProfile::getSimpleMotifScoring() {
//	vector<int> vscores(colorProfile.size(), 0);
    vector<double> vscores(colorProfile.size(), 0.0);
    // for each profile...
    for(int i = 0; i < (int)colorProfile.size(); i++) {

	COLOR_DEBUG_PRINT(i<<":\t"<<colorProfile[i]<<"\n\t"<<refSeqs[i]);
	// compute by summing matches
	for(int j = 0; j < (int)colorProfile[i].size(); j++) {
	    //if(colorProfile[i][j] != '.') vscores[i] += 1;
	    if(colorProfile[i][j] != '.') vscores[i] += 1;
	    // mismatch, check if higher chance of mutation
	    else if(colorProfile[i][j] == '.' && j >= 2) {
		string lmer2Back(refSeqs[i].substr(j-2, 4));
		string lmer1Back(refSeqs[i].substr(j-1, 4));
		if(matchesWRCYMotif(lmer2Back)) vscores[i] += 3;
		else if(matchesRGYWMotif(lmer1Back)) vscores[i] += 3;
		else if(matchesWANMotif(lmer1Back)) vscores[i] += 2;
		else {}

		COLOR_DEBUG_PRINT(i<<"\t"<<j<<"\t"<<lmer2Back<<"\t"<<matchesWRCYMotif(lmer2Back));
	    }
	    else {}
	}
    }
    return vscores;
}

/**
 * checks if matches WRCY motif [A/T][A/G]C[C/T]
 */
bool ColorProfile::matchesWRCYMotif(string lmer) {
    if(lmer.size() != 4) return false;
    if( !(lmer[0] == 'A' || lmer[0] == 'T') ) return false;
    if( !(lmer[1] == 'A' || lmer[1] == 'G') ) return false;
    if(lmer[2] != 'C') return false;
    if( !(lmer[3] == 'C' || lmer[3] == 'T') ) return false;
    return true;
}
/**
 * checks if matches RGYW motif [A/G]G[C/T][A/T]
 */
bool ColorProfile::matchesRGYWMotif(string lmer) {
    if(lmer.size() != 4) return false;
    if( !(lmer[0] == 'A' || lmer[0] == 'G') ) return false;
    if( lmer[1] != 'G' ) return false;
    if( !(lmer[2] == 'C' || lmer[2] == 'T') ) return false;
    if( !(lmer[3] == 'A' || lmer[3] == 'T') ) return false;
    return true;
}
/**
 * checks if matches WAN motif
 */
bool ColorProfile::matchesWANMotif(string lmer) {
    if(lmer.size() != 4) return false;
    if( !(lmer[0] == 'A' || lmer[0] == 'T') ) return false;
    if(lmer[1] != 'A') return false;
    return true;
}
