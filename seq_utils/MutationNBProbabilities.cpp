/*
 * MutationNBProbabilities.cpp
 *
 *  Created on: Jan 31, 2014
 *      Author: sbonisso
 */

#include "MutationNBProbabilities.h"

bool MutationNBProbabilities::_instanceExists;
MutationNBProbabilities* MutationNBProbabilities::_nbInstance;

MutationNBProbabilities& MutationNBProbabilities::getInstance() {
    if(!_instanceExists) {
	_nbInstance = new MutationNBProbabilities();
	_instanceExists = true;
    }
    return (*_nbInstance);
}

MutationNBProbabilities::MutationNBProbabilities() { lmer_len = 0; }

MutationNBProbabilities::~MutationNBProbabilities() { this->_instanceExists = false; }
/**
 *
 */
void MutationNBProbabilities::setParamDir(std::string paramDirStr) {
    this->paramDir = paramDirStr;
    this->readInData();
    MUTNBPROB_DEBUG_PRINT("SET PARAM DIR:\t"<<paramDir);
    MUTNBPROB_DEBUG_PRINT("AGAA:\t"<<this->lmerProbs["match"]["AGAA"]
			  <<"\t"<<this->lmerProbs["mutation"]["AGAA"]);
    MUTNBPROB_DEBUG_PRINT("91:\t"<<this->posProbs["match"][91]
			  <<"\t"<< this->posProbs["mutation"][91]);
    MUTNBPROB_DEBUG_PRINT("priors\t"<<this->priorProbs["match"]
			  <<"\t"<<priorProbs["mutation"]);
}
/**
 *
 */
void MutationNBProbabilities::readInData() {	
    std::string lmerPath = paramDir + "lmerCounts.tab";
    std::string posPath = paramDir + "posCounts.tab";
    std::string priorPath = paramDir + "priorProbs.tab";

    this->fillLMerHash(lmerPath);
    this->fillPosHash(posPath);
    this->fillPriors(priorPath);
}
/**
 *
 */
double MutationNBProbabilities::getProb(std::string lmer, int pos) {
    double mutProb = posProbs["mutation"][pos]*priorProbs["mutation"];
    double matchProb = posProbs["match"][pos]*priorProbs["match"];
    
    // is '....' l-mer, signalling a position only prob
    if(lmer[0]=='.') {	}
    else {
	double lmerMut =lmerProbs["mutation"][lmer];
	double lmerMatch = lmerProbs["match"][lmer];
	
	// if l-mer is not present, leave it out
	if(!(lmerMut == 0 && lmerMatch == 0)) {
	    mutProb *= lmerMut;
	    matchProb *= lmerMatch;
	}
    }

    double currProb = mutProb;
    double denomProb = mutProb + matchProb;
    
    return (currProb / denomProb);
}

double MutationNBProbabilities::getProb(unsigned int lmer_i, int pos) {
    double mutProb = pos_mut_[pos]*prior_mut_;
    double matchProb = pos_nomut_[pos]*prior_nomut_;
    
    // is '....' l-mer, signalling a position only prob
    //if(lmer_i < 0) {	}
    if(lmer_i == 999999) {}
    else {
	//double lmerMut =lmerProbs["mutation"][lmer];
	//double lmerMatch = lmerProbs["match"][lmer];
	double lmerMut = lmer_mut_[lmer_i];
	double lmerMatch = lmer_nomut_[lmer_i];
	
	// if l-mer is not present, leave it out
	if(!(lmerMut == 0 && lmerMatch == 0)) {
	    mutProb *= lmerMut;
	    matchProb *= lmerMatch;
	}
    }

    double currProb = mutProb;
    double denomProb = mutProb + matchProb;

    return (currProb / denomProb);
}


double MutationNBProbabilities::getScore(std::string lmer, int pos) { return this->getProb(lmer, pos); }

int MutationNBProbabilities::getLMerLen() { return this->lmer_len; }
/**
 *
 */
void MutationNBProbabilities::fillPriors(std::string &priorPath) {
    std::string line;
    std::ifstream in(priorPath.c_str());
    while(getline(in, line)) {
	std::istringstream iss(line);
	std::string lmer;
	double matchP; double mutP;
	iss >> lmer >> matchP >> mutP;
	priorProbs["match"] = matchP;
	priorProbs["mutation"] = mutP;
	prior_nomut_ = matchP;
	prior_mut_ = mutP;
    }
}
/**
 *
 */
void MutationNBProbabilities::fillLMerHash(std::string &filePath) {
    std::string line;
    std::ifstream in(filePath.c_str());
    int totalMut = 0;
    int totalMatch=0;
    while(getline(in, line)) {
	std::istringstream iss(line);
	if(line[0] == '#') continue;	// skip header - should be commented out
	std::string lmer;
	int matchCount; int mutCount;
	iss >> lmer >> matchCount >> mutCount;
	// to avoid any -inf from 0's
	if(matchCount == 0 || mutCount == 0) { matchCount++; mutCount++; }	
	lmerCounts["match"][lmer] = matchCount;
	lmerCounts["mutation"][lmer] = mutCount;
	lmer_len = (int)lmer.size();
	totalMatch += matchCount;
	totalMut += mutCount;
    }
    
    string max_lmer;
    max_lmer.assign(lmer_len, 'T');
    max_val_ = Encoding::encode_kmer(max_lmer);
    
    lmer_nomut_.resize(lmerCounts["match"].size());
    lmer_mut_.resize(lmerCounts["mutation"].size());
    for(std::unordered_map<std::string,int>::iterator itr = lmerCounts["match"].begin(); 
	itr != lmerCounts["match"].end();
	itr++) {
	std::string lmer = itr->first;
	int numMatch = itr->second;
	lmerProbs["match"][lmer] = ((double)numMatch/((double)totalMatch));
	unsigned int lmer_i = Encoding::encode_kmer(lmer);
	lmer_nomut_[lmer_i] = ((double)numMatch/((double)totalMatch));
    }
    for(std::unordered_map<std::string,int>::iterator itr = lmerCounts["mutation"].begin(); 
	itr != lmerCounts["mutation"].end(); 
	itr++) {
	std::string lmer = itr->first;
	int numMut = itr->second;
	lmerProbs["mutation"][lmer] = ((double)numMut/((double)totalMut));
	unsigned int lmer_i = Encoding::encode_kmer(lmer);
	lmer_mut_[lmer_i] = ((double)numMut/((double)totalMut));
    }

}
/**
 *
 */
void MutationNBProbabilities::fillPosHash(std::string &filePath) {
    std::string line;
    std::ifstream in(filePath.c_str());
    int totalMut = 0;
    int totalMatch=0;
//	std::cout<<"FILL POS HASH"<<std::endl;
    while(getline(in, line)) {
	std::istringstream iss(line);
	if(line[0] == '#') continue;	// skip header - should be commented out
	int pos;
	int matchCount; int mutCount;
	iss >> pos >> matchCount >> mutCount;
	posCounts["match"][pos] = matchCount;
	posCounts["mutation"][pos] = mutCount;
	totalMatch += matchCount;
	totalMut += mutCount;
    }
    pos_nomut_.resize(400,0);
    pos_mut_.resize(400,0);
    for(std::unordered_map<int,int>::iterator itr = posCounts["match"].begin(); 
	itr != posCounts["match"].end();
	itr++) {
	int pos = itr->first;
	int numMatch = itr->second;
	posProbs["match"][pos] = ((double)numMatch/((double)totalMatch));
	pos_nomut_[pos] = ((double)numMatch/((double)totalMatch));
    }
    for(std::unordered_map<int,int>::iterator itr = posCounts["mutation"].begin(); 
	itr != posCounts["mutation"].end(); 
	itr++) {
	int pos = itr->first;
	int numMut = itr->second;
	posProbs["mutation"][pos] = ((double)numMut/((double)totalMut));
	pos_mut_[pos] = ((double)numMut/((double)totalMut));
    }
}
