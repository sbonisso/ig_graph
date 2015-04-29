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

//	std::cout<<mutProb<<" = "<<posProbs["mutation"][pos]<<" * "<<priorProbs["mutation"]<<std::endl;
//	std::cout<<matchProb<<" = "<<posProbs["match"][pos]<<" * "<<priorProbs["match"]<<std::endl;

    // is '....' l-mer, signalling a position only prob
    if(lmer[0]=='.') {	}
    else {
	double lmerMut =lmerProbs["mutation"][lmer];
	double lmerMatch = lmerProbs["match"][lmer];


	// if l-mer is not present, leave it out
	if(!(lmerMut == 0 && lmerMatch == 0)) {
	    mutProb *= lmerMut;
	    matchProb *= lmerMatch;

//			std::cout<<mutProb<<" = "<<lmerMut<<" * "<<posProbs["mutation"][pos]<<" * "<<priorProbs["mutation"]<<std::endl;
//			std::cout<<matchProb<<" = "<<lmerMatch<<" * "<<posProbs["match"][pos]<<" * "<<priorProbs["match"]<<std::endl;
	}
    }

    double currProb = mutProb;
    double denomProb = mutProb + matchProb;

//	std::cout<<(currProb / denomProb)<<" = "<<currProb<<" / "<<denomProb<<std::endl;

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
//			double matchProb = (double)((double)matchCount / (double)(matchCount+mutCount));
//			double mutProb = (double)((double)mutCount / (double)(matchCount+mutCount));
//			cout<<lmer<<"\t"<<matchCount<<"\t"<<matchProb<<"\t"<<mutCount<<"\t"<<mutProb<<endl;
	if(matchCount == 0 || mutCount == 0) { matchCount++; mutCount++; }	// to avoid any -inf from 0's
	lmerCounts["match"][lmer] = matchCount;
	lmerCounts["mutation"][lmer] = mutCount;
	lmer_len = (int)lmer.size();
	totalMatch += matchCount;
	totalMut += mutCount;
//			cout<<"FILL LMER\t"<<lmer<<"\t"<<matchCount<<"\t"<<mutCount<<"\t"<<totalMatch<<endl;
//			lmerProbs["match"][lmer] = matchProb;
//			lmerProbs["mutation"][lmer] = mutProb;
    }
    for(std::map<std::string,int>::iterator itr = lmerCounts["match"].begin(); 
	itr != lmerCounts["match"].end();
	itr++) {
	std::string lmer = itr->first;
	int numMatch = itr->second;
	lmerProbs["match"][lmer] = ((double)numMatch/((double)totalMatch));
    }
    for(std::map<std::string,int>::iterator itr = lmerCounts["mutation"].begin(); 
	itr != lmerCounts["mutation"].end(); 
	itr++) {
	std::string lmer = itr->first;
	int numMut = itr->second;
	lmerProbs["mutation"][lmer] = ((double)numMut/((double)totalMut));
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
//			double matchProb = (double)((double)matchCount / (double)(matchCount+mutCount));
//			double mutProb = (double)((double)mutCount / (double)(matchCount+mutCount));
//			posProbs["match"][pos] = matchProb;
//			posProbs["mutation"][pos] = mutProb;
	posCounts["match"][pos] = matchCount;
	posCounts["mutation"][pos] = mutCount;
	totalMatch += matchCount;
	totalMut += mutCount;
//		std::cout<<totalMatch<<"\t"<<totalMut<<"\t"<<line<<std::endl;
//			cout<<pos<<"\t"<<matchCount<<"\t"<<matchProb<<"\t"<<mutCount<<"\t"<<mutProb<<"\t"<<totalMatch<<endl;
//			cout<<pos<<"\t"<<matchCount<<"\t"<<mutCount<<"\t"<<totalMatch<<endl;
    }
    for(std::map<int,int>::iterator itr = posCounts["match"].begin(); 
	itr != posCounts["match"].end();
	itr++) {
	int pos = itr->first;
	int numMatch = itr->second;
	posProbs["match"][pos] = ((double)numMatch/((double)totalMatch));
//		if(pos == 91) std::cout<<"match"<<"\t"<<pos<<"\t"<<numMatch<<" / "<<totalMatch<<" = "<<posProbs["match"][pos]<<std::endl;
    }
    for(std::map<int,int>::iterator itr = posCounts["mutation"].begin(); 
	itr != posCounts["mutation"].end(); 
	itr++) {
	int pos = itr->first;
	int numMut = itr->second;
	posProbs["mutation"][pos] = ((double)numMut/((double)totalMut));
//		if(pos == 91) std::cout<<"mutation"<<"\t"<<pos<<"\t"<<numMut<<" / "<<totalMatch<<" = "<<posProbs["mutation"][pos]<<std::endl;
    }
}
