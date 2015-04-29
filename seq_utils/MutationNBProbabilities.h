/*
 * MutationNBProbabilities.h
 *
 *  Created on: Jan 24, 2014
 *      Author: sbonisso
 */

#ifndef MUTATIONNBPROBABILITIES_H_
#define MUTATIONNBPROBABILITIES_H_

#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "Utils.h"

#ifdef DEBUG
#define MUTNBPROB_DEBUG_PRINT(EXPR) DEBUG_PRINT("[MUT_NB_PROB]", EXPR)
#else
#define MUTNBPROB_DEBUG_PRINT(x) do {} while (0)
#endif

class MutationNBProbabilities {
private:
	MutationNBProbabilities();
	~MutationNBProbabilities();
	
	// hashes to hold class conditional probabilities
	std::map<std::string,std::map<std::string, int> > lmerCounts;
	std::map<std::string,std::map<int, int> > posCounts;

	std::map<std::string,std::map<std::string, double> > lmerProbs;
	std::map<std::string,std::map<int, double> > posProbs;
	std::map<std::string,double> priorProbs;

	std::string paramDir;
	int lmer_len;

	// to read in data from files
	void readInData();
	void fillLMerHash(std::string &filePath);
	void fillPosHash(std::string &filePath);
	void fillPriors(std::string &priorPath);

	// to control singleton
	static bool _instanceExists;
	static MutationNBProbabilities* _nbInstance;

public:

	static MutationNBProbabilities& getInstance();
	void setParamDir(std::string paramDirStr);
	int getLMerLen();

	double getProb(std::string lmer, int pos);
	double getScore(std::string lmer, int pos);
};

#endif /* MUTATIONNBPROBABILITIES_H_ */
