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
#include <unordered_map>
#include <vector>

#include "seq_utils/Encoding.hpp"

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
	std::unordered_map<std::string,std::unordered_map<std::string, int> > lmerCounts;
	std::unordered_map<std::string,std::unordered_map<int, int> > posCounts;

	std::unordered_map<std::string,std::unordered_map<std::string, double> > lmerProbs;
	std::unordered_map<std::string,std::unordered_map<int, double> > posProbs;
	std::unordered_map<std::string,double> priorProbs;

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
	
	std::vector<double> pos_mut_;
	std::vector<double> pos_nomut_;
	std::vector<double> lmer_mut_;
	std::vector<double> lmer_nomut_;
	double prior_mut_;
	double prior_nomut_;
	
public:

	static MutationNBProbabilities& getInstance();
	void setParamDir(std::string paramDirStr);
	int getLMerLen();

	double getProb(std::string lmer, int pos);
	double getScore(std::string lmer, int pos);

	double getProb(unsigned int lmer_i, int pos);
	int getMaxVal() { return max_val_; }
	
	int max_val_;

};

#endif /* MUTATIONNBPROBABILITIES_H_ */
