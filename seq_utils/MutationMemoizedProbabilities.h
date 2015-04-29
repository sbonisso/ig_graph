/*
 * MutationMemoizedProbabilities.h
 *
 *  Created on: Jan 27, 2014
 *      Author: sbonisso
 */

#ifndef MUTATIONMEMOIZEDPROBABILITIES_H_
#define MUTATIONMEMOIZEDPROBABILITIES_H_

#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>

class MutationMemoizedProbabilities {
private:
	const static int MAX_LEN = 500;
	const static int MAX_COLOR = 500;

	const static int LMER_LEN =   4; //8; //6; //4;
	const static int LMER_SHIFT = 1; //3; //2; //1;

	double pMutMatrix[MAX_COLOR][MAX_LEN];
	double pNullVect[MAX_LEN];
	int currColor;

//	static MutationMemoizedProbabilities *mp;

	MutationMemoizedProbabilities() {
		currColor = 0;
		addNullProbs();
	}
//	MutationMemoizedProbabilities(int nColors) {
//		numColors = nColors;
//	}
	MutationMemoizedProbabilities(MutationMemoizedProbabilities const&);
	void operator=(MutationMemoizedProbabilities const&);
	/**
	 * get prob of mutation given position, irrespective of l-mer
	 */
	void addNullProbs() {
		string nullLmerStr(LMER_LEN, '.');
		for(int i = 0; i < MAX_LEN; i++) { pNullVect[i] = 0.0; }
		for(int i = 0; i < 300; i++) {
			pNullVect[i] = (MutationNBProbabilities::getInstance()).getProb(nullLmerStr, i);
		}
	}

public:

	static MutationMemoizedProbabilities& getInstance() {
		static MutationMemoizedProbabilities *mP = new MutationMemoizedProbabilities();
		return *mP;
//		if(!mp) {
//			mp = new MutationMemoizedProbabilities();
//		}
//		return *mp;
	}
	/**
	 * add a color (reference sequence) to the memoized probability matrix
	 */
	void addColor(std::string seq) {
		// zero out row
		for(int i = 0; i < MAX_LEN; i++) { pMutMatrix[currColor][i] = 0.0; }

		for(int pos = LMER_SHIFT; pos < (int)seq.size()-(LMER_LEN-LMER_SHIFT); pos++) {
			std::string lmer = seq.substr(pos-LMER_SHIFT, LMER_LEN);
			// use other class to compute probabilities from model
			pMutMatrix[currColor][pos] = ((MutationNBProbabilities::getInstance()).getProb(lmer, pos));
//			pMutMatrix[currColor][pos] = (MutationModel::getInstance()).getProb(lmer, pos);

//			pMutMatrix[currColor][pos] = (MutationProbabilities::getInstance()).getProb(lmer, pos);
//			cout<<currColor<<"\t"<<lmer<<"\t"<<pos<<"\t"<<pMutMatrix[currColor][pos]<<endl;
		}
		currColor++;
	}

	int getNumColors() { return currColor; }
	/**
	 * get probability of mutation given l-mer, given pos, of a specified reference sequence index
	 */
	double getProb(int colorIndex, int pos) { return pMutMatrix[colorIndex][pos]; }
	/**
	 * get probability of mutation given pos, irrespective of l-mer
	 */
	double getProbNull(int pos) { return pNullVect[pos]; }
};


#endif /* MUTATIONMEMOIZEDPROBABILITIES_H_ */
