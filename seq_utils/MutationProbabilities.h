/*
 * MutationProbabilities.h
 *
 *  Created on: Nov 4, 2013
 *      Author: sbonisso
 */

#ifndef MUTATIONPROBABILITIES_H_
#define MUTATIONPROBABILITIES_H_

#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>

class MutationProbabilities {
private:
	static const int NUM_LMERS = 256;
	static const int NUM_POS = 400;

	MutationProbabilities() {}
	MutationProbabilities(std::string filePath) {
		std::string line;
		std::ifstream in(filePath.c_str());
		while(getline(in, line)) {
			std::istringstream iss(line);
			std::string lmer;
			int pos;
			double prob;
			iss >> lmer >> pos >> prob;
//			pMutMap[lmer][pos] = prob;
//			cout<<lmer<<"\t"<<pos<<"\t"<<prob<<endl;
			pMutMatrix[getLMerIndex(lmer)][pos] = prob;
		}
	}
	MutationProbabilities(MutationProbabilities const&);
	void operator=(MutationProbabilities const&);
//	~MutationProbabilities() {}

	// map of [pos,l-mer] -> prob mutation
	//std::unordered_map<std::string, std::unordered_map<int, double> > pMutMap;
//	std::map<std::string,int> lmerToIndex;
	double pMutMatrix[NUM_LMERS][NUM_POS];


public:

	static MutationProbabilities& getInstance() {
		//static MutationProbabilities *mP = new MutationProbabilities("/tmp/pMat_mut.tab");
		static MutationProbabilities *mP = new MutationProbabilities("/tmp/pMat_mut_pos4.tab");
//		static MutationProbabilities *mP = new MutationProbabilities("pmat.tab"); // was this one "/tmp/pMat_mut_pos4.tab"
//		static MutationProbabilities *mP = new MutationProbabilities("/tmp/pMat_mut_pos4_noAmp.tab");
//		static MutationProbabilities *mP = new MutationProbabilities("/tmp/pMat_mut_pos4_dependent.tab");
		return *mP;
	}

//	double getProb(std::string lmer, int pos) { return pMutMap[lmer][pos]; }
	double getProb(std::string lmer, int pos) { return pMutMatrix[getLMerIndex(lmer)][pos]; }

	double getScore(std::string lmer, int pos) { return this->getProb(lmer, pos); }

	int powInt(int base, int power) { return power == 0 ? 1 : powInt(base, power-1)*base; }

	int getLMerIndex(string lmer) {
		int index = 0;
		int powI = (int)lmer.size()-1;
		for(int i = 0; i < (int)lmer.size(); i++) {
			if(lmer[i] == 'A') index += 0*powInt(4, powI);
			else if(lmer[i] == 'C') index += 1*powInt(4, powI);
			else if(lmer[i] == 'G') index += 2*powInt(4, powI);
			else if(lmer[i] == 'T') index += 3*powInt(4, powI);
			else {}
//			cout<<i<<"\t"<<powI<<"\t"<<lmer[i]<<"\t"<<index<<endl;
			powI--;
		}
		return index;
	}
};


#endif /* MUTATIONPROBABILITIES_H_ */
