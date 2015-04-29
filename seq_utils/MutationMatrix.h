/*
 * MutationMatrix.h
 *
 *  Created on: Nov 15, 2013
 *      Author: sbonisso
 */

#ifndef MUTATIONMATRIX_H_
#define MUTATIONMATRIX_H_

#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>

class MutationMatrix {
private:
	static const int NUM_LMERS = 256;

	MutationMatrix() {}
	MutationMatrix(std::string filePath) {
		std::string line;
		std::ifstream in(filePath.c_str());
		while(getline(in, line)) {
			std::istringstream iss(line);
//			std::string lmer;
			//int pos;
			std::string lmer1;
			std::string lmer2;
			double val;
//			iss >> lmer >> pos >> val;
			iss >> lmer1 >> lmer2 >> val;
			mutMatrix[getLMerIndex(lmer1)][getLMerIndex(lmer2)] = val;
		}
	}
	MutationMatrix(MutationProbabilities const&);
	void operator=(MutationProbabilities const&);


	// map of [pos,l-mer] -> score
	double mutMatrix[NUM_LMERS][NUM_LMERS];


public:

	static MutationMatrix& getInstance() {
		static MutationMatrix *mP = new MutationMatrix("/tmp/loMat_4mers_v2.tab");
		return *mP;
	}

//	double getProb(std::string lmer, int pos) { return pMutMatrix[getLMerIndex(lmer)][pos]; }
	double getScore(std::string lmer1, std::string lmer2) { return mutMatrix[getLMerIndex(lmer1)][getLMerIndex(lmer2)]; }

	int powInt(int base, int power) { return power == 0 ? 1 : powInt(base, power-1)*base; }

	int getLMerIndex(std::string lmer) {
		int index = 0;
		int powI = (int)lmer.size()-1;
		for(int i = 0; i < (int)lmer.size(); i++) {
			if(lmer[i] == 'A') index += 0*powInt(4, powI);
			else if(lmer[i] == 'C') index += 1*powInt(4, powI);
			else if(lmer[i] == 'G') index += 2*powInt(4, powI);
			else if(lmer[i] == 'T') index += 3*powInt(4, powI);
			else {}
			powI--;
		}
		return index;
	}
};



#endif /* MUTATIONMATRIX_H_ */
