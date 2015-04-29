/*
 * MutationScoring.h
 *
 *  Created on: Jan 24, 2014
 *      Author: sbonisso
 */

#ifndef MUTATIONSCORING_H_
#define MUTATIONSCORING_H_

#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>

class MutationScoring {
private:

	//virtual ~MutationScoring() {}

public:

	// will be singleton class, so to enforce that use this getInstance
	virtual static MutationScoring& getInstance() = 0;

	// getScore will yield the score/probability of mutation of an l-mer and position
	virtual double getScore(std::string lmer, int pos) = 0;

	// utility method for computing a power of an integer
	int powInt(int base, int power) { return power == 0 ? 1 : powInt(base, power-1)*base; }

	// hashing for 4-mers
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



#endif /* MUTATIONSCORING_H_ */
