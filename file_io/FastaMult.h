/*
 * FastaMult.h
 *	Templated class to additionally parse multiplicity information
 *	contained in the ID line of a FASTA file
 *
 *  Created on: Jan 30, 2013
 *      Author: stef
 */

#ifndef FASTAMULT_H_
#define FASTAMULT_H_

#include <sstream>
#include <utility>
#include <string>

using namespace std;

template <class FP>
class FastaMult : public FP {
	int currMult;

public:
	FastaMult() : FP() {}
	FastaMult(string fastaFilePath);
	virtual ~FastaMult() {}

	int getCurrentMultiplicity() { return currMult; }
	virtual pair<string,string> getNextEntry();
};


template <class FP> FastaMult<FP>::FastaMult(string fastaFilePath) : FP(fastaFilePath) {
	currMult = -1;
}

template <class FP> pair<string,string> FastaMult<FP>::getNextEntry() {
	// get the entry
	pair<string,string> entry = FP::getNextEntry();
	string readID = entry.first;
	string seq = entry.second;
	// parse for multiplicity in ID string
	//int pos = readID.find_first_of("num=",0);
	int pos = this->currIDLine.find_first_of("num=",0);
	string numStr = this->currIDLine.substr(pos+4,readID.size()-pos);
	// convert to integer
	int multNum;
	stringstream(numStr)>>multNum;
	// save in object
	currMult = multNum;

	return entry;
}

#endif /* FASTAMULT_H_ */
