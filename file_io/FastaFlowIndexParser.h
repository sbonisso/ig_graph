/*
 * FastaFlowIndexParser.h
 *
 *  Created on: Jun 25, 2013
 *      Author: sbonisso
 */

#ifndef FASTAFLOWINDEXPARSER_H_
#define FASTAFLOWINDEXPARSER_H_

#include <sstream>
#include <utility>
#include <string>
#include <boost/regex.hpp>

using namespace std;

template <class FP>
class FastaFlowIndexParser : public FP {
	vector<int> currFlowgram;
	const int FLOW_VECT_SIZE = 500;
	void zeroOut();
	int currFlowLen;

public:
	FastaFlowIndexParser() : FP() { }
	FastaFlowIndexParser(string fastaFilePath);
	virtual ~FastaFlowIndexParser();

	int getCurrFlowLength() { return currFlowLen; }
	virtual pair<string,string> getNextEntry();
	vector<int> *getCurrentFlowgram();
};


template <class FP> FastaFlowIndexParser<FP>::FastaFlowIndexParser(string fastaFilePath) : FP(fastaFilePath) {
	currFlowgram.resize(FLOW_VECT_SIZE,0);
}
template <class FP> FastaFlowIndexParser<FP>::~FastaFlowIndexParser() {

}

template <class FP> vector<int> *FastaFlowIndexParser<FP>::getCurrentFlowgram() {
	return &currFlowgram;
}

template <class FP> void FastaFlowIndexParser<FP>::zeroOut() {
	for(int i = 0; i < FLOW_VECT_SIZE; i++) { currFlowgram[i] = 0; }
}

template <class FP> pair<string,string> FastaFlowIndexParser<FP>::getNextEntry() {
	// get the entry
	pair<string,string> entry = FP::getNextEntry();
	string readID = entry.first;
	string aryStr = entry.second;

	zeroOut();
	currFlowLen = 0;
	boost::regex re("(\\d|\\.|\\d)+");
	boost::sregex_token_iterator iter(aryStr.begin(), aryStr.end(), re, 0);
	boost::sregex_token_iterator endIter;
	for(int i =0  ; iter != endIter; ++iter, i++, currFlowLen++) {
		// convert to integer
		int flowIndex;
		stringstream(*iter)>>flowIndex;
		currFlowgram[i] = flowIndex;
	}

	return entry;
}



#endif /* FASTAFLOWINDEXPARSER_H_ */
