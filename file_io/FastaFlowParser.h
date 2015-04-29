/*
 * FastaFlowParser.h
 *
 *  Created on: Jun 25, 2013
 *      Author: sbonisso
 */

#ifndef FASTAFLOWPARSER_H_
#define FASTAFLOWPARSER_H_

#include <sstream>
#include <utility>
#include <string>
#include <boost/regex.hpp>

using namespace std;

template <class FP>
class FastaFlowParser : public FP {
	vector<double> currFlowgram;
	const int FLOW_VECT_SIZE = 500;
	void zeroOut();
	int currFlowLen;

public:
	FastaFlowParser() : FP() { currFlowLen = -1; }
	FastaFlowParser(string fastaFilePath);
	virtual ~FastaFlowParser();

	int getCurrFlowLength() { return currFlowLen; }
	virtual pair<string,string> getNextEntry();
	vector<double> *getCurrentFlowgram();
};

template <class FP> FastaFlowParser<FP>::FastaFlowParser(string fastaFilePath) : FP(fastaFilePath) {
	currFlowgram.resize(FLOW_VECT_SIZE,0.0);
}
template <class FP> FastaFlowParser<FP>::~FastaFlowParser() {

}

template <class FP> vector<double> *FastaFlowParser<FP>::getCurrentFlowgram() {
	return &currFlowgram;
}

template <class FP> void FastaFlowParser<FP>::zeroOut() {
	for(int i = 0; i < FLOW_VECT_SIZE; i++) { currFlowgram[i] = 0.0; }
}

template <class FP> pair<string,string> FastaFlowParser<FP>::getNextEntry() {
	// get the entry
	pair<string,string> entry = FP::getNextEntry();
	string readID = entry.first;
	string aryStr = entry.second;
	aryStr = aryStr.substr(1,aryStr.size()-2);

	zeroOut();
	currFlowLen = 0;
	boost::regex re("(\\d|\\.|\\d)+");
	boost::sregex_token_iterator iter(aryStr.begin(), aryStr.end(), re, 0);
	boost::sregex_token_iterator endIter;

	//for(int i =0  ; iter != endIter; ++iter, i++, currFlowLen++) {
	for(int i =0  ; iter != endIter; iter++, i++, currFlowLen++) {
		// convert to integer
		double flowVal;
		stringstream(*iter)>>flowVal;
		currFlowgram[i] = flowVal;
//		cout<<i<<"\t"<<flowVal<<"\t"<<currFlowgram[i]<<endl;
	}
	if(currFlowLen < (int)currFlowgram.size()) { currFlowgram.resize(currFlowLen, 0.0); }

	return entry;
}

#endif /* FASTAFLOWPARSER_H_ */
