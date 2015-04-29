/*
 * FastaRefID.h
 *
 *  Created on: Jan 30, 2013
 *      Author: stef
 */

#ifndef FASTAREFID_H_
#define FASTAREFID_H_

#include <sstream>
#include <utility>
#include <string>
#include <boost/regex.hpp>

using namespace std;

template <class FP>
class FastaRefID : public FP {
	string currRefID;

public:
	FastaRefID() : FP() {}
	FastaRefID(string fastaFilePath);
	virtual ~FastaRefID() {}

	string getCurrentRefID() { return currRefID; }
	virtual pair<string,string> getNextEntry();
//	pair<string,string> getNextEntry();
};


template <class FP> FastaRefID<FP>::FastaRefID(string fastaFilePath) : FP(fastaFilePath) {
	currRefID = "";
}

template <class FP> pair<string,string> FastaRefID<FP>::getNextEntry() {
	// get the entry
	pair<string,string> entry = FP::getNextEntry();
	string readID = entry.first;
	string seq = entry.second;

	// boost regex matching of reference ID
//	boost::regex re("reference=([\\w|\\d|\\.|\\:|\\,]+)");
	boost::regex re("reference=([\\w|\\d|\\.|\\:|\\,|\\*|\\-]+)");
	boost::smatch match;
	string::const_iterator start = this->currIDLine.begin();
	string::const_iterator end   = this->currIDLine.end();
	boost::regex_search(start, end, match, re);
	string stest(match[1].first, match[1].second);
	currRefID = stest;
//	cout<<stest<<endl;
//	cout<<this->currIDLine<<endl;
//	cout<<"pos: "<<pos<<"\tspacepos: "<<spacePos<<endl;

	return entry;
}



#endif /* FASTAREFID_H_ */
