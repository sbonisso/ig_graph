/*
 * FastaErrors.h
 *
 *  Created on: Mar 28, 2013
 *      Author: stef
 */

#ifndef FASTAERRORS_H_
#define FASTAERRORS_H_

#include <sstream>
#include <utility>
#include <string>
#include <vector>
#include <boost/regex.hpp>

using namespace std;

template <class FP>
class FastaErrors : public FP {
	string currRefID;

	vector<int> entryInsertions;
	vector<int> entryDeletions;
	int startPos;
	int endPos;

public:
	FastaErrors() : FP() {}
	FastaErrors(string fastaFilePath);
	virtual ~FastaErrors() {}

	virtual pair<string,string> getNextEntry();
//	pair<string,string> getNextEntry();
	vector<int> getEntryInsertions() { return this->entryInsertions; }
	vector<int> getEntryDeletions() { return this->entryDeletions; }
	int getStartPos() { return this->startPos; }
	int getEndPos() { return this->endPos; }
};


template <class FP> FastaErrors<FP>::FastaErrors(string fastaFilePath) : FP(fastaFilePath) {
	currRefID = "";
}

template <class FP> pair<string,string> FastaErrors<FP>::getNextEntry() {
	// get the entry
	pair<string,string> entry = FP::getNextEntry();
	string readIDNum = entry.first;
	string seq = entry.second;

	this->entryInsertions.clear();
	this->entryDeletions.clear();
	this->startPos = -1;
	this->endPos = -1;

	string readID = FP::getCurrIDLine().substr(1);

	// extract start/end positions of amplicon
	boost::regex startEndPosRE("position=complement\\((\\d+)\\.\\.(\\d+)\\)");
	boost::smatch seMatch;
	string::const_iterator seStart = readID.begin();
	string::const_iterator seEnd   = readID.end();
	boost::regex_search(seStart, seEnd, seMatch, startEndPosRE);
	string startPosStr(seMatch[1].first, seMatch[1].second);
	string endPosStr(seMatch[2].first, seMatch[2].second);
	this->startPos = atoi(startPosStr.c_str());
	this->endPos = atoi(endPosStr.c_str());

	// extract list of pos/error type
	boost::regex re("errors=((\\d+[-|+][A|C|G|T]?,?)+)");
	boost::smatch match;
	string::const_iterator start = readID.begin();
	string::const_iterator end   = readID.end();
	boost::regex_search(start, end, match, re);
	string stest(match[1].first, match[1].second);
//	cout<<"STEST:\t"<<stest<<endl;

	// parse the positions for deletion errors
	boost::regex posRE("(\\d+)-");
	string::const_iterator dstart = stest.begin();
	string::const_iterator dend   = stest.end();
	boost::smatch numMatch;
	boost::regex_search(dstart, dend, numMatch, posRE);

	string nums(numMatch[1].first, numMatch[1].second);
	boost::sregex_token_iterator iter(dstart, dend, posRE, 0);
	boost::sregex_token_iterator send;
	for( ; iter != send; ++iter ) {
		string tmp(*iter);
//		cout<<"tmp:\t"<<tmp<<endl;
		this->entryDeletions.push_back(atoi(tmp.c_str()));
	}
//	cout<<"DELETIONS:\t"<<this->entryDeletions<<endl;

	// now parse the positions for insertion errors
	boost::regex posREIns("(\\d+)\\+");
	string::const_iterator istart = stest.begin();
	string::const_iterator iend   = stest.end();
	boost::smatch numInsMatch;
	boost::regex_search(istart, iend, numInsMatch, posREIns);

	boost::sregex_token_iterator ins_iter(istart, iend, posREIns, 0);
	boost::sregex_token_iterator ins_end;
	for( ; ins_iter != ins_end; ++ins_iter ) {
		string tmp(*ins_iter);
		entryInsertions.push_back(atoi(tmp.c_str()));
	}
//	cout<<"INSERTIONS:\t"<<this->entryInsertions<<endl;

	return entry;
}

#endif /* FASTAERRORS_H_ */
