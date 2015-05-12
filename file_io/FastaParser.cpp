/*
 * FastaParser.cpp
 *
 *  Created on: Jan 10, 2013
 *      Author: stef
 */

#include "FastaParser.h"
#include <cstring>

/*
 * @param fastaFilePath A string denoting the path of the FASTA file
 */
FastaParser::FastaParser(string fastaFilePath) {
	this->filename = fastaFilePath;
}

/*
 * @return A boolean if the file was opened properly or not
 */
bool FastaParser::openFile() {
    this->fin.open(this->filename.c_str());
    //return true;
    return this->fin.is_open();
}

/*
 * @return A string of the next sequence in the FASTA file
 */
string FastaParser::getNextSequence() {
	string str;
	string tmpStr;
	do {
		getline(this->fin, tmpStr);
		if(tmpStr[0] != '>') str += tmpStr;
	} while(tmpStr[0] != '>');

	return str;
}

/*
 * @return A pair of strings representing the ID and sequence of the next entry of the FASTA file
 */
pair<string,string> FastaParser::getNextEntry() {
	pair<string,string> currEntry;
	currEntry.first = ""; currEntry.second = "";
	string str = ""; string tmpStr = "";
	do {
//		if(this->fin.eof()) { currEntry.second = str+tmpStr; return currEntry; }
		if(this->fin.eof()) { currEntry.second = str; return currEntry; }	// placed before getline for last line - so that don't duplicate
		getline(this->fin, tmpStr);		// get next line

		if(tmpStr[0] == '>') {		// if an entry id, parse it
			currIDLine = tmpStr; // save current ID line
			int next = tmpStr.find_first_of(" ");
			currEntry.first = tmpStr.substr(1,next);
		}
		else { str += tmpStr; }		// if a sequence, cat it
	} while( fin.peek() != '>');	// peek at char, if new entry, end (will always get 1st, end on 2nd)

	currEntry.second = str;
	return currEntry;
}

bool FastaParser::hasNextSequence() { return !this->fin.eof(); }

bool FastaParser::closeFile() {
	return true;
}


