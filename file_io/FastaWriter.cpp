/*
 * FastaWriter.cpp
 *
 *  Created on: Jan 10, 2013
 *      Author: stef
 */

#include "FastaWriter.hpp"

FastaWriter::FastaWriter(string filePath) {
	this->numEntries = 0;
	this->pathToFile = filePath;
	this->fOut.open(this->pathToFile.c_str(), ios::out);
}

FastaWriter::~FastaWriter() { this->closeFasta(); }

/*
 * returns 0 if went well, -1 if problem writing
 */
int FastaWriter::addEntry(string id, string seq) {
	if(!this->fOut) { return -1; }
	this->fOut<<">"<<id<<endl;
	this->fOut<<seq<<endl;
	return 0;
}

/*
 * closes file opened during constructor
 */
void FastaWriter::closeFasta() { this->fOut.close(); }
