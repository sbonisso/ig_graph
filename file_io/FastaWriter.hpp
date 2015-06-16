/*
 * FastaWriter.h
 *
 *  Created on: Jan 10, 2013
 *      Author: stef
 */

#ifndef FASTAWRITER_H_
#define FASTAWRITER_H_

#include <string>
#include <iostream>
#include <fstream>

using namespace std;

class FastaWriter {
public:
	FastaWriter(string);
	~FastaWriter();  // close file on destructor

	void closeFasta();
	int addEntry(string, string);
	string getPath() { return this->pathToFile; }
private:
	int numEntries;
	string pathToFile;
	ofstream fOut;
};


#endif /* FASTAWRITER_H_ */
