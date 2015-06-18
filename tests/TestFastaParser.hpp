#ifndef TESTFASTAPARSER_H_
#define TESTFASTAPARSER_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>

#include "prettyprint.hpp"

#include "file_io/FastaParser.hpp"
#include "file_io/FastaRefID.hpp"

using namespace std;

class TestFastaParser : public Test::Suite {
public:
TestFastaParser();
virtual ~TestFastaParser();

protected:
virtual void setup();      		// setup resources...
    
void test_parse_smab_file();
void test_parse_smab_file_ext();
};

#endif
