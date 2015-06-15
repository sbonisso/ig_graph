/*
 * TestEncoding.h
 *
 *  Created on: 
 *      Author: sbonisso
 */

#ifndef TESTENCODING_H_
#define TESTENCODING_H_


#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>

#include "seq_utils/Encoding.hpp"

using namespace std;

class TestEncoding : public Test::Suite
{
public:
    TestEncoding()
    {
	TEST_ADD(TestEncoding::test_lmer_encode);
	TEST_ADD(TestEncoding::test_lmer_decode);
	TEST_ADD(TestEncoding::test_lmer_update);
    }
    
private:

    void test_lmer_encode();
    void test_lmer_decode();
    void test_lmer_update();
};

#endif /* TESTMUTATIONPROBSINGLETON_H_ */
