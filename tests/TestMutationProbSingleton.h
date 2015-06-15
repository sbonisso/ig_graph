/*
 * TestMutationProbSingleton.h
 *
 *  Created on: Nov 4, 2013
 *      Author: sbonisso
 */

#ifndef TESTMUTATIONPROBSINGLETON_H_
#define TESTMUTATIONPROBSINGLETON_H_


#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>

//#include "../graphs/ColoredAntibodyGraphMultiple.h"
//#include "ParameterizedGraph.h"
//#include "../ParameterizedColorGraph.h"

#include "file_io/FastaParser.h"
#include "graphs/ColorProfile.h"
#include "seq_utils/MutationNBProbabilities.h"
#include "seq_utils/Encoding.hpp"

using namespace std;

class TestMutationProbSingleton : public Test::Suite
{
public:
    TestMutationProbSingleton()
    {
//		TEST_ADD(TestColoredAntibodyGraph::test_label_1)
//		TEST_ADD(TestMutationProbSingleton::test_nb_probs)
	TEST_ADD(TestMutationProbSingleton::test_lmer_indexer);
	//TEST_ADD(TestMutationProbSingleton::test_nb_probs_against_pMat);
	TEST_ADD(TestMutationProbSingleton::test_lmer_encode);
    }
    
private:

    void test_nb_probs_against_pMat();
    void test_lmer_indexer();
    void test_nb_probs();
    void test_lmer_encode();
};

#endif /* TESTMUTATIONPROBSINGLETON_H_ */
