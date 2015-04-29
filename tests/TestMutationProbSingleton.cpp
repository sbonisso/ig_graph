/*
 * TestMutationProbSingleton.cpp
 *
 *  Created on: Nov 4, 2013
 *      Author: sbonisso
 */



#include "TestMutationProbSingleton.h"

void TestMutationProbSingleton::test_lmer_indexer() {

//	cout<<(MutationProbabilities::getInstance()).getLMerIndex("TTTT")<<endl;
	TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("TTTT") == 255)

//	cout<<(MutationProbabilities::getInstance()).getLMerIndex("AAAA")<<endl;
	TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AAAA") == 0)

//	cout<<(MutationProbabilities::getInstance()).getLMerIndex("CTTT")<<endl;
//	cout<<(MutationProbabilities::getInstance()).getLMerIndex("AAAC")<<endl;
	TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AAAC") == 1)

//	cout<<(MutationProbabilities::getInstance()).getLMerIndex("AAAG")<<endl;
	TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AAAG") == 2)

//	cout<<(MutationProbabilities::getInstance()).getLMerIndex("AAAT")<<endl;
	TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AAAT") == 3)

//	cout<<(MutationProbabilities::getInstance()).getLMerIndex("AACA")<<endl;
	TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AACA") == 4)

//	cout<<(MutationProbabilities::getInstance()).getLMerIndex("AACT")<<endl;

//	cout<<(MutationProbabilities::getInstance()).getLMerIndex("ACCT")<<endl;
	TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AACT") == 7)


	cout<<"P(AACT|91) = "<<(MutationProbabilities::getInstance()).getProb("AACT", 91)<<endl;
}

void TestMutationProbSingleton::test_nb_probs() {
	string homeDir(getenv("HOME"));
	(MutationNBProbabilities::getInstance()).setParamDir(homeDir+"/.nb_params/6mer_amp/");

	// for 6-mers
	double val = (MutationNBProbabilities::getInstance()).getProb("AAAGAA", 91);
	cout<<"P(mut|AAACTT,91) = "<<val<<endl;
	TEST_ASSERT(abs(val - 0.939394) < 0.000001)
//	double val = (MutationNBProbabilities::getInstance()).getProb("AACT", 91, "mutation");
//	cout<<"SECOND TEST:\t"<<val<<endl;

	double posVal = (MutationNBProbabilities::getInstance()).getProb("......", 91);
	cout<<"P(mut |......,91) = "<<posVal<<endl;
	TEST_ASSERT(abs(posVal - 0.709924) < 0.000001)

	double val2 = (MutationNBProbabilities::getInstance()).getProb("TGAACT", 46);
	cout<<"P(mut|TGAACT,46) = "<<val2<<endl;

	double val3 = (MutationNBProbabilities::getInstance()).getProb("GCGAGA", 42);
	cout<<"P(mut|GCGAGA,43) = "<<val3<<endl;

	// this l-mer is not present, so should use position only, i.e., '......'
	double val4 = (MutationNBProbabilities::getInstance()).getProb("GCTAGG", 43);
	cout<<"P(mut|GCTAGG,43) = "<<val4<<endl;
//	TEST_ASSERT(abs(val4 - 0.0213071) < 0.000001)

	double val4Pos = (MutationNBProbabilities::getInstance()).getProb("......", 43);
	cout<<"P(mut |......,43) = "<<val4Pos<<endl;
//	TEST_ASSERT(abs(val4Pos - 0.0213071) < 0.000001)

	double val5 = (MutationNBProbabilities::getInstance()).getProb("TCCGTA", 215);
	cout<<"P(mut |TCCGTA,215) = "<<val5<<endl;
	//TEST_ASSERT(abs(val5 - 0.0213071) < 0.000001)

}

void TestMutationProbSingleton::test_nb_probs_against_pMat() {
	string homeDir(getenv("HOME"));
	(MutationNBProbabilities::getInstance()).setParamDir(homeDir+"/.nb_params/4mer_amp/");
	
	// for 6-mers
	double val1 = (MutationNBProbabilities::getInstance()).getProb("AGAA", 91);
	double val2 = (MutationProbabilities::getInstance()).getProb("AGAA", 91);
	cout<<"P(mut|AGAA,91) = "<<val1<<"\t"<<val2<<endl;
	TEST_ASSERT(abs(val1 - val2) < 0.000001)

	val1 = (MutationNBProbabilities::getInstance()).getProb("....", 91);
	val2 = (MutationProbabilities::getInstance()).getProb("....", 91);
	cout<<"P(mut|....,91) = "<<val1<<"\t"<<val2<<endl;
	TEST_ASSERT(abs(val1 - val2) < 0.000001)

//	double val2 = (MutationNBProbabilities::getInstance()).getProb("TGAACT", 46);
//	cout<<"P(mut|TGAACT,46) = "<<val2<<endl;
//
//	double val3 = (MutationNBProbabilities::getInstance()).getProb("GCGAGA", 42);
//	cout<<"P(mut|GCGAGA,43) = "<<val3<<endl;
//
//	// this l-mer is not present, so should use position only, i.e., '......'
//	double val4 = (MutationNBProbabilities::getInstance()).getProb("GCTAGG", 43);
//	cout<<"P(mut|GCTAGG,43) = "<<val4<<endl;
////	TEST_ASSERT(abs(val4 - 0.0213071) < 0.000001)
//
//	double val4Pos = (MutationNBProbabilities::getInstance()).getProb("......", 43);
//	cout<<"P(mut |......,43) = "<<val4Pos<<endl;
////	TEST_ASSERT(abs(val4Pos - 0.0213071) < 0.000001)
//
//	double val5 = (MutationNBProbabilities::getInstance()).getProb("TCCGTA", 215);
//	cout<<"P(mut |TCCGTA,215) = "<<val5<<endl;
//	//TEST_ASSERT(abs(val5 - 0.0213071) < 0.000001)

}