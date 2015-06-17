/*
 * TestMutationProbSingleton.cpp
 *
 *  Created on: Nov 4, 2013
 *      Author: sbonisso
 */

#include "TestMutationProbSingleton.hpp"

void TestMutationProbSingleton::test_lmer_indexer() {

	// TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("TTTT") == 255)

	// TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AAAA") == 0)

	// TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AAAC") == 1)

	// TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AAAG") == 2)

	// TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AAAT") == 3)

	// TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AACA") == 4)

	// TEST_ASSERT((MutationProbabilities::getInstance()).getLMerIndex("AACT") == 7)
}

void TestMutationProbSingleton::test_nb_probs() {
	string homeDir(getenv("HOME"));
	(MutationNBProbabilities::getInstance()).setParamDir(homeDir+"/.nb_params/6mer_amp/");

	// for 6-mers
	double val = (MutationNBProbabilities::getInstance()).getProb("AAAGAA", 91);
	//cout<<"P(mut|AAACTT,91) = "<<val<<endl;
	TEST_ASSERT(abs(val - 0.939394) < 0.000001);

	double posVal = (MutationNBProbabilities::getInstance()).getProb("......", 91);
	//cout<<"P(mut |......,91) = "<<posVal<<endl;
	TEST_ASSERT(abs(posVal - 0.709924) < 0.000001);
	
	//double val2 = (MutationNBProbabilities::getInstance()).getProb("TGAACT", 46);
	//cout<<"P(mut|TGAACT,46) = "<<val2<<endl;

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
	
	// // for 6-mers
	// double val1 = (MutationNBProbabilities::getInstance()).getProb("AGAA", 91);
	// double val2 = (MutationProbabilities::getInstance()).getProb("AGAA", 91);
	// TEST_ASSERT(abs(val1 - val2) < 0.000001)

	// double val1 = (MutationNBProbabilities::getInstance()).getProb("....", 91);
	// double val2 = (MutationProbabilities::getInstance()).getProb("....", 91);
	// TEST_ASSERT(abs(val1 - val2) < 0.000001);
}

void TestMutationProbSingleton::test_lmer_encode() {
    string homeDir(getenv("HOME"));
    (MutationNBProbabilities::getInstance()).setParamDir(homeDir+"/.nb_params/4mer_amp/");
	
    double v3 = (MutationNBProbabilities::getInstance()).getProb("AAAA", 91);
    
    //unsigned int agaa_i = Encoding::encode_kmer("AGAA");
    //double v1 = (MutationNBProbabilities::getInstance()).getProb(agaa_i, 91);
    unsigned int aaaa_i = Encoding::encode_kmer("AAAA");
    double v2 = (MutationNBProbabilities::getInstance()).getProb(aaaa_i, 91);
    TEST_ASSERT(abs(v2 - v3) < 0.00001);
    
    //double val1 = (MutationNBProbabilities::getInstance()).getProb("AGAA", 91);
    //double val2 = (MutationProbabilities::getInstance()).getProb("AGAA", 91);
    // cout<<"\nAGAA\t"<<agaa_i<<"\t"<<v1<<"\t"<<val1<<endl;
    // cout<<"AAAA\t"<<aaaa_i<<"\t"<<v2<<"\t"<<v3<<endl;
    // cout<<"max val:\t"<<(MutationNBProbabilities::getInstance()).getMaxVal()<<endl;

    //CATG	258	0.23085	0.463547
    unsigned int catg_i = Encoding::encode_kmer("CATG");
    double catg_v1 = (MutationNBProbabilities::getInstance()).getProb(catg_i, 258);
    double catg_v2 = (MutationNBProbabilities::getInstance()).getProb(catg_i, 258);
    //cout<<"\nCATG\t"<<catg_i<<"\t"<<catg_v1<<"\t"<<catg_v2<<endl;
    TEST_ASSERT(abs(catg_v1 - catg_v2) < 0.00001);
}
