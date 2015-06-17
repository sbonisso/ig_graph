#include "TestMutationNBModel.hpp"

void TestMutationNBModel::test_model_similarity_along_pos() {
    string homeDir(getenv("HOME"));
    (MutationNBProbabilities::getInstance()).setParamDir(homeDir+"/.nb_params/4mer_amp/");
    //(MutationNBModel::getInstance()).setParamDir("data/model.json");
    (MutationNBModel::getInstance()).setParamDir("data/model.bin");
    
    // double val_1 = (MutationNBProbabilities::getInstance()).getProb("TCCG", 215);
    // double val_2 = (MutationNBModel::getInstance()).getProb("TCCG", 215);
    
    for(int i = 30; i < 250; i++) {
	double v1 = (MutationNBProbabilities::getInstance()).getProb("TCCG", i);
	double v2 = (MutationNBModel::getInstance()).getProb("TCCG", i);
	TEST_ASSERT(abs(v1-v2) < 0.00001);
    }
}


void TestMutationNBModel::test_model_similarity_across_lmers() {
    string homeDir(getenv("HOME"));
    (MutationNBProbabilities::getInstance()).setParamDir(homeDir+"/.nb_params/4mer_amp/");
    //(MutationNBModel::getInstance()).setParamDir("data/model.json");
    (MutationNBModel::getInstance()).setParamDir("data/model.bin");
    
    int pos = 91;
    vector<string> lmers { "AAAA", "ACTG", "ATGT", "ACCC",
	    "CATG", "CCCC", "CGTG", 
	    "GATG", "GCTG", "GTGG", "GGGG",
	    "TATA", "TACG", "TCAT", "TTTT"}; 
    for(int i = 0; i < (int)lmers.size(); i++) {
	double v1 = (MutationNBProbabilities::getInstance()).getProb(lmers[i], pos);
	double v2 = (MutationNBModel::getInstance()).getProb(lmers[i], pos);
	TEST_ASSERT(abs(v1 - v2) < 0.00001);
    }    
}
