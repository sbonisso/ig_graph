#include "catch.hpp"
#include <stdlib.h>

#include <utility>
#include <string>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>

#include "file_io/FastaParser.hpp"

#include "seq_utils/MutationNBProbabilities.h"
#include "seq_utils/MutationNBModel.hpp"
#include "seq_utils/Encoding.hpp"

TEST_CASE("MutationNBModel similarity pos", "[mut_model,pos]") {
    string homeDir(getenv("HOME"));
    (MutationNBProbabilities::getInstance()).setParamDir(homeDir+"/.nb_params/4mer_amp/");
    (MutationNBModel::getInstance()).setParamDir("data/model.bin");
    
    // double val_1 = (MutationNBProbabilities::getInstance()).getProb("TCCG", 215);
    // double val_2 = (MutationNBModel::getInstance()).getProb("TCCG", 215);
    
    for(int i = 30; i < 250; i++) {
	double v1 = (MutationNBProbabilities::getInstance()).getProb("TCCG", i);
	double v2 = (MutationNBModel::getInstance()).getProb("TCCG", i);
	REQUIRE(abs(v1-v2) < 0.00001);
    }
}

TEST_CASE("MutationNBModel similarity lmer", "[mut_model,lmer]") {
    string homeDir(getenv("HOME"));
    (MutationNBProbabilities::getInstance()).setParamDir(homeDir+"/.nb_params/4mer_amp/");
    (MutationNBModel::getInstance()).setParamDir("data/model.bin");
    
    int pos = 91;
    vector<string> lmers { "AAAA", "ACTG", "ATGT", "ACCC",
	    "CATG", "CCCC", "CGTG", 
	    "GATG", "GCTG", "GTGG", "GGGG",
	    "TATA", "TACG", "TCAT", "TTTT"}; 
    for(int i = 0; i < (int)lmers.size(); i++) {
	double v1 = (MutationNBProbabilities::getInstance()).getProb(lmers[i], pos);
	double v2 = (MutationNBModel::getInstance()).getProb(lmers[i], pos);
	REQUIRE(abs(v1 - v2) < 0.00001);
    }
}
