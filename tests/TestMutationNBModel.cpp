#include "catch.hpp"
#include <stdlib.h>

#include <utility>
#include <string>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/external/rapidjson/filestream.h>

#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>

#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/memory.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/array.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/bitset.hpp>
#include <cereal/types/polymorphic.hpp>

#include "file_io/FastaParser.hpp"

#include "seq_utils/MutationNBModel.hpp"
#include "seq_utils/Encoding.hpp"

TEST_CASE("MutationNBModel similarity pos", "[mut_model,pos]") {
    (MutationNBModel::getInstance()).setParamDir("data/model.bin");
    
    vector<double> pos_vect(250,0);
    std::ifstream is("tests/baselines/mut_pos.json");
    cereal::JSONInputArchive ar(is);
    ar(pos_vect);
     
    for(int i = 30; i < 250; i++) {
	double v2 = (MutationNBModel::getInstance()).getProb("TCCG", i);
	double v1 = pos_vect[i];
	REQUIRE(abs(v1-v2) < 0.00001);
    }
}

TEST_CASE("MutationNBModel similarity lmer", "[mut_model,lmer]") {
    (MutationNBModel::getInstance()).setParamDir("data/model.bin");

    map<string,double> lmer_prob;
    std::ifstream is("tests/baselines/mut_lmer_pos91.json");
    cereal::JSONInputArchive ar(is);
    ar(lmer_prob);
    
    int pos = 91;
    vector<string> lmers { "AAAA", "ACTG", "ATGT", "ACCC",
	    "CATG", "CCCC", "CGTG", 
	    "GATG", "GCTG", "GTGG", "GGGG",
	    "TATA", "TACG", "TCAT", "TTTT"}; 
    for(int i = 0; i < (int)lmers.size(); i++) {
	double v2 = (MutationNBModel::getInstance()).getProb(lmers[i], pos);	
	double v1 = lmer_prob[lmers[i]];
	REQUIRE(abs(v1 - v2) < 0.00001);
	lmer_prob[lmers[i]] = v1;
    }
}
