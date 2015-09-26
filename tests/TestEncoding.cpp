#include "catch.hpp"
#include <stdlib.h>

#include <utility>
#include <string>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>

#include "seq_utils/Encoding.hpp"

using namespace std;

/**
 * test encoding k-mers
 */
TEST_CASE("test Encoding encode kmers", "[encode]") {
    unsigned int aaaa_i = Encoding::encode_kmer("AAAA");
    REQUIRE(aaaa_i == 0);
    //
    unsigned int catg_i = Encoding::encode_kmer("CATG");
    REQUIRE(catg_i == 177);
}
/**
 * test decoding k-mers
 */
TEST_CASE("test Encoding decode kmers", "[decode]") {
    unsigned int aaaa_i = Encoding::encode_kmer("AAAA");
    string aaaa_s = Encoding::decode_kmer(aaaa_i, 4);
    REQUIRE(aaaa_s == "AAAA");
    //
    unsigned int catg_i = Encoding::encode_kmer("CATG");
    string catg_s = Encoding::decode_kmer(catg_i, 4);
    REQUIRE(catg_s == "CATG");
}
/**
 * test updating encoded k-mers
 */
TEST_CASE("test Encoding updating kmers", "[update]") {
    string seq = "ACGACTACTTTGACTGCAT";
    int len = (int)seq.size();
    int l = 4;
    
    unsigned int lmer_i = Encoding::encode_kmer(seq.substr(0,l));
    
    for(int i = 1; i < len-l; i++) {
	string lmer( seq.substr(i,l) );
	lmer_i = Encoding::update_kmer(lmer_i, lmer[l-1], l);
	
	string curr_kmer = Encoding::decode_kmer(lmer_i, l);
	REQUIRE(curr_kmer == lmer);
    }
}
