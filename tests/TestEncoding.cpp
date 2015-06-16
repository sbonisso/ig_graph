#include "TestEncoding.hpp"

/**
 * test encoding k-mers
 */
void TestEncoding::test_lmer_encode() {
    unsigned int aaaa_i = Encoding::encode_kmer("AAAA");
    TEST_ASSERT(aaaa_i == 0);
    //
    unsigned int catg_i = Encoding::encode_kmer("CATG");
    TEST_ASSERT(catg_i == 177);
    
}
/**
 * test decoding k-mers
 */
void TestEncoding::test_lmer_decode() {
    unsigned int aaaa_i = Encoding::encode_kmer("AAAA");
    string aaaa_s = Encoding::decode_kmer(aaaa_i, 4);
    TEST_ASSERT(aaaa_s == "AAAA");
    //
    unsigned int catg_i = Encoding::encode_kmer("CATG");
    string catg_s = Encoding::decode_kmer(catg_i, 4);
    TEST_ASSERT(catg_s == "CATG");
}
/**
 * test updating encoded k-mers
 */
void TestEncoding::test_lmer_update() {
    string seq = "ACGACTACTTTGACTGCAT";
    int len = (int)seq.size();
    int l = 4;
    
    unsigned int lmer_i = Encoding::encode_kmer(seq.substr(0,l));
    
    for(int i = 1; i < len-l; i++) {
	string lmer( seq.substr(i,l) );
	lmer_i = Encoding::update_kmer(lmer_i, lmer[l-1], l);
	
	string curr_kmer = Encoding::decode_kmer(lmer_i, l);
	TEST_ASSERT(curr_kmer == lmer);
    }
}
