#include "TestCreateProfile.h"
#include "../prettyprint.hpp"

TestCreateProfile::TestCreateProfile() { 
    TEST_ADD(TestCreateProfile::test_cp_preds_1);
    TEST_ADD(TestCreateProfile::test_cp_scores_1);
    
    TEST_ADD(TestCreateProfile::test_cp_preds_2);
    TEST_ADD(TestCreateProfile::test_cp_scores_2);

    TEST_ADD(TestCreateProfile::test_max_n);
    
    TEST_ADD(TestCreateProfile::test_cdr3_1);
    TEST_ADD(TestCreateProfile::test_cdr3_2);
}

TestCreateProfile::~TestCreateProfile() {}
/**
 * setup tests with human refs
 */
void TestCreateProfile::setup() {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    string d_fasta = "data/igh_refs_simple/human_IGHD.fa";
    string j_fasta = "data/igh_refs_simple/human_IGHJ.fa";
    //
    cab = new CanonicalAntibodyGraph(21);
    cab->addVReferences(v_fasta);
    cab->addDReferences(d_fasta);
    cab->addJReferences(j_fasta);   
}

void TestCreateProfile::tear_down() { delete cab; }
/**
 *
 */
void TestCreateProfile::test_cp_preds_1() {
    //ifstream in("tests/data/trunc_seq.txt");
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    CreateProfile cp(cab, 2);
    //
    //ColorProfileMatrix cp_mat = cp.getColorProfile(seq);
    cp.compute(seq);
    //
    vector<string> vpred = cp.getPredictedV();
    vector<string> dpred = cp.getPredictedD();
    vector<string> jpred = cp.getPredictedJ();
    //
    vector<string> vtruth = {"IGHV1-18*01", "IGHV1-18*03"};
    vector<string> dtruth = {"?", "?"}; //{"IGHD1-1*01", "IGHD1-20*01"};
    vector<string> jtruth = {"IGHJ1*01", "IGHJ2*01"};
    //
    for(int i = 0; i < 2; i++) { TEST_ASSERT(vpred[0] == vtruth[0]); }
    for(int i = 0; i < 2; i++) { TEST_ASSERT(dpred[0] == dtruth[0]); }
    for(int i = 0; i < 2; i++) { TEST_ASSERT(jpred[0] == jtruth[0]); }
}
/**
 *
 */
void TestCreateProfile::test_cp_scores_1() {
    //ifstream in("tests/data/trunc_seq.txt");
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);    
    //
    CreateProfile cp(cab);
    //
    //ColorProfileMatrix cp_mat = cp.getColorProfile(seq);
    cp.compute(seq);
    //
    vector<double> vscore = cp.getVScores();
    vector<double> dscore = cp.getDScores();
    vector<double> jscore = cp.getJScores();
    //   
    vector<int> v_truth = {296, 295, 295, 265, 269, 267, 266, 267, 246, 268, 266, 250, 251, 264, 263, 263, 254, 255, 243, 241, 244, 242, 244, 243, 244, 245, 243, 243, 244, 243, 271, 272, 240, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 198, 199, 0, 202, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0, 0, 0, 204, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 143, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 141, 139, 140, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 200, 28, 26, 26, 0, 0, 0, 0, 0, 0, 0, 0, 170, 166, 171, 0, 169, 169, 170, 28, 24, 0, 28, 28, 28, 28, 28, 29, 29, 29, 24, 0, 0, 0, 0, 0, 0, 0, 133, 135, 0, 0, 167, 168, 0, 0, 0, 29, 29, 0, 164, 137, 137, 139, 138, 135, 0, 0, 132, 137, 137, 138, 29, 29, 29, 27, 29, 167, 168, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 231, 234, 233, 232};//with tip prop
    //vector<int> v_truth = {296, 295, 295, 206, 252, 234, 249, 208, 226, 210, 225, 106, 107, 247, 246, 241, 213, 214, 243, 241, 244, 242, 244, 243, 244, 245, 243, 243, 244, 243, 254, 255, 218, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 34, 32, 0, 34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 0, 0, 0, 0, 29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30, 28, 24, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 29, 23, 21, 21, 0, 0, 0, 0, 0, 0, 0, 0, 32, 32, 34, 0, 32, 32, 32, 23, 24, 0, 23, 23, 23, 23, 23, 29, 29, 29, 24, 0, 0, 0, 0, 0, 0, 0, 29, 27, 0, 0, 32, 32, 0, 0, 0, 29, 29, 0, 34, 29, 29, 29, 29, 24, 0, 0, 24, 27, 27, 27, 29, 29, 29, 27, 29, 32, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 120, 213, 212, 120};
    vector<double> d_truth = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<double> j_truth = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    //
    for(int i = 0; i < (int)vscore.size(); i++) {
	TEST_ASSERT_MSG(vscore[i] == v_truth[i], "V scores are incorrect");
    }
    for(int i = 0; i < (int)dscore.size(); i++) {
	TEST_ASSERT_MSG(dscore[i] == d_truth[i], "D scores are incorrect");
    }
    for(int i = 0; i < (int)jscore.size(); i++) {
	TEST_ASSERT_MSG(jscore[i] == j_truth[i], "J scores are incorrect");
    }
}
/**
 *
 */
void TestCreateProfile::test_cp_preds_2() {
    ifstream in("tests/data/ten_ab.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string read_id;  string seq = ""; 
    for(int i = 0; i < 2; i++) {
	getline(in, read_id);
	getline(in, seq);
    }
    //
    CreateProfile cp(cab, 2);
    //
    //ColorProfileMatrix cp_mat = cp.getColorProfile(seq);
    cp.compute(seq);
    //
    vector<string> vpred = cp.getPredictedV();
    vector<string> dpred = cp.getPredictedD();
    vector<string> jpred = cp.getPredictedJ();
    //
    vector<string> v_truth = {"IGHV1-45*02", "IGHV1-45*01"};
    vector<string> d_truth = {"IGHD2-2*01", "IGHD2-2*02"};
    vector<string> j_truth = {"IGHJ3*02", "IGHJ3*01"};
    //
    for(int i = 0; i < 2; i++) { TEST_ASSERT(vpred[0] == v_truth[0]); }
    for(int i = 0; i < 2; i++) { TEST_ASSERT(dpred[0] == d_truth[0]); }
    for(int i = 0; i < 2; i++) { TEST_ASSERT(jpred[0] == j_truth[0]); }
}
/**
 *
 */
void TestCreateProfile::test_cp_scores_2() {
    ifstream in("tests/data/ten_ab.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string read_id;  string seq = ""; 
    for(int i = 0; i < 2; i++) {
	getline(in, read_id);
	getline(in, seq);
    }
    //
    CreateProfile cp(cab);
    //
    //ColorProfileMatrix cp_mat = cp.getColorProfile(seq);
    cp.compute(seq);
    //
    vector<double> vscore = cp.getVScores();
    vector<double> dscore = cp.getDScores();
    vector<double> jscore = cp.getJScores();
    //
    vector<double> v_truth = {249, 248, 250, 256, 258, 256, 257, 256, 244, 252, 251, 292, 293, 257, 256, 257, 243, 243, 248, 231, 248, 248, 249, 231, 249, 249, 247, 247, 248, 248, 252, 253, 241, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 220, 221, 220, 219};
    // vector<double> d_truth = {0, 0, 0, 0, 0, 27, 27, 26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};    
    // vector<double> j_truth = {0, 0, 49, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0};    
    //vector<double> v_truth = {107, 107, 107, 223, 224, 222, 223, 223, 217, 224, 233, 292, 293, 239, 238, 239, 219, 219, 230, 213, 229, 229, 231, 213, 231, 230, 228, 228, 229, 230, 234, 235, 212, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 96, 96, 96, 96};
    vector<double> d_truth = {0, 0, 0, 0, 0, 25, 25, 23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<double> j_truth = {0, 0, 49, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0};    
    //
    for(int i = 0; i < (int)vscore.size(); i++) {
    	TEST_ASSERT_MSG(vscore[i] == v_truth[i], "V scores are incorrect");
    }
    for(int i = 0; i < (int)dscore.size(); i++) {
    	TEST_ASSERT_MSG(dscore[i] == d_truth[i], "D scores are incorrect");
    }
    for(int i = 0; i < (int)jscore.size(); i++) {
    	TEST_ASSERT_MSG(jscore[i] == j_truth[i], "J scores are incorrect");
    }
}
/**
 * test finding maximal n indeces
 */
void TestCreateProfile::test_max_n() {
    vector<double> test_scores = {39, 0, 0, 0, 43, 42, 43, 48, 46, 0, 0, 0, 0};
    
    CreateProfile cp(cab);
    vector<int> max_index = cp.getNMaxIndex(3, test_scores);
    
    TEST_ASSERT(max_index[0] == 7);
    TEST_ASSERT(max_index[1] == 8);
    TEST_ASSERT(max_index[2] == 4);    
}
/**
 * test for making sure will properly display no found CDR3
 * for a truncated sequence
 */
void TestCreateProfile::test_cdr3_1() {
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);            
    //
    CreateProfile cp(cab, 2);
    //    
    cp.compute(seq);
    string cdr3_seq = cp.getCDR3();
    TEST_ASSERT_MSG(cdr3_seq == "?", "CDR3 seq incorrect");    
}
/**
 * test for capturing the CDR3 of a sequence
 */
void TestCreateProfile::test_cdr3_2() {
    ifstream in("tests/data/ten_ab.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string read_id;  string seq = ""; 
    for(int i = 0; i < 2; i++) {
	getline(in, read_id);
	getline(in, seq);
    }
    //
    CreateProfile cp(cab, 2);
    //    
    cp.compute(seq);
    string cdr3_seq = cp.getCDR3();
    TEST_ASSERT_MSG(cdr3_seq == "GCAAGCCCATTGTAGTAGTACCAGCTGCTATCCCTGATG", 
    		    "CDR3 seq incorrect");    
}
