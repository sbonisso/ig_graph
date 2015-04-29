#include "TestCreateProfile.h"
#include "../prettyprint.hpp"

TestCreateProfile::TestCreateProfile() { 
    TEST_ADD(TestCreateProfile::test_cp_preds_1);
    TEST_ADD(TestCreateProfile::test_cp_scores_1);
    
    TEST_ADD(TestCreateProfile::test_cp_preds_2);
    TEST_ADD(TestCreateProfile::test_cp_scores_2);

    TEST_ADD(TestCreateProfile::test_max_n);
}

TestCreateProfile::~TestCreateProfile() {}

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

void TestCreateProfile::test_cp_preds_1() {
    //ifstream in("tests/data/trunc_seq.txt");
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    CreateProfile cp(cab, 2);
    //
    ColorProfileMatrix cp_mat = cp.getColorProfile(seq);
    //
    vector<string> vpred = cp.getPredictedV();
    vector<string> dpred = cp.getPredictedD();
    vector<string> jpred = cp.getPredictedJ();
    //
    vector<string> vtruth = {"IGHV1-18*01", "IGHV1-18*03"};
    vector<string> dtruth = {"IGHD1-1*01", "IGHD1-20*01"};
    vector<string> jtruth = {"IGHJ1*01", "IGHJ2*01"};
    //
    for(int i = 0; i < 2; i++) { TEST_ASSERT(vpred[0] == vtruth[0]); }
    for(int i = 0; i < 2; i++) { TEST_ASSERT(dpred[0] == dtruth[0]); }
    for(int i = 0; i < 2; i++) { TEST_ASSERT(jpred[0] == jtruth[0]); }
}

void TestCreateProfile::test_cp_scores_1() {
    //ifstream in("tests/data/trunc_seq.txt");
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);    
    //
    CreateProfile cp(cab);
    //
    ColorProfileMatrix cp_mat = cp.getColorProfile(seq);
    //
    vector<double> vscore = cp.getVScores();
    vector<double> dscore = cp.getDScores();
    vector<double> jscore = cp.getJScores();
    //
    //vector<int> v_truth = {275, 274, 274, 181, 226, 208, 223, 183, 201, 185, 205, 83, 84, 221, 220, 216, 191, 192, 222, 220, 223, 221, 223, 222, 223, 224, 222, 222, 223, 222, 229, 230, 198, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 6, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 6, 3, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 3, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 7, 7, 8, 0, 7, 7, 7, 3, 3, 0, 3, 3, 3, 3, 3, 8, 8, 8, 3, 0, 0, 0, 0, 0, 0, 0, 8, 7, 0, 0, 6, 6, 0, 0, 0, 8, 8, 0, 8, 8, 8, 8, 8, 3, 0, 0, 3, 6, 7, 7, 8, 8, 8, 6, 8, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 95, 187, 186, 95}; // no tip prop
     vector<int> v_truth = {296, 295, 295, 265, 269, 267, 266, 267, 246, 268, 266, 250, 251, 264, 263, 263, 254, 255, 243, 241, 244, 242, 244, 243, 244, 245, 243, 243, 244, 243, 271, 272, 240, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 198, 199, 0, 202, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0, 0, 0, 204, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 143, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 141, 139, 140, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 200, 28, 26, 26, 0, 0, 0, 0, 0, 0, 0, 0, 170, 166, 171, 0, 169, 169, 170, 28, 24, 0, 28, 28, 28, 28, 28, 29, 29, 29, 24, 0, 0, 0, 0, 0, 0, 0, 133, 135, 0, 0, 167, 168, 0, 0, 0, 29, 29, 0, 164, 137, 137, 139, 138, 135, 0, 0, 132, 137, 137, 138, 29, 29, 29, 27, 29, 167, 168, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 231, 234, 233, 232};//with tip prop
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
    ColorProfileMatrix cp_mat = cp.getColorProfile(seq);
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
    ColorProfileMatrix cp_mat = cp.getColorProfile(seq);
    //
    vector<double> vscore = cp.getVScores();
    vector<double> dscore = cp.getDScores();
    vector<double> jscore = cp.getJScores();
    //
    // vector<double> v_truth = {84, 84, 84, 200, 201, 199, 200, 200, 194, 201, 210, 272, 273, 216, 215, 216, 199, 199, 207, 193, 206, 206, 208, 193, 208, 207, 205, 205, 206, 207, 211, 212, 192, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 73, 73, 73, 73};
    // vector<double> d_truth = {0, 0, 0, 0, 0, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    // vector<double> j_truth = {0, 0, 15, 29, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<double> v_truth = {249, 248, 250, 256, 258, 256, 257, 256, 244, 252, 251, 292, 293, 257, 256, 257, 243, 243, 248, 231, 248, 248, 249, 231, 249, 249, 247, 247, 248, 248, 252, 253, 241, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 220, 221, 220, 219};
    vector<double> d_truth = {0, 0, 0, 0, 0, 27, 27, 26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};    
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

void TestCreateProfile::test_max_n() {
    vector<double> test_scores = {39, 0, 0, 0, 43, 42, 43, 48, 46, 0, 0, 0, 0};
    
    CreateProfile cp(cab);
    vector<int> max_index = cp.getNMaxIndex(3, test_scores);
    
    TEST_ASSERT(max_index[0] == 7);
    TEST_ASSERT(max_index[1] == 8);
    TEST_ASSERT(max_index[2] == 4);    
}
