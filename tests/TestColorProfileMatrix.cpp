#include "TestColorProfileMatrix.h"
#include "../prettyprint.hpp"

TestColorProfileMatrix::TestColorProfileMatrix() { 
    TEST_ADD(TestColorProfileMatrix::test_cp_matrix_scores_1);    
    TEST_ADD(TestColorProfileMatrix::test_cp_matrix_scores_2);

    TEST_ADD(TestColorProfileMatrix::test_cp_matrix_partition_1);
}

TestColorProfileMatrix::~TestColorProfileMatrix() {}

void TestColorProfileMatrix::setup() {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    string d_fasta = "data/igh_refs_simple/human_IGHD.fa";
    string j_fasta = "data/igh_refs_simple/human_IGHJ.fa";
    //
    cab = new CanonicalAntibodyGraph(21);
    cab->addVReferences(v_fasta);
    cab->addDReferences(d_fasta);
    cab->addJReferences(j_fasta);   
}

void TestColorProfileMatrix::tear_down() { delete cab; }

void TestColorProfileMatrix::test_cp_matrix_scores_1() {
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
    vector<vector<int> > v_vect = cp_mat.getVColorProfile();
    vector<int> vscore = cp_mat.getSimpleScoring(v_vect);
    //
    vector<vector<int> > d_vect = cp_mat.getDColorProfile();
    vector<int> dscore = cp_mat.getSimpleScoring(d_vect);
    //
    vector<vector<int> > j_vect = cp_mat.getJColorProfile();
    vector<int> jscore = cp_mat.getSimpleScoring(j_vect);
    //
    vector<int> v_truth = {296, 295, 295, 265, 269, 267, 266, 267, 246, 268, 266, 250, 251, 264, 263, 263, 254, 255, 243, 241, 244, 242, 244, 243, 244, 245, 243, 243, 244, 243, 271, 272, 240, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 198, 199, 0, 202, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 0, 0, 0, 0, 204, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 143, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 141, 139, 140, 140, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 200, 28, 26, 26, 0, 0, 0, 0, 0, 0, 0, 0, 170, 166, 171, 0, 169, 169, 170, 28, 24, 0, 28, 28, 28, 28, 28, 29, 29, 29, 24, 0, 0, 0, 0, 0, 0, 0, 133, 135, 0, 0, 167, 168, 0, 0, 0, 29, 29, 0, 164, 137, 137, 139, 138, 135, 0, 0, 132, 137, 137, 138, 29, 29, 29, 27, 29, 167, 168, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 231, 234, 233, 232};//with tip prop
    //vector<int> v_truth = {296, 295, 295, 206, 252, 234, 249, 208, 226, 210, 225, 106, 107, 247, 246, 241, 213, 214, 243, 241, 244, 242, 244, 243, 244, 245, 243, 243, 244, 243, 254, 255, 218, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 34, 32, 0, 34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 0, 0, 0, 0, 29, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 30, 28, 24, 30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 29, 23, 21, 21, 0, 0, 0, 0, 0, 0, 0, 0, 32, 32, 34, 0, 32, 32, 32, 23, 24, 0, 23, 23, 23, 23, 23, 29, 29, 29, 24, 0, 0, 0, 0, 0, 0, 0, 29, 27, 0, 0, 32, 32, 0, 0, 0, 29, 29, 0, 34, 29, 29, 29, 29, 24, 0, 0, 24, 27, 27, 27, 29, 29, 29, 27, 29, 32, 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 120, 213, 212, 120};
    vector<double> d_truth = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<double> j_truth = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    //
    for(int i = 0; i < (int)v_truth.size(); i++) {
	TEST_ASSERT_MSG(vscore[i] == v_truth[i], "V scores are incorrect");
    }
    for(int i = 0; i < (int)d_truth.size(); i++) {
	TEST_ASSERT_MSG(dscore[i] == d_truth[i], "D scores are incorrect");
    }
    for(int i = 0; i < (int)j_truth.size(); i++) {
	TEST_ASSERT_MSG(jscore[i] == j_truth[i], "J scores are incorrect");
    }
}

void TestColorProfileMatrix::test_cp_matrix_scores_2() {
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
    vector<vector<int> > v_vect = cp_mat.getVColorProfile();
    vector<int> vscore = cp_mat.getSimpleScoring(v_vect);
    //
    vector<vector<int> > d_vect = cp_mat.getDColorProfile();
    vector<int> dscore = cp_mat.getSimpleScoring(d_vect);
    //
    vector<vector<int> > j_vect = cp_mat.getJColorProfile();
    vector<int> jscore = cp_mat.getSimpleScoring(j_vect);    
    //
     vector<double> v_truth = {249, 248, 250, 256, 258, 256, 257, 256, 244, 252, 251, 292, 293, 257, 256, 257, 243, 243, 248, 231, 248, 248, 249, 231, 249, 249, 247, 247, 248, 248, 252, 253, 241, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 220, 221, 220, 219};
    // vector<double> d_truth = {0, 0, 0, 0, 0, 27, 27, 26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};    
    // vector<double> j_truth = {0, 0, 49, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0}; // without early bail on color prop on tips
    //vector<double> v_truth = {107, 107, 107, 223, 224, 222, 223, 223, 217, 224, 233, 292, 293, 239, 238, 239, 219, 219, 230, 213, 229, 229, 231, 213, 231, 230, 228, 228, 229, 230, 234, 235, 212, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 96, 96, 96, 96};
    vector<double> d_truth = {0, 0, 0, 0, 0, 25, 25, 23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    vector<double> j_truth = {0, 0, 49, 50, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    //
    for(int i = 0; i < (int)v_truth.size(); i++) {
    	TEST_ASSERT_MSG(vscore[i] == v_truth[i], "V scores are incorrect");
    }
    for(int i = 0; i < (int)d_truth.size(); i++) {
    	TEST_ASSERT_MSG(dscore[i] == d_truth[i], "D scores are incorrect");
    }
    for(int i = 0; i < (int)j_truth.size(); i++) {
    	TEST_ASSERT_MSG(jscore[i] == j_truth[i], "J scores are incorrect");
    }
}


void TestColorProfileMatrix::test_cp_matrix_partition_1() {
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
    vector<vector<int> > v_vect = cp_mat.getVColorProfile();
    vector<pair<int,int> > v_part=cp_mat.getPartitions(v_vect);
    vector<int> vscore = cp_mat.getSimpleScoring(v_vect);
    //
    vector<vector<int> > d_vect = cp_mat.getDColorProfile();
    vector<pair<int,int> > d_part=cp_mat.getPartitions(d_vect);
    vector<int> dscore = cp_mat.getSimpleScoring(d_vect);
    //
    vector<vector<int> > j_vect = cp_mat.getJColorProfile();
    vector<pair<int,int> > j_part=cp_mat.getPartitions(j_vect);
    vector<int> jscore = cp_mat.getSimpleScoring(j_vect);
    //
    // cout<<endl<<v_part<<endl<<d_part<<endl<<j_part<<endl;
    // cout<<d_part[5]<<endl<<d_part[6]<<endl;
    // cout<<d_vect[5]<<endl;
    // TEST_ASSERT_MSG(d_part[5].first == 296, "D start is incorrect!");
    // TEST_ASSERT_MSG(d_part[5].second == 318, "D end is incorrect!");
    // TEST_ASSERT_MSG(d_part[5].first == 291, "D start is incorrect!");
    // TEST_ASSERT_MSG(d_part[5].second == 321, "D end is incorrect!");
   

    // cout<<j_part[2]<<endl<<j_part[3]<<endl;
    // cout<<j_vect[2]<<endl<<j_vect[3]<<endl;
    TEST_ASSERT_MSG(j_part[2].first == 322, "J start is incorrect!");
    TEST_ASSERT_MSG(j_part[2].second == 371, "J end is incorrect!");
    // TEST_ASSERT_MSG(j_part[2].first == 322, "J start is incorrect!");
    // TEST_ASSERT_MSG(j_part[2].second == 371, "J end is incorrect!");
}

