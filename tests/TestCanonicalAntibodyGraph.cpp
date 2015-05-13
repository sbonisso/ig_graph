#include "TestCanonicalAntibodyGraph.h"
#include "../prettyprint.hpp"

TestCanonicalAntibodyGraph::TestCanonicalAntibodyGraph() { 
    TEST_ADD(TestCanonicalAntibodyGraph::test_paint_1);
    TEST_ADD(TestCanonicalAntibodyGraph::test_paint_2);
    TEST_ADD(TestCanonicalAntibodyGraph::test_paint_2_color_propagate);
    TEST_ADD(TestCanonicalAntibodyGraph::test_paint_2_simple);
    //
    TEST_ADD(TestCanonicalAntibodyGraph::test_paint_1_color_propagate_tip_j_start);
    TEST_ADD(TestCanonicalAntibodyGraph::test_paint_1_color_propagate_tip_j_end);
    TEST_ADD(TestCanonicalAntibodyGraph::test_paint_2_color_propagate_tip1);
    TEST_ADD(TestCanonicalAntibodyGraph::test_paint_2_color_propagate_tip2);
}

TestCanonicalAntibodyGraph::~TestCanonicalAntibodyGraph() {}

void TestCanonicalAntibodyGraph::setup() {}

void TestCanonicalAntibodyGraph::test_paint_1() {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(v_fasta);
    // read in seq from file
    //ifstream in("tests/data/trunc_seq.txt");
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    vector<int> v = cab.getVPainting("IGHV1-18*01", seq);
    for(int i = 0; i < (int)v.size(); i++) {
	TEST_ASSERT_MSG(v[i] == i, "error in comparison to reference");
    }
}

void TestCanonicalAntibodyGraph::test_paint_2() {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(v_fasta);
    // read in seq from file
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    seq[0] = 'A';
    seq[100] = 'A';
    seq[150] = 'T';
    vector<int> v = cab.getVPainting("IGHV1-18*01", seq);    
    
    int k1 = 21-1;
    TEST_ASSERT(v[0] == -1);
    for(int i = 100-k1; i < 100; i++) { TEST_ASSERT(v[i] > -1); }    
    for(int i = 150-k1; i < 150; i++) { TEST_ASSERT(v[i] > -1); }
}
/**
 *
 */
void TestCanonicalAntibodyGraph::test_paint_2_color_propagate() {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(v_fasta);
    // read in seq from file
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);
        
    seq[0] = 'A';
    seq[100] = 'A';
    seq[150] = 'T';
    vector<int> v = cab.getVPainting("IGHV1-18*01", seq);
    // int k1 = 21-1;
    
    // TEST_ASSERT(v[0] == -1); 
    
    // for(int i = 100-k1; i < 100+1; i++) {
    // 	if(i == 100) { TEST_ASSERT(v[i] == -1); }
    // 	else { TEST_ASSERT(v[i] == i); }
    // }
    
    // for(int i = 150-k1; i < 150+1; i++) { 
    // 	if(i == 150) { TEST_ASSERT(v[i] == -1); }
    // 	else { TEST_ASSERT(v[i] == i); }
    // }
}
/**
 *
 */
void TestCanonicalAntibodyGraph::test_paint_2_simple() {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(v_fasta);
    
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    seq[0] = 'A';
    seq[100] = 'A';
    seq[150] = 'T';    
    // test with normal color propagation
    string ref_id = "IGHV1-18*01";
    vector<int> v = cab.getPainting(ref_id, seq, cab.v_kmers_, cab.k_, true);
    TEST_ASSERT_MSG(v[99] > -1, "Incorrect value at 99");
    TEST_ASSERT_MSG(v[100] == -1, "Incorrect value at 100");
    TEST_ASSERT_MSG(v[101] > -1, "Incorrect value at 101");
    // test non-propagation of color
    vector<int> v2 = cab.initPainting(ref_id, seq, cab.v_kmers_, cab.k_);
    TEST_ASSERT_MSG(v2[98] == -1, "Incorrect value at 98");
    TEST_ASSERT_MSG(v2[99] == -1, "Incorrect value at 99");
    TEST_ASSERT_MSG(v2[100] == -1, "Incorrect value at 100");
    TEST_ASSERT_MSG(v2[101] > -1, "Incorrect value at 101");
    
}
/**
 * tests color propagation of tip at start of read
 */
void TestCanonicalAntibodyGraph::test_paint_2_color_propagate_tip1() {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(v_fasta);
    
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    seq[3] = 'A';
    //
    string ref_id = "IGHV1-18*01";
    vector<int> v = cab.initPainting(ref_id, seq, cab.v_kmers_, cab.k_);
    cab.propagateColorTip(ref_id, seq, cab.v_kmers_, v, true);
    //
    TEST_ASSERT_MSG(v[0] >= 0, "Incorrect value at 0");
    TEST_ASSERT_MSG(v[1] >= 0, "Incorrect value at 1");
    TEST_ASSERT_MSG(v[2] >= 0, "Incorrect value at 2");
    TEST_ASSERT_MSG(v[3] < 0, "Incorrect value at 3");
}
/**
 * tests color propagation of tip at end of read
 */
void TestCanonicalAntibodyGraph::test_paint_2_color_propagate_tip2() {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(v_fasta);
    
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    seq[3] = 'A';
    //
    string ref_id = "IGHV1-18*01";
    vector<int> v = cab.initPainting(ref_id, seq, cab.v_kmers_, cab.k_);
    cab.propagateColorTip(ref_id, seq, cab.v_kmers_, v, true);
    //
    TEST_ASSERT_MSG(v[0] >= 0, "Incorrect value at 0");
    TEST_ASSERT_MSG(v[1] >= 0, "Incorrect value at 1");
    TEST_ASSERT_MSG(v[2] >= 0, "Incorrect value at 2");
    TEST_ASSERT_MSG(v[3] < 0, "Incorrect value at 3");
}
/**
 * tests color propagation of tip at start of J reference
 */
void TestCanonicalAntibodyGraph::test_paint_1_color_propagate_tip_j_start() {
    //string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    string j_fasta = "data/igh_refs_simple/human_IGHJ.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(j_fasta);
    //
    ifstream in("tests/data/ten_ab.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string read_id;  string seq = ""; 
    for(int i = 0; i < 2; i++) {
	getline(in, read_id);
	getline(in, seq);
    }
    //
    seq[333] = 'T';    
    //
    string ref_id = "IGHJ3*02";
    vector<int> v = cab.initPainting(ref_id, seq, cab.v_kmers_, cab.k_);
    cab.propagateColorTip(ref_id, seq, cab.v_kmers_, v, true);
    //
    TEST_ASSERT_MSG(v[331] == 9, "Incorrect value at 331");
    TEST_ASSERT_MSG(v[332] == 10, "Incorrect value at 332");
    TEST_ASSERT_MSG(v[333] < 0, "Incorrect value at 333");
}
/**
 * tests color propagation of tip at end of J reference
 */
void TestCanonicalAntibodyGraph::test_paint_1_color_propagate_tip_j_end() {
    //string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    string j_fasta = "data/igh_refs_simple/human_IGHJ.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(j_fasta);
    //
    ifstream in("tests/data/ten_ab.fa");
    if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
    string read_id;  string seq = ""; 
    for(int i = 0; i < 2; i++) {
	getline(in, read_id);
	getline(in, seq);
    }
    //
    seq[365] = 'A';
    //
    string ref_id = "IGHJ3*02";
    vector<int> v = cab.initPainting(ref_id, seq, cab.v_kmers_, cab.k_);
    cab.propagateColorTip(ref_id, seq, cab.v_kmers_, v, true);
    //
    TEST_ASSERT_MSG(v[363] >= 0, "Incorrect value at 363");
    TEST_ASSERT_MSG(v[364] >= 0, "Incorrect value at 364");
    TEST_ASSERT_MSG(v[365] < 0, "Incorrect value at 365");
}
