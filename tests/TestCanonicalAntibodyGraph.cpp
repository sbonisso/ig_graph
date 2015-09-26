#include "catch.hpp"
#include <stdlib.h>

#include "prettyprint.hpp"

#include "file_io/FastaParser.hpp"
#include "file_io/FastaRefID.hpp"
#include "graphs/ReferenceMap.hpp"
#include "graphs/CanonicalAntibodyGraph.hpp"

#include "TestCABGraph.hpp"

TEST_CASE("CanonicalAntibodyGraph test_paint_1",
	  "[paint1]") {
    
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    //CanonicalAntibodyGraph cab(21);
    TestCABGraph cab(21);
    cab.addVReferences(v_fasta);
    // read in seq from file
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { REQUIRE(false); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    vector<int> v = cab.getVPainting("IGHV1-18*01", seq);
    for(int i = 0; i < (int)v.size(); i++) {
	//TEST_ASSERT_MSG(v[i] == i, "error in comparison to reference");
	REQUIRE(v[i] == i);
    }
}

TEST_CASE("CanonicalAntibodyGraph test_paint_2",
	  "[paint2]") { 
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(v_fasta);
    // read in seq from file
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { REQUIRE(false); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    seq[0] = 'A';
    seq[100] = 'A';
    seq[150] = 'T';
    vector<int> v = cab.getVPainting("IGHV1-18*01", seq);    
    
    int k1 = 21-1;
    REQUIRE(v[0] == -1);
    for(int i = 100-k1; i < 100; i++) { REQUIRE(v[i] > -1); }    
    for(int i = 150-k1; i < 150; i++) { REQUIRE(v[i] > -1); }
}
/**
 *
 */
TEST_CASE("test CanonicalAntibodyGraph test_paint_2_color_propagate",
	  "[paint2,color_propagate]") { 
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(v_fasta);
    // read in seq from file
    ifstream in("tests/data/trunc_seq.fa");
    //if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
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
//TEST_CASE("test CanonicalAntibodyGraph test_paint_2_simple",
//	  "[paint2,simple]")  { 

//     string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
//     string read_fasta = "tests/data/trunc_seq.fa";
//     // CanonicalAntibodyGraph cab(21);
//     // cab.addVReferences(v_fasta);
    
//     // ifstream in("tests/data/trunc_seq.fa");
//     // //if(!in.is_open()) { TEST_FAIL("file not open for reading"); }
//     // string seq;
//     // getline(in, seq); getline(in, seq);
//     // //
//     // seq[0] = 'A';
//     // seq[100] = 'A';
//     // seq[150] = 'T';    
//     // // test with normal color propagation
//     // string ref_id = "IGHV1-18*01";
//     // vector<int> v = cab.getPainting(ref_id, seq, cab.v_kmers_, cab.k_, true);
//     // REQUIRE(v[99] > -1); //"Incorrect value at 99");
//     // REQUIRE(v[100] == -1); //"Incorrect value at 100");
//     // REQUIRE(v[101] > -1);  // "Incorrect value at 101");
//     // // test non-propagation of color
//     // vector<int> v2 = cab.initPainting(ref_id, seq, cab.v_kmers_, cab.k_);
//     // REQUIRE(v2[98] == -1);  //"Incorrect value at 98");
//     // REQUIRE(v2[99] == -1);  // "Incorrect value at 99");
//     // REQUIRE(v2[100] == -1); // "Incorrect value at 100");
//     // REQUIRE(v2[101] > -1);  // "Incorrect value at 101");

//     TestCABGraph cab;
//     vector<int> v = cab.run_test(21, v_fasta, read_fasta, 100, 'A', "IGHV1-18*01");
//     CHECK(v[99] > -1); 
//     CHECK(v[100] == -1);
//     CHECK(v[101] > -1); 
//     //
//     TestCABGraph cab2;
//     vector<int> v2 = cab.run_test(21, v_fasta, read_fasta, 150, 'T', "IGHV1-18*01");
//     CHECK(v2[149] > -1); 
//     CHECK(v2[150] == -1);
//     CHECK(v2[151] > -1); 
// }
/**
 * tests color propagation of tip at start of read
 */
TEST_CASE("CanonicalAntibodyGraph test_paint_2_color_propagate_tip1",
	  "[paint2,color_propagate_tip1]") { 

    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    string read_fasta = "tests/data/trunc_seq.fa";
    //
    TestCABGraph cab;
    vector<int> v = cab.run_test(21, v_fasta, read_fasta, 3, 'A', "IGHV1-18*01");            
    //   
    REQUIRE(v[0] >= 0);
    REQUIRE(v[1] >= 0);
    REQUIRE(v[2] >= 0);
    REQUIRE(v[3] < 0);
}
/**
 * tests color propagation of tip at end of read
 */
TEST_CASE("CanonicalAntibodyGraph test_paint_2_color_propagate_tip2",
	  "[paint2,color_propagate_tip2]") { 

    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    string read_fasta = "tests/data/trunc_seq.fa";
    //
    TestCABGraph cab;
    vector<int> v = cab.run_test(21, v_fasta, read_fasta, 3, 'A', "IGHV1-18*01");    
    //    
    REQUIRE(v[0] >= 0);
    REQUIRE(v[1] >= 0);
    REQUIRE(v[2] >= 0);
    REQUIRE(v[3] < 0);
}
/**
 * tests color propagation of tip at start of J reference
 */
TEST_CASE("CanonicalAntibodyGraph test_paint_1_color_propagate_tip_j_start",
	  "[paint1,color_propagate_tip,j,start]") {
    string j_fasta = "data/igh_refs_simple/human_IGHJ.fa";
    string read_fasta = "tests/data/ten_ab.fa";
    //
    TestCABGraph cab;
    vector<int> v = cab.run_test(21, j_fasta, read_fasta, 333, 'T', "IGHJ3*02");
    //    
    REQUIRE(v[331] == 9);
    REQUIRE(v[332] == 10);
    REQUIRE(v[333] < 0);
}
/**
 * tests color propagation of tip at end of J reference
 */
TEST_CASE("test CanonicalAntibodyGraph test_paint_1_color_propagate_tip_j_end",
	  "[paint1,color_propagate_tip,j,end]") {
    string j_fasta = "data/igh_refs_simple/human_IGHJ.fa";
    string read_fasta = "tests/data/ten_ab.fa";    
    //
    TestCABGraph cab;
    vector<int> v = cab.run_test(21, j_fasta, read_fasta, 365, 'A', "IGHJ3*02");
    // 
    REQUIRE(v[363] >= 0);
    REQUIRE(v[364] >= 0);
    REQUIRE(v[365] < 0);
}
