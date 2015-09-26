#include "catch.hpp"
#include <stdlib.h>

#include <utility>
#include <string>
#include <ostream>
#include <sstream>

#include "prettyprint.hpp"
#include "file_io/FastaParser.hpp"
#include "file_io/FastaRefID.hpp"

#include "graphs/MultiKCanonAbGraph.hpp"
#include "graphs/CanonicalAntibodyGraph.hpp"

TEST_CASE("test MultiKCanonAbGraph multi k id ", "[multi-k]") {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    string d_fasta = "data/igh_refs_simple/human_IGHD.fa";
    string j_fasta = "data/igh_refs_simple/human_IGHJ.fa";
    //
    MultiKCanonAbGraph mk_cab(21, 21, 21);
    mk_cab.addVReferences(v_fasta);
    mk_cab.addDReferences(d_fasta);
    mk_cab.addJReferences(j_fasta);
    // read in seq from file
    ifstream in("tests/data/trunc_seq.fa");
    if(!in.is_open()) { REQUIRE(false); }
    string seq;
    getline(in, seq); getline(in, seq);
    //
    vector<int> v = mk_cab.getVPainting("IGHV1-18*01", seq);
    for(int i = 0; i < (int)v.size(); i++) {
	REQUIRE(v[i] == i); //, "error in comparison to reference");
    }
}

TEST_CASE("test MultiKCanonAbGraph multi k test1 ", "[multi-k]") {
    string v_fasta = "data/igh_refs_simple/human_IGHV.fa";
    string d_fasta = "data/igh_refs_simple/human_IGHD.fa";
    string j_fasta = "data/igh_refs_simple/human_IGHJ.fa";
    //
    MultiKCanonAbGraph mk_cab(21, 11, 21);
    mk_cab.addVReferences(v_fasta);
    mk_cab.addDReferences(d_fasta);
    mk_cab.addJReferences(j_fasta);
    //
    CanonicalAntibodyGraph cab(21);
    cab.addVReferences(v_fasta);
    cab.addDReferences(d_fasta);
    cab.addJReferences(j_fasta);       
    // read in seq from file // IGHV1-45*02,IGHD2-2*01,IGHJ3*02
    ifstream in("tests/data/ten_ab.fa");
    if(!in.is_open()) { REQUIRE(false); }
    string read_id;  string seq = ""; 
    for(int i = 0; i < 2; i++) {
	getline(in, read_id);
	getline(in, seq);
    }
    //
    vector<int> mk_v_v = mk_cab.getVPainting("IGHV1-45*02", seq);
    vector<int> mk_d_v = mk_cab.getDPainting("IGHD2-2*01", seq);
    vector<int> mk_j_v = mk_cab.getJPainting("IGHJ3*02", seq);
    //
    vector<int> cab_v_v = cab.getVPainting("IGHV1-45*02", seq);
    vector<int> cab_d_v = cab.getDPainting("IGHD2-2*01", seq);
    vector<int> cab_j_v = cab.getJPainting("IGHJ3*02", seq);
    //
    int index = 5;
    for(int i = 296; i < 308; i++) {
	REQUIRE(mk_d_v[i] == index);
	index++;
    }
    REQUIRE(mk_v_v == cab_v_v);
    REQUIRE(mk_j_v == cab_j_v);
}
