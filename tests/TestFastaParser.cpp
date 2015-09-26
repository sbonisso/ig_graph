#include "catch.hpp"
#include <stdlib.h>

#include <utility>
#include <string>
#include <utility>
#include <string>
#include <ostream>

#include "prettyprint.hpp"

#include "file_io/FastaParser.hpp"
#include "file_io/FastaRefID.hpp"

using namespace std;

TEST_CASE("test FastaParser parse smAb FASTA ", "[parse]") {
    std::string readFasta = "tests/data/ten_ab.fa";
    FastaRefID<FastaParser> fp(readFasta);
    if(!fp.openFile()) {
	std::cerr<<readFasta<<" not opened properly!\n";
	exit(1);
    }
    
    std::vector<std::string> prefix_seqs{ "CAGGTACAGC",
	    "CAGATGCAGC",
	    "CAGGTGCAGC",
	    "CAGGTGCAGC",
	    "CAGGTGCAGC",
	    "GAGGTGCAGC",
	    "CAGCTGCAGC",
	    "CAGGTGCAGC",
	    "CAGGTGCAGC",
	    "GAGGTGCAGC"};
    int i = 0;
    while(fp.hasNextSequence()) {
    	pair<string,string> entry = fp.getNextEntry();
    	string readID = fp.getCurrIDLine().substr(1); 
    	string seq = entry.second;
	REQUIRE( prefix_seqs[i] == seq.substr(0,10) );
	i++;
    }
    fp.closeFile();
}

TEST_CASE("test FastaParser parse file ext", "[parse]") {
    std::string readFasta = "tests/data/ten_ab_ext_ids.fa";
    std::string readFasta_o = "tests/data/ten_ab.fa";
    FastaRefID<FastaParser> fp(readFasta);
    FastaRefID<FastaParser> fp_o(readFasta_o);
    if(!fp.openFile()) {
	std::cerr<<readFasta<<" not opened properly!\n";
	exit(1);
    }
    if(!fp_o.openFile()) {
	std::cerr<<readFasta<<" not opened properly!\n";
	exit(1);
    }
    
    int i = 0;
    while(fp.hasNextSequence()) {
    	pair<string,string> entry = fp.getNextEntry();
	pair<string,string> entry_o = fp_o.getNextEntry();
	string readID = fp.getID();
    	string seq = entry.second;
	REQUIRE(readID == fp_o.getCurrIDLine().substr(1));
	i++;
    }
    fp.closeFile();
}
