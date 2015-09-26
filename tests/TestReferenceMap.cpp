#include "catch.hpp"
#include <stdlib.h>

#include "prettyprint.hpp"
#include "graphs/ReferenceMap.hpp"

TEST_CASE("test ReferenceMap on basic sequences", 
	  "[test1]") {
    
    ReferenceMap rm;
    bool retval = rm.addReference("V1", "ATCGAA");
    REQUIRE(retval == true);
    
    rm.addKmerToReference("V1", "ATCG", 1);
    rm.addKmerToReference("V1", "TCGA", 2);
    rm.addKmerToReference("V1", "CGAA", 3); 
    
    int t1 = rm.getIndex("V1", "TCGA");
    REQUIRE(t1 == 2);
    
    int f1 = rm.getIndex("V1", "AAAA");
    REQUIRE(f1 == -1);
}
