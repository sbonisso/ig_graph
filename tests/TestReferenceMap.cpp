#include "TestReferenceMap.h"
#include "../prettyprint.hpp"

TestReferenceMap::TestReferenceMap() {
    TEST_ADD(TestReferenceMap::test_ref_map_1);
}

TestReferenceMap::~TestReferenceMap() {}

void TestReferenceMap::setup() {}

void TestReferenceMap::test_ref_map_1() {
    ReferenceMap rm;
    bool retval = rm.addReference("V1", "ATCGAA");
    TEST_ASSERT(retval == true);
    
    rm.addKmerToReference("V1", "ATCG", 1);
    rm.addKmerToReference("V1", "TCGA", 2);
    rm.addKmerToReference("V1", "CGAA", 3); 
    
    int t1 = rm.getIndex("V1", "TCGA");
    TEST_ASSERT(t1 == 2);
    
    int f1 = rm.getIndex("V1", "AAAA");
    TEST_ASSERT(f1 == -1);
}

