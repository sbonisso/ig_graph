#include <iostream>
#include <ostream>

// #include "TestColoredAntibodyGraph.h"
// #include "TestColorProfileScoring.h"
//#include "TestPartitionedReads.h"

#include "TestReferenceMap.h"
#include "TestCanonicalAntibodyGraph.h"
#include "TestMultiKCanonAbGraph.h"

#include "TestMutationProbSingleton.h"
#include "TestColorProfileMatrix.h"
#include "TestCreateProfile.h"

using namespace std;

int main(int argc, char **argv) {

//	TestColoredAntibodyGraph sts;
//	TestMutationProbSingleton sts;
//	TestColorProfileScoring sts;
//    TestPartitionedReads sts;
    TestReferenceMap sts;
    Test::TextOutput output(Test::TextOutput::Verbose);
        
    // how to add additional tests to be run at once
//	sts.add(auto_ptr<Test::Suite>(new TestFastaParser));
//	sts.add(auto_ptr<Test::Suite>(new TestColoredAntibodyGraph));
//	sts.add(auto_ptr<Test::Suite>(new TestColorProfileScoring));
//    sts.add(auto_ptr<Test::Suite>(new TestReferenceMap));
    
    //sts.add(auto_ptr<Test::Suite>(new TestMutationProbSingleton));
    sts.add(auto_ptr<Test::Suite>(new TestCanonicalAntibodyGraph));
    sts.add(auto_ptr<Test::Suite>(new TestMultiKCanonAbGraph));
    sts.add(auto_ptr<Test::Suite>(new TestColorProfileMatrix));
    sts.add(auto_ptr<Test::Suite>(new TestCreateProfile));
    
return sts.run(output, false); // Note the 'false' parameter
//return sts.run(output, true);
}
