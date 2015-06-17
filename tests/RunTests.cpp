#include <iostream>
#include <ostream>

#include "TestReferenceMap.hpp"
#include "TestCanonicalAntibodyGraph.hpp"
#include "TestMultiKCanonAbGraph.hpp"

#include "TestMutationProbSingleton.hpp"
#include "TestMutationNBModel.hpp"
#include "TestColorProfileMatrix.hpp"
#include "TestCreateProfile.hpp"

#include "TestDClassify.hpp"
#include "TestEncoding.hpp"

using namespace std;

int main(int argc, char **argv) {

    TestReferenceMap sts;
    Test::TextOutput output(Test::TextOutput::Verbose);
        
    // how to add additional tests to be run at once
//	sts.add(auto_ptr<Test::Suite>(new TestFastaParser));
//	sts.add(auto_ptr<Test::Suite>(new TestColoredAntibodyGraph));
//	sts.add(auto_ptr<Test::Suite>(new TestColorProfileScoring));
//    sts.add(auto_ptr<Test::Suite>(new TestReferenceMap));
    
    sts.add(auto_ptr<Test::Suite>(new TestMutationProbSingleton));
    sts.add(auto_ptr<Test::Suite>(new TestCanonicalAntibodyGraph));
    sts.add(auto_ptr<Test::Suite>(new TestMultiKCanonAbGraph));
    sts.add(auto_ptr<Test::Suite>(new TestColorProfileMatrix));
    sts.add(auto_ptr<Test::Suite>(new TestCreateProfile));
    sts.add(auto_ptr<Test::Suite>(new TestDClassify));
    sts.add(auto_ptr<Test::Suite>(new TestEncoding));

    sts.add(auto_ptr<Test::Suite>(new TestMutationNBModel));
    
    return sts.run(output, false); // Note the 'false' parameter
}
