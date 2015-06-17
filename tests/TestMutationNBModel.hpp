#ifndef TESTMUTATIONNBMODEL_H_
#define TESTMUTATIONNBMODEL_H_

#include "cpptest.h"
#include <utility>
#include <string>
#include <ostream>
#include <stdlib.h>
#include <stdio.h>

#include "file_io/FastaParser.hpp"

#include "seq_utils/MutationNBProbabilities.h"
#include "seq_utils/MutationNBModel.hpp"
#include "seq_utils/Encoding.hpp"

using namespace std;

class TestMutationNBModel : public Test::Suite
{
public:
    TestMutationNBModel()
    {
TEST_ADD(TestMutationNBModel::test_model_similarity_along_pos);
TEST_ADD(TestMutationNBModel::test_model_similarity_across_lmers);
}
    
private:

void test_model_similarity_along_pos();
void test_model_similarity_across_lmers();
};

#endif /* TESTMUTATIONNBMODEL_H_ */
