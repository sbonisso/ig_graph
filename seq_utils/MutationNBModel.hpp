#ifndef MUTATIONNBMODEL_H_
#define MUTATIONNBMODEL_H_

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>

#include "seq_utils/Encoding.hpp"

#include <cereal/cereal.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/external/rapidjson/filestream.h>

#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/map.hpp>

#ifdef DEBUG
#define MUTNBMOD_DEBUG_PRINT(EXPR) DEBUG_PRINT("[MUT_NB_MOD]", EXPR)
#else
#define MUTNBMOD_DEBUG_PRINT(x) do {} while (0)
#endif

class MutationNBModel {
private:
    MutationNBModel();
    ~MutationNBModel();
	
    // hashes to hold class conditional probabilities
    std::map<std::string,std::map<std::string, int> > lmerCounts;
    std::map<std::string,std::map<int, int> > posCounts;
    std::map<std::string,double> priorProbs;
    
    std::map<std::string,std::map<std::string, double> > lmerProbs;
    std::map<std::string,std::map<int, double> > posProbs;
    
    std::string paramDir;
    int lmer_len;
    
    // to read in data from files
    void readInData();
    void fillLMerHash();
    void fillPosHash();
    void fillPriors();

    // to control singleton
    static bool _instanceExists;
    static MutationNBModel* _nbInstance;
	
    std::vector<double> pos_mut_;
    std::vector<double> pos_nomut_;
    std::vector<double> lmer_mut_;
    std::vector<double> lmer_nomut_;
    double prior_mut_;
    double prior_nomut_;
	
public:

    static MutationNBModel& getInstance();
    void setParamDir(std::string paramDirStr);
    int getLMerLen();

    double getProb(std::string lmer, int pos);
    double getScore(std::string lmer, int pos);

    double getProb(unsigned int lmer_i, int pos);
    int getMaxVal() { return max_val_; }
    
    int max_val_;
};

#endif /* MUTATIONNBMODEL_H_ */
