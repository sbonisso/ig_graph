#include "MutationNBModel.hpp"

bool MutationNBModel::_instanceExists;
MutationNBModel* MutationNBModel::_nbInstance;

MutationNBModel& MutationNBModel::getInstance() {
    if(!_instanceExists) {
	_nbInstance = new MutationNBModel();
	_instanceExists = true;
    }
    return (*_nbInstance);
}

MutationNBModel::MutationNBModel() { lmer_len = 0; }

MutationNBModel::~MutationNBModel() { this->_instanceExists = false; }
/**
 *
 */
void MutationNBModel::setParamDir(std::string paramDirStr) {
    this->paramDir = paramDirStr;
    this->readInData();
    MUTNBMOD_DEBUG_PRINT("SET PARAM DIR:\t"<<paramDir);
    MUTNBMOD_DEBUG_PRINT("AGAA:\t"<<this->lmerProbs["match"]["AGAA"]
			 <<"\t"<<this->lmerProbs["mutation"]["AGAA"]);
    MUTNBMOD_DEBUG_PRINT("91:\t"<<this->posProbs["match"][91]
			 <<"\t"<< this->posProbs["mutation"][91]);
    MUTNBMOD_DEBUG_PRINT("priors\t"<<this->priorProbs["match"]
			 <<"\t"<<priorProbs["mutation"]);
}
/**
 *
 */
void MutationNBModel::readInData() {	
    std::ifstream is(paramDir);
    cereal::BinaryInputArchive ar(is);
    ar(lmerCounts, posCounts, priorProbs);
    is.close();
    
    lmer_len = 4;
    
    this->fillLMerHash();
    this->fillPosHash();
    this->fillPriors();
}
/**
 *
 */
double MutationNBModel::getProb(std::string lmer, int pos) {
    double mutProb = posProbs["mutation"][pos]*priorProbs["mutation"];
    double matchProb = posProbs["match"][pos]*priorProbs["match"];
    
    // is '....' l-mer, signalling a position only prob
    if(lmer[0]=='.') {	}
    else {
	double lmerMut =lmerProbs["mutation"][lmer];
	double lmerMatch = lmerProbs["match"][lmer];
	
	// if l-mer is not present, leave it out
	if(!(lmerMut == 0 && lmerMatch == 0)) {
	    mutProb *= lmerMut;
	    matchProb *= lmerMatch;
	}
    }

    double currProb = mutProb;
    double denomProb = mutProb + matchProb;
    
    return (currProb / denomProb);
}

double MutationNBModel::getProb(unsigned int lmer_i, int pos) {
    double mutProb = pos_mut_[pos]*prior_mut_;
    double matchProb = pos_nomut_[pos]*prior_nomut_;
    
    // is '....' l-mer, signalling a position only prob
    if(lmer_i == 999999) {}
    else {
	double lmerMut = lmer_mut_[lmer_i];
	double lmerMatch = lmer_nomut_[lmer_i];
	// if l-mer is not present, leave it out
	if(!(lmerMut == 0 && lmerMatch == 0)) {
	    mutProb *= lmerMut;
	    matchProb *= lmerMatch;
	}
    }

    double currProb = mutProb;
    double denomProb = mutProb + matchProb;

    return (currProb / denomProb);
}


double MutationNBModel::getScore(std::string lmer, int pos) { return this->getProb(lmer, pos); }

int MutationNBModel::getLMerLen() { return this->lmer_len; }
/**
 *
 */
void MutationNBModel::fillPriors() {
    prior_nomut_ = priorProbs["match"];
    prior_mut_ = priorProbs["mutation"];
}
/**
 *
 */
void MutationNBModel::fillLMerHash() {
    int totalMut = 0;
    int totalMatch=0;
    for(auto const &vals : lmerCounts["match"]) { lmer_len = vals.first.size(); }
    for(auto const &vals : lmerCounts["match"]) { totalMatch += vals.second; }
    for(auto const &vals : lmerCounts["mutation"]) { totalMut += vals.second; }
    
    std::string max_lmer;
    max_lmer.assign(lmer_len, 'T');
    max_val_ = Encoding::encode_kmer(max_lmer);
    
    lmer_nomut_.resize(lmerCounts["match"].size());
    lmer_mut_.resize(lmerCounts["mutation"].size());
    for(auto const &vals : lmerCounts["match"]) {
	std::string lmer = vals.first;
	int numMatch = vals.second;
	lmerProbs["match"][lmer] = ((double)numMatch/((double)totalMatch));
	unsigned int lmer_i = Encoding::encode_kmer(lmer);
	lmer_nomut_[lmer_i] = ((double)numMatch/((double)totalMatch));
    }
    for(auto const &vals : lmerCounts["mutation"]) {
	std::string lmer = vals.first;
	int numMut = vals.second;
	lmerProbs["mutation"][lmer] = ((double)numMut/((double)totalMut));
	unsigned int lmer_i = Encoding::encode_kmer(lmer);
	lmer_mut_[lmer_i] = ((double)numMut/((double)totalMut));
    }

}
/**
 *
 */
void MutationNBModel::fillPosHash() { 
    int totalMut = 0;
    int totalMatch=0;
    for(auto const& vals : posCounts["match"]) { totalMatch += vals.second; }
    for(auto const& vals : posCounts["mutation"]) { totalMut += vals.second; }
    
    pos_nomut_.resize(400,0);
    pos_mut_.resize(400,0);
    for(auto const &vals : posCounts["match"]) {
	int pos = vals.first;
	int numMatch = vals.second;
	posProbs["match"][pos] = ((double)numMatch/((double)totalMatch));
	pos_nomut_[pos] = ((double)numMatch/((double)totalMatch));
    }
    for(auto const &vals : posCounts["mutation"]) {
	int pos = vals.first;
	int numMut = vals.second;
	posProbs["mutation"][pos] = ((double)numMut/((double)totalMut));
	pos_mut_[pos] = ((double)numMut/((double)totalMut));
    }
}
