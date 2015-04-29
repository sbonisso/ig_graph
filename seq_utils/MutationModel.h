/*
 * MutationModel.h
 *
 *  Created on: Dec 9, 2013
 *      Author: sbonisso
 */

#ifndef MUTATIONMODEL_H_
#define MUTATIONMODEL_H_

#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>

#include "opencv2/core/core.hpp"
#include "opencv2/ml/ml.hpp"

using namespace cv;

class MutationModel {
private:

	MutationModel() {}
	MutationModel(std::string filePath) {
		rtrees.load(filePath.c_str());	// load model from XML file
		// create mapping that works for the trained model
		lblMap["A"] = 4;
		lblMap["C"] = 2;
		lblMap["G"] = 3;
		lblMap["T"] = 1;
		lblMap["match"] = 5;
		lblMap["mutation"] = 6;
		tstX.create(1, 9, CV_32F);
	}

	CvRTrees rtrees;	// RF model to be loaded
	//Mat tstX(1,9, CV_32F);
	Mat tstX;
	map<string,int> lblMap;

	MutationModel(MutationModel const&);
	void operator=(MutationModel const&);

public:

	static MutationModel& getInstance() {
		static MutationModel *mP = new MutationModel("/tmp/rf_model.xml");
		return *mP;
	}

	double getProb(std::string lmer, int pos) {
		for(int i = 0; i < (int)lmer.size(); i++) {
			string ch(lmer.substr(i,1));
			tstX.at<float>(0,i) = lblMap[ch];
		}
		tstX.at<float>(0,8) = (float)pos;
		float predVal = rtrees.predict_prob(tstX);
		//return (predVal == 0) ? 0.0001 : predVal;
		if(predVal == 0) return 0.0001;
		else if(predVal == 1) return 1-0.0001;
		else return predVal;
	}

	double getScore(std::string lmer, int pos) { return this->getProb(lmer, pos); }
};


#endif /* MUTATIONMODEL_H_ */
