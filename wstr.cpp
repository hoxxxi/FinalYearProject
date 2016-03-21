/*
 * wstr.cpp
 *
 *  Created on: 21 Mar 2016
 *      Author: yordan
 */

#include "defs.h"

WStr::~WStr(){
	str.clear();
	vector<int>().swap(str);
	prob.clear();
	vector<double>().swap(prob);
	bpt.clear();
	vector<vector<double>>().swap(bpt);
	WP.clear();
	vector<unsigned int>().swap(WP);
	BP.clear();
	vector<unsigned int>().swap(BP);
	FP.clear();
	vector<double>().swap(FP);
}
