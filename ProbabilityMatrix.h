/*
 *  Created on: 29 Feb 2016
 *  Author: Yordan Petrov Yordanov
 */

#include <map>
#include <cmath>
#include <string>
#include <sstream>
#include <iostream>

#ifndef PROBABILITYMATRIX_H_
#define PROBABILITYMATRIX_H_

using namespace std;

class ProbabilityMatrix {
	int splitPoint;
	string sequence;
	string score;
	int size;
	double** matrix;
	map<char, int> baseMapping;
public:
	ProbabilityMatrix(string sequenceIn, string scoreIn, int splitPoint);
	virtual ~ProbabilityMatrix();
	string getSequence();
	string getScore();
	int getSize();
	double** getMatrix();
	void printMatrix();
	void setZeroMatrix();
	void applyBigram(int windowSize);
	void applyQualityScore(int qsCoefficient);
};

#endif /* PROBABILITYMATRIX_H_ */
