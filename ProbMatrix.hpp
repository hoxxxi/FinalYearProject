/*
 * ProbabilityMatrix.hpp
 *
 *  Created on: 8 Mar 2016
 *      Author: yordan
 */

#include "ProbMatrix.cpp"

#include <string>

#ifndef PROBABILITYMATRIX_HPP_
#define PROBABILITYMATRIX_HPP_

struct ProbabilityMatrix
{
	string sequence;
	string score;
	double** matrix;
	ProbabilityMatrix(string sequenceIn, string scoreIn)
	{
		sequence = sequenceIn;
		score = scoreIn;
		matrix = new double*[4];
		for(int i = 0; i < 4; ++i)
			matrix[i] = new double[sequence.size()];
	}
	void setInitialMatrix();
	void calculateDistributedMatrixProbability(int qualityScoreCoefficient);
	void printMatrix();
	double** getMatrix();
};

ProbabilityMatrix joinedMatrix(ProbabilityMatrix* left,ProbabilityMatrix* right, int overlapping);
#endif /* PROBABILITYMATRIX_HPP_ */
