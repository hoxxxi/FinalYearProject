/*
 * ProbabilityMatrix.cpp
 *
 *  Created on: 10 Mar 2016
 *      Author: yordan
 */

#include "ProbabilityMatrix.h"

ProbabilityMatrix::ProbabilityMatrix(string sequenceIn, string scoreIn) {
	sequence = sequenceIn;
	score = scoreIn;
	size = sequenceIn.size();
	matrix = new double*[size];
	for(int i = 0; i < size; i++)
		matrix[i] = new double[4];
	baseMapping['A'] = 0;
	baseMapping['C'] = 1;
	baseMapping['G'] = 2;
	baseMapping['T'] = 3;
}
ProbabilityMatrix::~ProbabilityMatrix() {
	// TODO Auto-generated destructor stub
}
string ProbabilityMatrix::getSequence() {
	return sequence;
}
string ProbabilityMatrix::getScore() {
	return score;
}
int ProbabilityMatrix::getSize() {
	return size;
}
double** ProbabilityMatrix::getMatrix() {
	return matrix;
}

void ProbabilityMatrix::printMatrix() {
	for(int i = 0; i < size; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			cout<<matrix[i][j]<<'\t';
		}
		cout<<endl;
	}
}

void ProbabilityMatrix::setZeroMatrix() {
	for(int i = 0; i < size; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrix[i][j]=0;
		}
	}
}

void ProbabilityMatrix::applyBigram(int windowSize)
{
	windowSize = (windowSize-1)/2;

	//Chargaff's rule
	matrix[0][baseMapping['A']]=0.292;
	matrix[0][baseMapping['C']]=0.206;
	matrix[0][baseMapping['G']]=0.201;
	matrix[0][baseMapping['T']]=0.301;

	for(int currentPosition = 1; currentPosition < size; currentPosition++)
	{
		//Adjust window
		int start = 0;
		int end = size-1;
		if(!(windowSize*2+1>size))
		{
			int windowPrefix = windowSize;
			int windowSuffix = windowSize;

			if(currentPosition-windowSize<0)
			{
				windowPrefix-=windowSize-currentPosition;
				windowSuffix+=windowSize-currentPosition;
			}
			else if(currentPosition+windowSize>=size)
			{
				windowSuffix-=windowSize-(sequence.length()-1-currentPosition);
				windowPrefix+=windowSize-windowSuffix;

			}

			start = currentPosition-windowPrefix;
			end = currentPosition+windowSuffix;
		}
		//Define combination counter for the window
		int count[4][4]={};
		for(int k = start;k<end;k++)
		{
			++count[baseMapping[sequence[k]]][baseMapping[sequence[k+1]]];
		}
		//Check total count of possible combinations
		int sum = 0;
		for(int l = 0; l < 4; l++){
			sum+= count[baseMapping[sequence[currentPosition-1]]][l];
		}
		matrix[currentPosition][baseMapping['A']]= (double) count[baseMapping[sequence[currentPosition-1]]][baseMapping['A']]/sum;
		matrix[currentPosition][baseMapping['C']]= (double) count[baseMapping[sequence[currentPosition-1]]][baseMapping['C']]/sum;
		matrix[currentPosition][baseMapping['G']]= (double) count[baseMapping[sequence[currentPosition-1]]][baseMapping['G']]/sum;
		matrix[currentPosition][baseMapping['T']]= (double) count[baseMapping[sequence[currentPosition-1]]][baseMapping['T']]/sum;
	}
}

void ProbabilityMatrix::applyQualityScore(int qsCoefficient) {
	for(int i= 0; i < size; i++)
	{
		//Take the int value of the ASCII char as double and use it in the conversion formula
		double value = 1-pow(10.0,(33.0-((double)((int)(score[i]))))/10.0);
		double notValue = 1-value;
		for(int k= 0; k < 4; k++)
		{
			double temp = matrix[i][k];
			if(baseMapping[sequence[i]]==k) {
				matrix[i][k] = (temp+(qsCoefficient*value))/(qsCoefficient+1.0);
			}
			else {
				matrix[i][k] = (temp+(qsCoefficient*notValue/3.0))/(qsCoefficient+1.0);
			}
		}
	}
}

