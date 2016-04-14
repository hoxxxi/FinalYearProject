/*
 * ProbabilityMatrix.cpp
 *
 *  Created on: 10 Mar 2016
 *      Author: yordan
 */

#include "ProbabilityMatrix.h"

ProbabilityMatrix::ProbabilityMatrix(string sequenceIn, string scoreIn,  int splitPointIn) {
	splitPoint = splitPointIn;
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
	baseMapping['N'] = 4;
}
ProbabilityMatrix::~ProbabilityMatrix() {
	for (int i = 0; i < size; i++ )
		delete[] matrix[i];
	delete[] matrix;
	//cout<<"Probability matrix destroyed: "<<sequence<<endl;
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
	cout<<endl;
	double roundedValue;
	for(int i = 0; i < size; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			roundedValue =((int) (matrix[i][j]*100))/100.0;
			cout<<roundedValue<<'\t';
		}
		cout<<endl;
	}
}

void ProbabilityMatrix::setZeroMatrix() {
	for(int i = 0; i < size; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			matrix[i][j]=0.0;
		}
	}
}
void ProbabilityMatrix::applyBigram(int windowSize)
{
	if(windowSize<2)
	{
		for(int i = 0; i < size; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					matrix[i][j]=0.25;
				}
			}
	}
	else if(windowSize >= size)
	{
		double count[5][5]={};
		double horizontalSum[5]={};
		for(int k = 0;k<size-1;k++)
		{
			++count[baseMapping[sequence[k]]][baseMapping[sequence[k+1]]];
		}

		for(int edno = 0;edno <5;edno++)
		{
			int temp = 0;
			for(int dve = 0;dve<5;dve++)
			{
				temp+=count[edno][dve];
			}
			horizontalSum[edno]=temp;
		}
		for(int currentPosition = 1; currentPosition < size; currentPosition++)
		{
			matrix[currentPosition][baseMapping['A']]= count[baseMapping[sequence[currentPosition-1]]][baseMapping['A']]/horizontalSum[baseMapping[sequence[currentPosition-1]]] +
					count[baseMapping[sequence[currentPosition-1]]][baseMapping['N']]/horizontalSum[baseMapping[sequence[currentPosition-1]]]/4.0;
			matrix[currentPosition][baseMapping['C']]= count[baseMapping[sequence[currentPosition-1]]][baseMapping['C']]/horizontalSum[baseMapping[sequence[currentPosition-1]]] +
					count[baseMapping[sequence[currentPosition-1]]][baseMapping['N']]/horizontalSum[baseMapping[sequence[currentPosition-1]]]/4.0;
			matrix[currentPosition][baseMapping['G']]= count[baseMapping[sequence[currentPosition-1]]][baseMapping['G']]/horizontalSum[baseMapping[sequence[currentPosition-1]]] +
					count[baseMapping[sequence[currentPosition-1]]][baseMapping['N']]/horizontalSum[baseMapping[sequence[currentPosition-1]]]/4.0;
			matrix[currentPosition][baseMapping['T']]= count[baseMapping[sequence[currentPosition-1]]][baseMapping['T']]/horizontalSum[baseMapping[sequence[currentPosition-1]]] +
					count[baseMapping[sequence[currentPosition-1]]][baseMapping['N']]/horizontalSum[baseMapping[sequence[currentPosition-1]]]/4.0;
		}
		//Initial position is set by an inverse bigram
		double verticalSum = 0;
		for(int v = 0;v<5;v++)
		{
			verticalSum+=count[v][baseMapping[sequence[1]]];
		}

		matrix[0][baseMapping['A']]=count[baseMapping['A']][baseMapping[sequence[1]]]/verticalSum+
				count[baseMapping['N']][baseMapping[sequence[1]]]/verticalSum/4.0;
		matrix[0][baseMapping['C']]=count[baseMapping['C']][baseMapping[sequence[1]]]/verticalSum+
				count[baseMapping['N']][baseMapping[sequence[1]]]/verticalSum/4.0;
		matrix[0][baseMapping['G']]=count[baseMapping['G']][baseMapping[sequence[1]]]/verticalSum+
				count[baseMapping['N']][baseMapping[sequence[1]]]/verticalSum/4.0;
		matrix[0][baseMapping['T']]=count[baseMapping['T']][baseMapping[sequence[1]]]/verticalSum+
				count[baseMapping['N']][baseMapping[sequence[1]]]/verticalSum/4.0;

	}
	else {
		for(int currentPosition = 1; currentPosition < size; currentPosition++)
		{
			//Adjust window
			int windowPrefix = windowSize/2;
			int windowSuffix = windowSize/2-((windowSize-1)%2);
			int start = currentPosition-windowPrefix;
			int end = currentPosition+windowSuffix+1;
			if(currentPosition-windowPrefix<0)
			{
				start+=windowPrefix-currentPosition;
				end+=windowPrefix-currentPosition;
			}
			else if(currentPosition+windowSuffix+1>size)
			{
				start-=currentPosition+windowSuffix+1-size;
				end-=currentPosition+windowSuffix+1-size;
			}
			//Define combination counter for the window
			double count[5][5]={};
			for(int k = start;k<end-1;k++)
			{
				++count[baseMapping[sequence[k]]][baseMapping[sequence[k+1]]];
			}
			//Check total count of possible combinations
			double sum = 0;
			for(int l = 0; l < 5; l++){
				sum+= count[baseMapping[sequence[currentPosition-1]]][l];
			}
			matrix[currentPosition][baseMapping['A']]= count[baseMapping[sequence[currentPosition-1]]][baseMapping['A']]/sum +
					count[baseMapping[sequence[currentPosition-1]]][baseMapping['N']]/sum/4.0;
			matrix[currentPosition][baseMapping['C']]= count[baseMapping[sequence[currentPosition-1]]][baseMapping['C']]/sum +
					count[baseMapping[sequence[currentPosition-1]]][baseMapping['N']]/sum/4.0;
			matrix[currentPosition][baseMapping['G']]= count[baseMapping[sequence[currentPosition-1]]][baseMapping['G']]/sum +
					count[baseMapping[sequence[currentPosition-1]]][baseMapping['N']]/sum/4.0;
			matrix[currentPosition][baseMapping['T']]= count[baseMapping[sequence[currentPosition-1]]][baseMapping['T']]/sum +
					count[baseMapping[sequence[currentPosition-1]]][baseMapping['N']]/sum/4.0;

			//Initial position is set by an inverse bigram
			if(currentPosition==1)
			{
				double verticalSum = 0;
				for(int v = 0;v<5;v++)
				{
					verticalSum+=count[v][baseMapping[sequence[1]]];
				}

				matrix[0][baseMapping['A']]=count[baseMapping['A']][baseMapping[sequence[1]]]/verticalSum+
						count[baseMapping['N']][baseMapping[sequence[1]]]/verticalSum/4.0;
				matrix[0][baseMapping['C']]=count[baseMapping['C']][baseMapping[sequence[1]]]/verticalSum+
						count[baseMapping['N']][baseMapping[sequence[1]]]/verticalSum/4.0;
				matrix[0][baseMapping['G']]=count[baseMapping['G']][baseMapping[sequence[1]]]/verticalSum+
						count[baseMapping['N']][baseMapping[sequence[1]]]/verticalSum/4.0;
				matrix[0][baseMapping['T']]=count[baseMapping['T']][baseMapping[sequence[1]]]/verticalSum+
						count[baseMapping['N']][baseMapping[sequence[1]]]/verticalSum/4.0;
			}
		}
	}
}

void ProbabilityMatrix::applyQualityScore(int qsCoefficient) {
	for(int i= 0; i < size; i++)
	{
		//Take the int value of the ASCII char as double and use it in the conversion formula
		double value = 1-pow(10.0,(33.0-((double)((int)(score[i]))))/10.0);

		if(sequence[i]=='N' || score[i] == (char) 34)
		{
			value = 0.25;
		}
		double notValue = 1-value;
		double maxProb = 0.0;
		char maxValue = 'N';
		for(int k= 0; k < 4; k++)
		{
			double temp = matrix[i][k];
			if(baseMapping[sequence[i]]==k) {
				matrix[i][k] = (temp+(qsCoefficient*value))/(qsCoefficient+1.0);
			}
			else {
				matrix[i][k] = (temp+(qsCoefficient*notValue/3.0))/(qsCoefficient+1.0);
			}
			if(matrix[i][k]>maxProb && matrix[i][k]>0.25)
			{
				maxProb = matrix[i][k];
				switch(k)
				{
				case 0: maxValue = 'A'; break;
				case 1: maxValue = 'C'; break;
				case 2: maxValue = 'G'; break;
				case 3: maxValue = 'T'; break;
				}
			}
		}
		sequence[i] = maxValue;
	}
}
