//============================================================================
// Name        : NGramGrap.cpp
// Author      : Yordan
// Version     : 1.1
// Copyright   : Your copyright notice
// Description : Di Gram Graph
//============================================================================

#include <iostream>
#include <map>
#include <string>

using namespace std;

struct NGramGraph{
	string sequence;
	double** matrix;
	map<char, int> baseMapping;

	NGramGraph(string sequenceIn)
	{
		sequence = sequenceIn;
		matrix = new double*[4];
		for(int i = 0; i < 4; ++i)
			matrix[i] = new double[sequenceIn.size()];
		baseMapping['A'] = 0;
		baseMapping['C'] = 1;
		baseMapping['G'] = 2;
		baseMapping['T'] = 3;
	}

	void calculateMatrix()
	{
		//Chargaff's rule
		matrix[baseMapping['A']][0]=0.292;
		matrix[baseMapping['C']][0]=0.206;
		matrix[baseMapping['G']][0]=0.201;
		matrix[baseMapping['T']][0]=0.301;

		//N-gram Graph distribution
		for(unsigned int n = 1;n<sequence.length();n++)
		{
			int currentPosition = n;
			unsigned int windowSize = 200;//Sample language is 401 positions
			int start = 0;
			int end = sequence.size()-1;
			if(!(windowSize*2+1>sequence.size()))
			{
				int windowPrefix = windowSize;
				int windowSuffix = windowSize;

				if(currentPosition-windowSize<0)
				{
					windowPrefix-=windowSize-currentPosition;
					windowSuffix+=windowSize-currentPosition;
				}
				else if(currentPosition+windowSize>=sequence.length())
				{
					windowSuffix-=windowSize-(sequence.length()-1-currentPosition);
					windowPrefix+=windowSize-windowSuffix;

				}

				start = currentPosition-windowPrefix;
				end = currentPosition+windowSuffix;
			}
			int count[4][4]={};
			for(int k = start;k<end;k++)
			{
				++count[baseMapping[sequence[k]]][baseMapping[sequence[k+1]]];
			}

			int sum=0;
			for(int l = 0;l<4;l++){
				sum+=count[baseMapping[sequence[n-1]]][l];
			}
			matrix[baseMapping['A']][n]= (double) count[baseMapping[sequence[n-1]]][baseMapping['A']]/sum; // ((int) (((double) count[baseMapping[sequence[n-1]]][baseMapping['A']]/sum)*1000))/1000.0;
			matrix[baseMapping['C']][n]= (double) count[baseMapping[sequence[n-1]]][baseMapping['C']]/sum; // ((int) (((double) count[baseMapping[sequence[n-1]]][baseMapping['C']]/sum)*1000))/1000.0;
			matrix[baseMapping['G']][n]= (double) count[baseMapping[sequence[n-1]]][baseMapping['G']]/sum; // ((int) (((double) count[baseMapping[sequence[n-1]]][baseMapping['G']]/sum)*1000))/1000.0;
			matrix[baseMapping['T']][n]= (double) count[baseMapping[sequence[n-1]]][baseMapping['T']]/sum; // ((int) (((double) count[baseMapping[sequence[n-1]]][baseMapping['T']]/sum)*1000))/1000.0;
		}
	}
	void printMatrix() {
		for(unsigned int i = 0;i<sequence.size();i++)
			cout << sequence[i] << "\t";
		cout << "\n";

		for (int i = 0; i < 4; ++i)
		{
			for (unsigned int j = 0; j < sequence.size(); ++j)
			{
				std::cout << matrix[i][j] << '\t';
			}
			std::cout << std::endl;
		}
	}
};
//int main() {
//	NGramGraph mg("TTAATTGGTAAATAA");
//	mg.calculateMatrix();
//	mg.printMatrix();
//}
