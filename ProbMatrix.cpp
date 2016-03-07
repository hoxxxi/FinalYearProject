//============================================================================
// Name        : ProbMatrix.cpp
// Author      : Yordan
// Version     : 1.1
// Copyright   : Your copyright notice
// Description : Probability Matrix Combination of quality score and ngram
//============================================================================


#include <fstream>
#include <string>
#include <iostream>
#include <sstream>

#include "NGramGraph.cpp"
#include "QualityScore.cpp"

using namespace std;

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
	void setInitialMatrix()
	{
		NGramGraph graph(sequence);
		graph.calculateMatrix();
		matrix = graph.matrix;
	}
	void calculateDistributedMatrixProbability(int qualityScoreCoefficient)
	{
		QualityScore qs(sequence,score);
		qs.calculateMatrix();

		for(unsigned int i =0;i<sequence.size();i++)
		{
			for(int k = 0;k<4;k++)
			{
				double temp = matrix[k][i];
				matrix[k][i] = (temp+(qualityScoreCoefficient*qs.matrix[k][i]))/(qualityScoreCoefficient+1);
			}
		}
	}
	void printMatrix()
	{
		std::stringstream s;
		for (unsigned int j = 0; j < sequence.size(); ++j)
		{
			for (int i = 0; i < 4; ++i)
			{
				s<<matrix[i][j]<< '\t';
			}
			s<<"\n";
		}
		cout << s.str() <<endl;
		ofstream myfile;
		myfile.open ("output.txt");
		myfile << s.str();
		myfile.close();
	}
	double** getMatrix()
	{
		return matrix;
	}
};
//int main() {
//
//	string sequenceID;
//	string sequence;
//	string score;
//	ifstream file ("left.txt", ios::in);
//	if (file.is_open())
//	{
//		getline (file,sequenceID);
//		getline (file,sequence);
//		getline (file,sequenceID);
//		getline (file,score);
//		file.close();
//	}
//
//	ProbabilityMatrix x(sequence,score);
//	x.setInitialMatrix();
//	x.calculateDistributedMatrixProbability(100);
//	x.printMatrix();
//}
