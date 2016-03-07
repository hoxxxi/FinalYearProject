//============================================================================
// Name        : QualityScore.cpp
// Author      : Yordan
// Version     : 1.1
// Copyright   : Your copyright notice
// Description : Quality score calculator
//============================================================================

#include <iostream>
#include <map>
#include <string>
#include <cmath>

using namespace std;

struct QualityScore
{
	string DNAsequence;
	string QSsequence;
	double** matrix;
	map<char, int> baseMapping;

	QualityScore(string DNAsequenceIn, string QSsequenceIn) {
		DNAsequence = DNAsequenceIn;
		QSsequence = QSsequenceIn;
		matrix = new double*[4];
		for(int i = 0; i < 4; ++i)
			matrix[i] = new double[DNAsequenceIn.size()];
		baseMapping['A']= 0;
		baseMapping['C']= 1;
		baseMapping['G']= 2;
		baseMapping['T']= 3;
	}
	double qsValue(char in)
	{
		return (double)((int)in);
	}
	void calculateMatrix() {
		for(unsigned int i =0;i<QSsequence.length();i++)
		{
			double temp = qsValue(QSsequence[i]);
			double value = 1-pow(10.0,(33.0-temp)/10.0);
			double notValue = 1-value;
			for(int k = 0;k<4;k++)
			{
				if(baseMapping[DNAsequence[i]]==k)
				{
					matrix[k][i] =  value;
				}
				else
					matrix[k][i] = notValue/3;
			}
		}
	}
	void printMatrix()
	{
		for(unsigned int i = 0;i<DNAsequence.size();i++)
			cout << DNAsequence[i] << "\t";
		cout << "\n";

		for (int i = 0; i < 4; ++i)
		{
			for (unsigned int j = 0; j < DNAsequence.size(); ++j)
			{
				cout << matrix[i][j] << '\t';
			}
			cout << std::endl;
		}
	}
};
//int main() {
//	QualityScore mg("TTAATTGGTAAATAA",
//					"#$%&()*,-./0123");
//	mg.calculateMatrix();
//	mg.printMatrix();
//}
