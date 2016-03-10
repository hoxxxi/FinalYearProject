#include <iostream>
#include <algorithm>
#include <fstream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>

#include "defs.h"
#include "global.h"
#include "BorderTable.cpp"
#include "StringMatching.cpp"

using namespace std;

int main (int argc, char **argv)
{
	string alphabet = DNA;
	unsigned int sigma = alphabet.size();
	int mod = 0;
	string weightedStringFile = "/home/yordan/workspace/FinalYearProject/data/text0.txt";
	double z = 100000;
	double ** y;			//weighted string
	unsigned int n;			//length of y

	clock_t start;
	clock_t finish;

	string line;
	int lineCounter;
	vector<ProbabilityMatrix> leftProbabilityMatrixes;
	vector<ProbabilityMatrix> rightProbabilityMatrixes;

	string left_sequenceID;
	string left_sequence;
	string left_score;

	ifstream left_file ("data/left.txt", ios::in);
	if (left_file.is_open())
	{
		lineCounter = 0;
		while (getline(left_file, line))
		{
			switch (lineCounter%5)
			{
			  case 0: left_sequenceID = line; break;
			  case 1: left_sequence = line; break;
			  case 2: break;
			  case 3:
			  {
				  left_score = line;

				  //Add to vector
				  ProbabilityMatrix left_probability_matrix(left_sequence,left_score);
				  leftProbabilityMatrixes.push_back(left_probability_matrix);
			  } break;
			  default: break;
			}
			lineCounter++;
		}

		left_file.close();
	}

	string right_sequenceID;
	string right_sequence;
	string right_score;

	ifstream right_file ("data/right.txt", ios::in);
	if (right_file.is_open())
	{
		lineCounter = 0;
		while (getline(right_file, line))
		{
			switch (lineCounter%5)
			{
			  case 0: right_sequenceID = line; break;
			  case 1: right_sequence = line; break;
			  case 2: break;
			  case 3:
			  {
				  right_score = line;

				  // Convert read from right file to reverse DNA inverse
				  char tempSequence;
				  char tempScore;
				  for(unsigned int i = 0; i<((right_sequence.size()/2)+(right_sequence.size()%2)); i++)
				  {
					  tempSequence = right_sequence[i];
					  right_sequence[i] = dnaInverse(right_sequence[right_sequence.size()-1-i]);
					  right_sequence[right_sequence.size()-1-i] = dnaInverse(tempSequence);

					  tempScore = right_score[i];
					  right_score[i] = right_score[right_score.size()-1-i];
					  right_score[right_score.size()-1-i] = tempScore;
				  }

				  //Add to vector
				  ProbabilityMatrix right_probability_matrix(right_sequence,right_score);
				  rightProbabilityMatrixes.push_back(right_probability_matrix);
			  } break;
			  default: break;
			}
			lineCounter++;
		}
		right_file.close();
	}

	for(unsigned int vi = 0; vi < leftProbabilityMatrixes.size() && vi < rightProbabilityMatrixes.size();vi++)
	{
		StringMatching in(&leftProbabilityMatrixes.at(vi), &rightProbabilityMatrixes.at(vi));

		ProbabilityMatrix resultingMatrix = in.joinedMatrix();
		resultingMatrix.setInitialMatrix();
		resultingMatrix.calculateDistributedMatrixProbability(100);

		n = resultingMatrix.sequence.size();
		y = new double * [n];
		for ( unsigned int i = 0; i < n; i++ )
			y[i] = new double [4];
		for(int mac=0;mac<n;mac++)
		{
			for(int coc=0;coc<4;coc++)
			{
				y[mac][coc]=resultingMatrix.getMatrix()[coc][mac];
			}
		}

		cout<<endl;
		cout<<resultingMatrix.sequence<<"\n";
		cout<<resultingMatrix.score<<"\n";
		resultingMatrix.printMatrix();

		start = clock();
		if ( mod == 0 )
		{
			string empty;
			if ( ! ( preparation ( empty, y, n, z, alphabet, mod ) ) )
			{
				return 0;
			}
			else
			{
				for ( unsigned int i = 0; i < n; i++ )
					delete[] y[i];
				delete[] y;
				unsigned int * WP = new unsigned int [n];
				wptable ( sigma, z, WP );

				unsigned int * borderArray = new unsigned int[n];
				borderArray = computerBorderArray(WP, n);

				finish = clock();
				double passtime = (	double ) ( finish - start ) / CLOCKS_PER_SEC;
				cout << "Elapsed time is " << passtime << endl;
	#if 1
				/*print*/
				cout << "\nWeighted Prefix Table:\n";
				for ( unsigned int i = 0; i < n; i++ )
				{
					cout << WP[i] << ' ';
				}

				cout << "\nWeighted Border Table:\n";
				for(int r = 0; r<n;r++)
				{
					cout<<borderArray[r]<<" ";
				}
				cout<<endl;
	#endif
			}
		}
	}

//	/* read input Weighted String */
//	ifstream weighted ( weightedStringFile );
//	if ( weighted.fail() )
//	{
//		cout << "Error: Cannot open the weighted string file!" << endl;
//		return 0;
//	}
//	else
//	{
//		double temp;
//		long int num = 0;
//		unsigned int column = alphabet.size();
//		while ( !weighted.eof() )
//		{
//			weighted >> temp;
//			num ++;
//		}
//		unsigned int row = num / column;
//		y = new double * [row];
//		for ( unsigned int i = 0; i < row; i++ )
//			y[i] = new double [column];
//		weighted.clear();
//		weighted.seekg( 0, ios::beg );
//		for ( unsigned int i = 0; i < row; i++ )
//		{
//			for ( unsigned int j = 0; j < column; j++ )
//			{
//				weighted >> y[i][j];
//			}
//		}
//		n = row;
//		weighted.close();
//	}
	return 0;
}
