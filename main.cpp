#include <iostream>
#include <algorithm>
#include <fstream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>

#include "defs.h"
#include "global.h"
#include "Read.h"
#include "ProbabilityMatrix.h"
#include "BorderTable.cpp"

using namespace std;

int main (int argc, char **argv)
{
	cout<<"hello";
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

	vector<Read> leftVector;
	vector<Read> rightVector;

	string left_sequenceID;
	string left_sequence;
	string left_score;

	string right_sequenceID;
	string right_sequence;
	string right_score;

	//Read left file
	ifstream left_file ("left.txt", ios::in);
	if (left_file.is_open())
	{
		lineCounter = 0;
		while (getline(left_file, line))
		{
			switch (lineCounter%5)
			{
			  case 0: left_sequenceID = line; break;
			  case 1: left_sequence = line; break;
			  case 3:
			  {
				  left_score = line;

				  //Add to vector
				  Read left_read(left_sequence,left_score);
				  leftVector.push_back(left_read);
			  } break;
			  default: break;
			}
			lineCounter++;
		}
		left_file.close();
	}

	//Read right file
	ifstream right_file ("right.txt", ios::in);
	if (right_file.is_open())
	{
		lineCounter = 0;
		while (getline(right_file, line))
		{
			switch (lineCounter%5)
			{
			  case 0: right_sequenceID = line; break;
			  case 1: right_sequence = line; break;
			  case 3:
			  {
				  right_score = line;
				  //Add to vector
				  Read right_read(right_sequence,right_score);
				  right_read.calculateReadInverse();
				  rightVector.push_back(right_read);
			  } break;
			  default: break;
			}
			lineCounter++;
		}
		right_file.close();
	}

	for(unsigned int vi = 0; vi < leftVector.size() && vi < rightVector.size();vi++)
	{
		ProbabilityMatrix resultingMatrix(
				rightVector.at(vi).getSequence()+
				leftVector.at(vi).getSequence(),
				rightVector.at(vi).getScore()+
				leftVector.at(vi).getScore());

		resultingMatrix.applyBigram(401); // window size 401
		resultingMatrix.applyQualityScore(100); // QS:BiGram = 100:1

		n = resultingMatrix.getSize();
		y = new double * [n];
		for ( unsigned int i = 0; i < n; i++ )
			y[i] = new double [4];
		for(int mac=0;mac<n;mac++)
		{
			for(int coc=0;coc<4;coc++)
			{
				y[mac][coc]=resultingMatrix.getMatrix()[mac][coc];
			}
		}

		cout<<endl;
		cout<<resultingMatrix.getSequence()<<"\n";
		cout<<resultingMatrix.getScore()<<"\n";
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
