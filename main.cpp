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
#include "BorderTable.h"

using namespace std;

int main (int argc, char **argv)
{
	string alphabet = DNA;
	unsigned int sigma = alphabet.size();
	int mod = 0;
	double z = 10; //The larger the threshold the higher the value of overlap
	double z1 = 10; //The lower the threshold the higher the tolerance to overlap
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
			switch (lineCounter%4)
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
			switch (lineCounter%4)
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

	stringstream output;
	for(unsigned int vi = 0; vi < leftVector.size() && vi < rightVector.size();vi++)
	{
		ProbabilityMatrix resultingMatrix(
				rightVector.at(vi).getSequence()+
				leftVector.at(vi).getSequence(),
				rightVector.at(vi).getScore()+
				leftVector.at(vi).getScore());

		start = clock();

		resultingMatrix.applyBigram(401); // window size 401
		resultingMatrix.applyQualityScore(1); // QS:BiGram = 100:1

		cout<<resultingMatrix.getSequence()<<"\n";
		cout<<resultingMatrix.getScore()<<"\n";
		resultingMatrix.printMatrix();

		n = resultingMatrix.getSize();
		y = resultingMatrix.getMatrix();

		string empty;
		unsigned int * borderArray;
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

			borderArray = new unsigned int[n];
			borderArray = computerBorderArray(WP, n);

			finish = clock();
			double passtime = (	double ) ( finish - start ) / CLOCKS_PER_SEC;
			cout << "Elapsed time is " << passtime << endl;

			/*print*/
			cout << "Weighted Prefix Table:\n";
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
		}

		int shortestReadLength = min(rightVector.at(vi).getSize(), leftVector.at(vi).getSize());
		int overlap = borderArray[n-1];
		if(shortestReadLength<overlap)
		{
			int window;
			for(int i = 0; i < shortestReadLength; i++)
			{
				window = leftVector.at(vi).getSize()-i;
				if(leftVector.at(vi).getSequence().substr(leftVector.at(vi).getSize()-window,window)==rightVector.at(vi).getSequence().substr(0,window))
					break;
			}
			overlap = 0;
			for(int k = 0; k<window;k++)
			{
				double averageOfScoreValues = (1-pow(10.0,(33.0-((double)((int)(leftVector.at(vi).getScore()[leftVector.at(vi).getSize()-1-k]))))/10.0) +
						1-pow(10.0,(33.0-((double)((int)(leftVector.at(vi).getScore()[k]))))/10.0)) / 2;
				if(averageOfScoreValues > 1.0-1/z1)
					overlap++;
			}
		}

		string joinedString = resultingMatrix.getSequence().substr(rightVector.at(vi).getSize())+
				resultingMatrix.getSequence().substr(overlap,rightVector.at(vi).getSize());
		cout<< joinedString<<endl;

		output<<joinedString+"\n";
		cout<<endl;
	}
	ofstream file;
	file.open("output.txt");
	file<<output.str();
	file.close();
	return 0;
}
