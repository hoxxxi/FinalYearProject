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
	double z = 20; //The larger the threshold the higher the value of overlap TODO
	double ** y;			//weighted string
	unsigned int n;			//length of y
	clock_t start;
	clock_t finish;
	string line;
	int lineCounter;
	vector<Read> leftVector;
	vector<Read> rightVector;
	string left_sequence;
	string left_score;
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

	//Create output file
	ofstream file;
	file.open("output.txt");
	file.close();

	start = clock();
	for(unsigned int vi = 0; vi < leftVector.size() && vi < rightVector.size();vi++)
	{
		ProbabilityMatrix resultingMatrix(
				rightVector.at(vi).getSequence()+"N"+
				leftVector.at(vi).getSequence(),
				rightVector.at(vi).getScore()+"!"+
				leftVector.at(vi).getScore(),rightVector.at(vi).size());

		resultingMatrix.applyBigram(401); // window size 401
		resultingMatrix.applyQualityScore(100); // QS:BiGram = 100:1

		string empty;
		n = resultingMatrix.getSize();
		y = resultingMatrix.getMatrix();

		if ( ! ( preparation ( empty, y, n, z, alphabet, mod ) ) )
		{
			return 0;
		}
		//TODO ADD NORMAL PREFIX TABLE
		else
		{
			unsigned int * WPT = new unsigned int [n];
			wptable ( sigma, z, WPT );

			unsigned int * BT = (unsigned int *) calloc(n, sizeof(unsigned int));
			borderTable ( WPT, n, BT );

			/*print*/
			cout<<resultingMatrix.getSequence()<<endl;
			cout<<resultingMatrix.getScore()<<endl;
			cout << "Weighted Prefix Table:"<<endl;
			for ( unsigned int i = 0; i < n; i++ )
			{
				cout << WPT[i] << ' ';
			}
			cout << "\nWeighted Border Table:"<<endl;
			for(int r = 0; r<n;r++)
			{
				cout<<BT[r]<<" ";
			}
			resultingMatrix.printMatrix();

			int shortestReadLength = min(rightVector.at(vi).size(), leftVector.at(vi).size());
			int overlap = BT[n-1];
			if(shortestReadLength<overlap) {
				cout<<"Overlap specified by border larger then read size"<<endl;
				return 1;
//				overlap = borderArray[leftVector.at(vi).getSize()-1+temp-rightVector.at(vi).getSize()];

			}
			string joinedString = resultingMatrix.getSequence().substr(rightVector.at(vi).size()+1)+
					resultingMatrix.getSequence().substr(overlap,leftVector.at(vi).size()-overlap);
			cout<<joinedString<<endl<<endl;

			//Write to file
			ofstream file;
			file.open("output.txt", fstream::out | fstream::app);
			file<<joinedString+"\n";
			file.close();

			//Clean up
			delete[] WPT;
			delete[] BT;
		}
	}
	finish = clock();
	double passtime = (	double ) ( finish - start ) / CLOCKS_PER_SEC;
	cout << "\nElapsed time is " << passtime << endl;
	return 0;
}
