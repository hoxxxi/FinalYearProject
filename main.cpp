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
#include "PrefixTable.h"
#include "BorderTable.h"

using namespace std;

int main (int argc, char **argv)
{
	string alphabet = DNA;
	unsigned int sigma = alphabet.size();
	int mod = 0;
	double z = 2;
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
	ifstream left_file ("left.fastq", ios::in); // /home/yordan/Desktop/
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
	ifstream right_file ("right.fastq", ios::in); // /home/yordan/Desktop/
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
				  cout<<"Right Case: "<<lineCounter/4<<endl;
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
				rightVector.at(vi).getSequence()+//"N"+
				leftVector.at(vi).getSequence(),
				rightVector.at(vi).getScore()+//"!"+
				leftVector.at(vi).getScore(),
				rightVector.at(vi).size());

		resultingMatrix.applyBigram(401); // window size 401
		resultingMatrix.applyQualityScore(100); // QS:BiGram = 100:1

		string empty;
		unsigned int * WPT = new unsigned int [resultingMatrix.getSize()];

		unsigned int * BT = (unsigned int *) calloc(resultingMatrix.getSize(), sizeof(unsigned int));

		if ( ! ( preparation ( empty, resultingMatrix.getMatrix(), resultingMatrix.getSize(), z, alphabet, mod ) ) )
		{
			calculatePrefixTable(resultingMatrix.getSequence(),resultingMatrix.getSize(),WPT);
		}
		else
		{
			wptable ( sigma, z, WPT );
		}

		borderTable ( WPT, resultingMatrix.getSize(), BT );

		int shortestReadLength = min(rightVector.at(vi).size(), leftVector.at(vi).size());
		int overlap = BT[resultingMatrix.getSize()-1];
		while(shortestReadLength<overlap) {
			overlap = BT[overlap-1];
		}
		string joinedString = resultingMatrix.getSequence().substr(rightVector.at(vi).size())+
				resultingMatrix.getSequence().substr(overlap,leftVector.at(vi).size()-overlap);

		/*print*/
		cout<<resultingMatrix.getSequence()<<endl;
		cout<<resultingMatrix.getScore()<<endl;
		cout << "Weighted Prefix Table:"<<endl;
		for ( unsigned int i = 0; i < resultingMatrix.getSize(); i++ )
		{
			cout << WPT[i] << ' ';
		}
		cout << "\nWeighted Border Table:"<<endl;
		for(int r = 0; r<resultingMatrix.getSize();r++)
		{
			cout<<BT[r]<<" ";
		}
		resultingMatrix.printMatrix();

		cout<<joinedString<<endl;

		//Write to file
		ofstream file;
		file.open("output.txt", fstream::out | fstream::app);
		file<<joinedString+"\n";
		file.close();

		//Clean up
		delete[] WPT;
		delete[] BT;
		leftVector.clear();
		rightVector.clear();

	}
	finish = clock();
	double passtime = (	double ) ( finish - start ) / CLOCKS_PER_SEC;
	cout << "\nElapsed time is " << passtime << endl;
	return 0;
}
