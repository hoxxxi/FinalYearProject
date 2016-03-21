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
	for(int zValue = 10; zValue< 100; zValue++){
		cout<<"\nTreshhold: "<<zValue<<endl;
	string alphabet = DNA;
	unsigned int sigma = alphabet.size();
	int mod = 0;
	double z = zValue;
	clock_t start;
	clock_t finish;

	ifstream left_file ("left.fastq", ios::in); // /home/yordan/Desktop/
	ifstream right_file ("right.fastq", ios::in);

	ofstream file;
	file.open("output.txt");
	file.close();
	
	if (left_file.is_open() && right_file.is_open())
	{
		string leftLine;
		string rightLine;
		int lineCounter = 0;
		string left_sequence;
		string left_score;
		string right_sequence;
		string right_score;
		start = clock();
		while (getline(left_file, leftLine) && getline(right_file, rightLine))
		{
			switch (lineCounter%4)
			{
			case 1: left_sequence = leftLine; right_sequence = rightLine; break;
			case 3:
				{

					left_score = leftLine;
					Read left_read(left_sequence,left_score);

					right_score = rightLine;
					Read right_read(right_sequence,right_score);
					right_read.calculateReadInverse();

					ProbabilityMatrix resultingMatrix(
							right_read.getSequence()+//"N"+
							left_read.getSequence(),
							right_read.getScore()+//"!"+
							left_read.getScore(),
							right_read.size());

					resultingMatrix.applyBigram(401); // window size 401
					resultingMatrix.applyQualityScore(100); // QS:BiGram = 100:1

					string empty;
					unsigned int * PT = new unsigned int [resultingMatrix.getSize()];
					unsigned int * BT = (unsigned int *) calloc(resultingMatrix.getSize(), sizeof(unsigned int));


					preparation ( empty, resultingMatrix.getMatrix(), resultingMatrix.getSize(), z, alphabet, mod );

					/*print*/
					cout<<"\n\nRead No. "<< (lineCounter/4)+1 <<endl;
//					cout<<resultingMatrix.getSequence().substr(right_sequence.size())+"\t"+resultingMatrix.getSequence().substr(0,left_sequence.size())<<endl;

					wptable ( sigma, z, PT ); // Weighted Prefix Table
					borderTable ( PT, resultingMatrix.getSize(), BT );

//					cout << "Weighted Prefix Table:"<<endl;
//					for ( unsigned int i = 0; i < resultingMatrix.getSize(); i++ )
//					{
//						cout << PT[i] << ' ';
//					}
					cout << "Weighted Border Table:"<<endl;
					for(int r = 0; r<resultingMatrix.getSize();r++)
					{
						cout<<BT[r]<<" ";
					}

					//Clean up
					delete[] PT;
					free (BT);
					PT = new unsigned int [resultingMatrix.getSize()];
					BT = (unsigned int *) calloc(resultingMatrix.getSize(), sizeof(unsigned int));
					calculatePrefixTable(resultingMatrix.getSequence(),resultingMatrix.getSize(), PT); // Normal Prefix Table
					borderTable ( PT, resultingMatrix.getSize(), BT );

//					cout << "\nNormal Prefix Table:"<<endl;
//					for ( unsigned int i = 0; i < resultingMatrix.getSize(); i++ )
//					{
//						cout << PT[i] << ' ';
//					}
//					cout << "\nNormal Border Table:"<<endl;
//					for(int r = 0; r<resultingMatrix.getSize();r++)
//					{
//						cout<<BT[r]<<" ";
//					}

//					int shortestReadLength = min(right_read.size(), left_read.size());
//					int overlap = BT[resultingMatrix.getSize()-1];
//
//					while(shortestReadLength<overlap) {
//						overlap = BT[overlap-1]; // Get the border of the border
//					}
//
//					string joinedString = resultingMatrix.getSequence().substr(right_read.size())+
//							resultingMatrix.getSequence().substr(overlap,left_read.size()-overlap);
//
//					//Write to file
//					ofstream file;
//					file.open("output.txt", fstream::out | fstream::app);
//					file<<joinedString+"\n";
//					file.close();

					//Clean up
					delete[] PT;
					free (BT);
				} break;
			default: break;
			}
			lineCounter++;
		}
		left_file.close();
		right_file.close();
	}

	finish = clock();
	double passtime = (	double ) ( finish - start ) / CLOCKS_PER_SEC;
	cout << "\n\nElapsed time is " << passtime << endl;}
	return 0;
}
