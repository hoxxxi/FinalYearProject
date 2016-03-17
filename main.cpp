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
	clock_t start;
	clock_t finish;

	ifstream left_file ("/home/yordan/Desktop/frag_1.fastq", ios::in);
	ifstream right_file ("/home/yordan/Desktop/frag_2.fastq", ios::in);

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
							right_read.getSequence()+"N"+
							left_read.getSequence(),
							right_read.getScore()+"!"+
							left_read.getScore(),right_read.size());

					resultingMatrix.applyBigram(401); // window size 401
					resultingMatrix.applyQualityScore(100); // QS:BiGram = 100:1

					string empty;
					if ( ! ( preparation ( empty, resultingMatrix.getMatrix(), resultingMatrix.getSize(), z, alphabet, mod ) ) )
					{
						return 0;
					}
//					//TODO ADD NORMAL PREFIX TABLE
//					else
//					{
//						unsigned int * WPT = new unsigned int [resultingMatrix.getSize()];
//						wptable ( sigma, z, WPT );
//
//						unsigned int * BT = (unsigned int *) calloc(resultingMatrix.getSize(), sizeof(unsigned int));
//						borderTable ( WPT, resultingMatrix.getSize(), BT );
//
//						/*print*/
//						cout<<"Read No. "<< lineCounter/4<<endl;
//						cout<<resultingMatrix.getSequence()<<endl;
//						cout<<resultingMatrix.getScore()<<endl;
//						cout << "Weighted Prefix Table:"<<endl;
//						for ( unsigned int i = 0; i < resultingMatrix.getSize(); i++ )
//						{
//							cout << WPT[i] << ' ';
//						}
//						cout << "\nWeighted Border Table:"<<endl;
//						for(int r = 0; r<resultingMatrix.getSize();r++)
//						{
//							cout<<BT[r]<<" ";
//						}
//						resultingMatrix.printMatrix();
//
//						int shortestReadLength = min(right_read.size(), left_read.size());
//						int overlap = BT[resultingMatrix.getSize()-1];
//
//						//Clean up
//						delete[] WPT;
//						delete[] BT;
//
//						if(shortestReadLength<overlap) {
//							cout<<"Overlap specified by border larger then read size"<<endl;
//							return 1;
//			//				overlap = borderArray[leftVector.at(vi).getSize()-1+temp-rightVector.at(vi).getSize()];
//
//						}
					int overlap = 0;
					string joinedString = resultingMatrix.getSequence().substr(right_read.size()+1)+
							resultingMatrix.getSequence().substr(overlap,left_read.size()-overlap);

					//Write to file
					ofstream file;
					file.open("output.txt", fstream::out | fstream::app);
					file<<joinedString+"\n";
					file.close();
//					}
//					cout<<"Left Line "<<leftLine<<endl;
//					cout<<"Right Line "<<rightLine<<endl;
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
	cout << "\nElapsed time is " << passtime << endl;
	return 0;
}
