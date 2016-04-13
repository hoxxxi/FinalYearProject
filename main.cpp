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
	TSwitch sw;
	string alphabet = DNA;
	unsigned int sigma = alphabet.size();
	int mod = 0;
	string left_file = "./data/left.fastq";
	string right_file = "./data/right.fastq";
	string output = "output.txt";
	double x = 10;
	double z = 10;
	int bigramWindow = 400;
	int qualityScoreCoeficient = 100;
	clock_t start;
	clock_t finish;

#if 1
	/* Decodes the arguments */
	unsigned int k = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
	if ( k < 7 || k > 15)
	{
		usage();
		return 1;
	}
	else
	{
		if ( sw.left_filename.size() == 0 )
		{
			cout << "Error: No left FASTQ input!" << endl;
			return 0;
		}
		else
		{
			left_file = sw.left_filename;
		}

		if ( sw.right_filename.size() == 0)
		{
			cout << "Error: No right FASTQ input!" << endl;
			return 0;
		}
		else
		{
			right_file = sw.right_filename;
		}

		if ( sw.output_filename.size() == 0 )
		{
			cout << "Error: Specify output file!" << endl;
		}
		else
		{
			output = sw.output_filename;
		}

		z = sw.z;
		x = sw.x;
		bigramWindow = sw.bigramWindow;
		qualityScoreCoeficient = sw.qualityScoreCoefficient;
	}
#endif

	ifstream leftFileStream (left_file, ios::in); // /home/yordan/Desktop/
	ifstream rightFileStream (right_file, ios::in);

	ofstream file;
	file.open(output.c_str());
	file.close();
	
	if (leftFileStream.is_open() && rightFileStream.is_open())
	{
		string leftLine;
		string rightLine;
		int lineCounter = 0;
		string left_sequence;
		string left_score;
		string right_sequence;
		string right_score;
		string lineLabel;
		start = clock();
		while (getline(leftFileStream, leftLine) && getline(rightFileStream, rightLine))
		{
			switch (lineCounter%4)
			{
			case 0: lineLabel = leftLine.substr(1,leftLine.size()-3); break;
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

					resultingMatrix.applyBigram(bigramWindow); // window size 401
					resultingMatrix.applyQualityScore(100); // QS:BiGram = 100:1

					string empty;
					unsigned int * PT = new unsigned int [resultingMatrix.getSize()];
					unsigned int * BT = (unsigned int *) calloc(resultingMatrix.getSize(), sizeof(unsigned int));

					if ( preparation ( empty, resultingMatrix.getMatrix(), resultingMatrix.getSize(), z, alphabet, mod ) ==  0 )
					{
						wptable ( sigma, z, PT ); // Weighted Prefix Table
					}
					else
					{
						calculatePrefixTable(resultingMatrix.getSequence(),resultingMatrix.getSize(), PT); // Normal Prefix Table
					}
					borderTable ( PT, resultingMatrix.getSize(), BT );

					/*print*/
#if 0
					cout<<"Read No. "<< (lineCounter/4)+1 <<" with label: "<<lineLabel<<endl;
					cout<<resultingMatrix.getSequence()<<endl;
					cout<<resultingMatrix.getScore()<<endl;
					cout << "Weighted Prefix Table:"<<endl;
					for ( unsigned int i = 0; i < resultingMatrix.getSize(); i++ )
					{
						cout << PT[i] << ' ';
					}
					cout << "\nWeighted Border Table:"<<endl;
					for(int r = 0; r<resultingMatrix.getSize();r++)
					{
						cout<<BT[r]<<" ";
					}
					resultingMatrix.printMatrix();
#endif

					int shortestReadLength = min(right_read.size(), left_read.size());
					int overlap = BT[resultingMatrix.getSize()-1];

					while(shortestReadLength<overlap) {
						overlap = BT[overlap-1]; // Get the border of the border
					}

					string joinedString = resultingMatrix.getSequence().substr(right_read.size())+
							resultingMatrix.getSequence().substr(overlap,left_read.size()-overlap);

					//Write to file
					if(((double) overlap)/((double)shortestReadLength)>=((double)x/100.0)){
						ofstream file;
						file.open(output.c_str(), fstream::out | fstream::app);
						file<<">"+lineLabel+"\n"+joinedString+"\n";
						file.close();
//						cout<<overlap<<"\t"<<shortestReadLength<<"\t"<<(double) overlap/shortestReadLength<<endl;
					}
					//Clean up
					delete[] PT;
					free (BT);
				} break;
			default: break;
			}
			lineCounter++;
		}
		leftFileStream.close();
		rightFileStream.close();
	}

	finish = clock();
	double passtime = (	double ) ( finish - start ) / CLOCKS_PER_SEC;
	cout << "\nElapsed time is " << passtime << endl;
	return 0;
}
