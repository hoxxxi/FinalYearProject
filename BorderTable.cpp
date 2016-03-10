/*
 * BorderTable.cpp
 *
 *  Created on: 19 Feb 2016
 *      Author: yordan
 */

#include "BorderTable.h"

//#include <fstream>
//#include <iostream>
//#include <sstream>
//#include <string>
//#include <vector>

using namespace std;

unsigned int* computerBorderArray(unsigned int* prefixArray, int n) {
	unsigned int* borderTable = new unsigned int[n];
	borderTable[0] = 0;
	int l = 0;

	for(int i = 1; i<=n-1; i++)
	{
		if(i+prefixArray[i]-1>=l)
		{
			for(int r = 0; r < i+prefixArray[i]-1-l; r++)
			{
				borderTable[i+prefixArray[i]-r-1]=prefixArray[i]-r;
			}
			l = i + prefixArray[i];
		}
	}
	return borderTable;
}
//int main()
//{
//	string line;
//	ifstream file ("output.txt", ios::in);
//	vector<int> prefixTable;
//	vector< vector<int> > prefixTableList;
//	if (file.is_open())
//	{
//		getline(file, line);
//		while (getline(file, line))
//		{
//			if (line.empty())
//			{
//				prefixTableList.push_back(prefixTable);
//				prefixTable.clear();
//				getline(file, line);
//			}
//			else
//			{
//				istringstream tmp(line);
//				int n;
//				tmp >> n;
////				prefixTable.push_back(n);
//			}
//		}
//		prefixTableList.push_back(prefixTable);
//	}
//	int test[] = {30, 5, 4, 11, 6, 5, 6, 10, 6, 6, 10, 7, 7, 8, 7, 8, 9, 6, 8, 6, 4, 6, 6, 5, 6, 5, 4, 3, 2, 1};
////	for(unsigned int i = 0;i< prefixTableList.size();i++)
////	{
//		int* borderArray = computerBorderArray(test/*&((prefixTableList.at(i))[0]) -- &prefixTableList[0][i]*/, 30/*prefixTableList[i].size()*/);
//
////		cout<<"\nWeighted Border Table: "<<i<<"\n";
//		for(int r = 0; r<30;r++)
//		{
//			cout<<borderArray[r]<<" ";
//		}
////		cout<<"Weighted Border Table: "<<i<<"\n";
////		for(unsigned int k = 0;i<prefixTableList[i].size();k++)
////		{
////			cout<<borderArray[k]<<" ";
////		}
////		cout<<endl;
////	}
//}
//
