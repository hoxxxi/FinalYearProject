/*
 *  Created on: 10 Mar 2016
 *  Author: Yordan Petrov Yordanov
 */

#ifndef READ_H_
#define READ_H_

#include <string>
#include <iostream>

using namespace std;

class Read {
	string sequence;
	string score;
public:
	Read(string sequenceIn, string scoreIn);
	virtual ~Read();
	string getSequence();
	string getScore();
	int size();
	void calculateReadInverse();
};

char dnaInverse(char in);

#endif /* READ_H_ */
