/*
 * Read.cpp
 *
 *  Created on: 10 Mar 2016
 *      Author: yordan
 */

#include "Read.h"

Read::Read(string sequenceIn, string scoreIn) {
	sequence = sequenceIn;
	score = scoreIn;
}

Read::~Read() {
//	cout<<"Read destroyed: "<<sequence<<endl;
}

string Read::getSequence() {
	return sequence;
}

string Read::getScore() {
	return score;
}

int Read::size() {
	return sequence.size();
}
void Read::calculateReadInverse() {
	// Convert read from right file to DNA inverse
	char tempSequence;
	char tempScore;
	for(unsigned int i = 0; i<((sequence.size()/2)+(sequence.size()%2)); i++)
	{
	  tempSequence = sequence[i];
	  sequence[i] = dnaInverse(sequence[sequence.size()-1-i]);
	  sequence[sequence.size()-1-i] = dnaInverse(tempSequence);

	  tempScore = score[i];
	  score[i] = score[score.size()-1-i];
	  score[score.size()-1-i] = tempScore;
	}
}

char dnaInverse(char in) {
	switch(in)
	{
	case 'A': return 'T';
	case 'C': return 'G';
	case 'G': return 'C';
	case 'T': return 'A';
	default: return 'N';
	}
}
