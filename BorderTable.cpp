/*
 *  Created on: 19 Feb 2016
 *  Author: Yordan Petrov Yordanov
 */

#include "BorderTable.h"

using namespace std;

unsigned int borderTable(unsigned int* prefixTablePtr, int n, unsigned int * borderTablePtr) {
	borderTablePtr[0] = 0;
	int l = 0;
	for(int i = 1; i<=n-1; i++) {
		if(i+prefixTablePtr[i]-1>=l) {
			for(int r = 0; r < i+prefixTablePtr[i]-1-l; r++) {
				borderTablePtr[i+prefixTablePtr[i]-r-1]=prefixTablePtr[i]-r;
			}
			l = i + prefixTablePtr[i];
		}
	}
	return 1;
}
