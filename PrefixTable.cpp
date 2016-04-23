/*
 *  Created on: 1 Feb 2016
 *  Author: Yordan Petrov Yordanov
 */

#include "PrefixTable.h"

using namespace std;

void calculatePrefixTable(string x, int m, unsigned int* pref) {
	pref[0]=m;
	int g = 0;
	int f = 0;
	for(int i = 1; i<m; i ++) {
		if(i<g && pref[i-f]!=g-i) {
			pref[i] = min((int) pref[i-f], g-i);
		}
		else {
			g= max(g, i);
			f= i;
			while(g<m && x[g] == x[g-f]) {
				g++;
			}
			pref[i]=g-i;
		}
	}
}
