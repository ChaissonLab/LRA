#ifndef CHAIN_H_
#define CHAIN_H_

#include <vector>
#include "Types.h"

typedef vector<unsigned int> chain;


class Primary_chain {
public:
	GenomePos qStart, qEnd, tStart, tEnd;
	vector<chain> chains;
	vector<bool> direction; 
	Primary_chain () {
		qStart = 0;
		qEnd = 0;
		tStart = 0;
		tEnd = 0;
	};
	~Primary_chain () {};
	Primary_chain (GenomePos &qS, GenomePos &qE, GenomePos &tS, GenomePos &tE, vector<unsigned int> &onechain) {
		qStart = qS;
		qEnd = qE;
		tStart = tS;
		tEnd = tE;
		chains.push_back(onechain);
	}
	int Overlaps (GenomePos &qS, GenomePos &qE, float rate);
};

int Primary_chain::Overlaps (GenomePos &qS, GenomePos &qE, float rate) {

	int ovp = 0;
	if (qS >= qStart and qS < qEnd) {
		ovp = min(qE, qEnd) - qS;
	}
	else if (qE > qStart and qE <= qEnd) {
		ovp = qE - max(qS, qStart);
	}
	else if (qS < qStart and qE > qEnd) {
		ovp = qEnd - qStart;
	}
	float denomA = qEnd - qStart;
	float denomB = qE - qS;
/*
	if (qS <= qStart and qE >= qEnd) {
		ovp = qEnd - qS;
	}
	else  if (qS >= qStart and qE <= qEnd) {
		ovp = qE - qStart;
	}
*/
	if (max(ovp/denomA, ovp/denomB) >= rate) {return true;}
	//if ((float)ovp/(qEnd - qStart) >= rate) {return true;}
	else {return false;}
}

#endif