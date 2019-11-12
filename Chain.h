#ifndef CHAIN_H_
#define CHAIN_H_

#include <vector>
#include "Types.h"


class CHain {
public:
	vector<unsigned int> ch;
	int direction;
	GenomePos qStart, qEnd, tStart, tEnd;
	int main; // If one chain is aligned to several chromosomes, then we need to split this chain. main refers to the index of the main alignment in vector chain
	CHain () {
		qStart = 0;
		qEnd = 0;
		tStart = 0;
		tEnd = 0;
		direction = 0;
		main = -1;
	};
	~CHain() {};
	CHain (vector<unsigned int> &onechain) {
		ch = onechain;
		main = -1;
	}	
	CHain (GenomePos &qS, GenomePos &qE, GenomePos &tS, GenomePos &tE, vector<unsigned int> &onechain) {
		qStart = qS;
		qEnd = qE;
		tStart = tS;
		tEnd = tE;
		ch = onechain;
		main = -1;
	}
	CHain (GenomePos &qS, GenomePos &qE, GenomePos &tS, GenomePos &tE, vector<unsigned int> &onechain, int idx) {
		qStart = qS;
		qEnd = qE;
		tStart = tS;
		tEnd = tE;
		ch = onechain;
		main = idx;
	}
	int Overlaps (GenomePos &qS, GenomePos &qE, float rate);
};


int CHain::Overlaps (GenomePos &qS, GenomePos &qE, float rate) {

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

	if (max(ovp/denomA, ovp/denomB) >= rate) {return true;}
	//if ((float)ovp/(qEnd - qStart) >= rate) {return true;}
	else {return false;}
}


class Primary_chain {
public:
	vector<CHain> chains;
	Primary_chain () {};
	~Primary_chain () {};
	Primary_chain (CHain chain) {
		chains.push_back(chain);
	}
};


#endif