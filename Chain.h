#ifndef CHAIN_H_
#define CHAIN_H_

#include <vector>
#include "Types.h"
#include "Clustering.h"


class CHain {
public:
	vector<unsigned int> ch;
	vector<bool> link;
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
	CHain (GenomePos &qS, GenomePos &qE, GenomePos &tS, GenomePos &tE, vector<unsigned int> &onechain, vector<bool> &lk) {
		qStart = qS;
		qEnd = qE;
		tStart = tS;
		tEnd = tE;
		ch = onechain;
		main = -1;
		link = lk;
		//value = val;
	}
	/*
	CHain (GenomePos &qS, GenomePos &qE, GenomePos &tS, GenomePos &tE, vector<unsigned int> &onechain, int idx) {
		qStart = qS;
		qEnd = qE;
		tStart = tS;
		tEnd = tE;
		ch = onechain;
		main = idx;
	}
	*/
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
	//float denomB = qE - qS;

	if (ovp/denomA >= rate) {return true;}
	//if (max(ovp/denomA, ovp/denomB) >= rate) {return true;}
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


class FinalChain {
public: 
	vector<Cluster> *ExtendClusters;
	vector<unsigned int> chain;
	vector<unsigned int> MatchStart;
	vector<int> ClusterIndex; // ClusterIndex[i] stores the index of the Cluster that anchor i comes from;
	vector<int> StartIndex; // StartIndex[i] stores the index of the start for Cluster inputchain[i];

	FinalChain (vector<Cluster> *clusters) {
		ExtendClusters = clusters;
	}

	void InitializeOtherParts (vector<unsigned int> & matchstart, int & totalMatch, vector<Fragment_Info> & Value) {
		MatchStart = matchstart;
		ClusterIndex.resize(totalMatch);
		StartIndex.resize(totalMatch);
		assert(ClusterIndex.size() == totalMatch);

		for (int v = 0; v < Value.size(); v ++) {
			ClusterIndex[v] = Value[v].clusterNum;
			StartIndex[v] = Value[v].matchstartNum;
		}
	}

	bool strand (int i) {
		return (*ExtendClusters)[ClusterIndex[chain[i]]].strand;
	}

	int ClusterNum (int i) {
		return ClusterIndex[chain[i]];
	}

	int size () {
		return chain.size();
	}

	void resize(int & m) {
		chain.resize(m);
	}

	GenomePair & operator[](int i) {
		int clusterNum = ClusterIndex[chain[i]];
		int matchstartNum = StartIndex[chain[i]];
		assert(chain[i] - MatchStart[matchstartNum] < (*ExtendClusters)[clusterNum].matches.size());
		return (*ExtendClusters)[clusterNum].matches[chain[i] - MatchStart[matchstartNum]];
	}

	int length (int i) {
		int clusterNum = ClusterIndex[chain[i]];
		int matchstartNum = StartIndex[chain[i]];		
		assert(chain[i] - MatchStart[matchstartNum] < (*ExtendClusters)[clusterNum].matchesLengths.size());
		return (*ExtendClusters)[clusterNum].matchesLengths[chain[i] - MatchStart[matchstartNum]];		
	}
};


class SplitChain {
public:
	vector<unsigned int> sptc;
	vector<bool> link;

	SplitChain (vector<unsigned int> &sp, vector<bool> &lk) {
		sptc = sp;
		link = lk;
	}

	int size() const {
		return sptc.size();
	}

	unsigned int & operator[] (int i) {
		return sptc[i];
	}
};
#endif