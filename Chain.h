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
	int main; // If one chain is aligned to several chromosomes, then we need to split this chain. main refers to the index of the main alignment in vector chain.
	int chromIndex;
	float value;
	int NumOfAnchors0;
	CHain () {
		qStart = 0;
		qEnd = 0;
		tStart = 0;
		tEnd = 0;
		direction = 0;
		main = -1;
		value = 0;
		NumOfAnchors0 = 0;
	};
	~CHain() {};
	CHain (vector<unsigned int> &onechain) {
		ch = onechain;
		main = -1;
	}	
	CHain (GenomePos &qS, GenomePos &qE, GenomePos &tS, GenomePos &tE, vector<unsigned int> &onechain, 
						vector<bool> &lk, float &val, int &numofanchors) {
		qStart = qS;
		qEnd = qE;
		tStart = tS;
		tEnd = tE;
		ch = onechain;
		main = -1;
		link = lk;
		value = val;
		NumOfAnchors0 = numofanchors;
	}
	int OverlapsOnQ (GenomePos &qS, GenomePos &qE, float rate);
	int OverlapsOnT (GenomePos &tS, GenomePos &tE, float rate);
};


int CHain::OverlapsOnQ (GenomePos &qS, GenomePos &qE, float rate) {

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


int CHain::OverlapsOnT (GenomePos &tS, GenomePos &tE, float rate) {

	int ovp = 0;
	if (tS >= tStart and tS < tEnd) {
		ovp = min(tE, tEnd) - tS;
	}
	else if (tE > tStart and tE <= tEnd) {
		ovp = tE - max(tS, tStart);
	}
	else if (tS < tStart and tE > tEnd) {
		ovp = tEnd - tStart;
	}
	float denomA = tEnd - tStart;
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
	vector<Cluster_SameDiag *> * ExtendClusters;
	vector<unsigned int> chain;
	vector<int> ClusterIndex; // ClusterIndex[i] stores the index of the Cluster that anchor i comes from;
	float SecondSDPValue;
	// vector<unsigned int> MatchStart;
	// vector<int> StartIndex; // StartIndex[i] stores the index of the start for Cluster inputchain[i];
	FinalChain() {}
	FinalChain (vector<Cluster_SameDiag *> * c) : ExtendClusters(c) {}
	~FinalChain() {};
	void Initialize(vector<unsigned int> &ch, vector<Fragment_Info> &Value) {
		chain.resize(ch.size());
		ClusterIndex.resize(ch.size());
		for (int i = 0; i < ch.size(); i++) {
			ClusterIndex[i] = Value[ch[i]].clusterNum;
			assert((*ExtendClusters)[ClusterIndex[i]]->matchStart != -1);
			chain[i] = ch[i] - (*ExtendClusters)[ClusterIndex[i]]->matchStart;
		}
	}
	int size () {
		return chain.size();
	}
	void clear() {
		chain.clear();
		ClusterIndex.clear();

	}
	int ClusterNum (int i) {
		return ClusterIndex[i];
	}
	bool strand (int i) {
		assert(ClusterIndex[i] < ExtendClusters->size());
		return (*ExtendClusters)[ClusterIndex[i]]->strand;
	}
	void resize(int m) {
		chain.resize(m);
		ClusterIndex.resize(m);
	}
	GenomePos qStart(int i) {
		int clusterNum = ClusterIndex[i];
		assert(chain[i] < (*ExtendClusters)[clusterNum]->size());
		return 	(*ExtendClusters)[clusterNum]->GetqStart(chain[i]);
	}
	GenomePos tStart(int i) {
		int clusterNum = ClusterIndex[i];
		assert(chain[i] < (*ExtendClusters)[clusterNum]->size());
		return 	(*ExtendClusters)[clusterNum]->GettStart(chain[i]);
	}
	int length (int i) {
		int clusterNum = ClusterIndex[i];
		assert(chain[i] < (*ExtendClusters)[clusterNum]->size());
		return (*ExtendClusters)[clusterNum]->length(chain[i]);		
	}
};

class UltimateChain {
public: 
	vector<Cluster> * clusters;
	vector<unsigned int> chain;
	vector<bool> link;
	vector<int> ClusterIndex; // ClusterIndex[i] stores the index of the Cluster that anchor i comes from;
	float SecondSDPValue;
	float FirstSDPValue;
	int NumOfAnchors0;
	UltimateChain() {}
	UltimateChain(vector<Cluster> * c, float &s) : clusters(c), SecondSDPValue(s) {}
	UltimateChain(vector<Cluster> * cl): clusters(cl) {}
	~UltimateChain() {};
	bool strand (int i) {
		assert(ClusterIndex[i] < clusters->size());
		return (*clusters)[ClusterIndex[i]].strand;
	}

	int ClusterNum (int i) {
		return ClusterIndex[i];
	}

	int size () {
		return chain.size();
	}

	void resize(int m) {
		chain.resize(m);
		ClusterIndex.resize(m);
	}

	void clear() {
		chain.clear();
		ClusterIndex.clear();
	}

	GenomePair & operator[](int i) {
		int clusterNum = ClusterIndex[i];
		assert(chain[i] < (*clusters)[clusterNum].matches.size());
		return (*clusters)[clusterNum].matches[chain[i]];
	}

	int length (int i) {
		int clusterNum = ClusterIndex[i];
		assert(chain[i] < (*clusters)[clusterNum].matches.size());
		return (*clusters)[clusterNum].matchesLengths[chain[i]];		
	}

	GenomePos & qStart(int &i) {
		return (*clusters)[ClusterIndex[i]].matches[chain[i]].first.pos;
	}
	GenomePos & tStart(int &i) {
		return (*clusters)[ClusterIndex[i]].matches[chain[i]].second.pos;
	}
	GenomePos qEnd(int &i) {
		return (*clusters)[ClusterIndex[i]].matches[chain[i]].first.pos + (*clusters)[ClusterIndex[i]].matchesLengths[chain[i]];
	}
	GenomePos tEnd(int &i) {
		return (*clusters)[ClusterIndex[i]].matches[chain[i]].second.pos + (*clusters)[ClusterIndex[i]].matchesLengths[chain[i]];
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