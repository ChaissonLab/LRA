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
	vector<bool> link;
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
	// void resize(int m) {
	// 	chain.resize(m);
	// 	ClusterIndex.resize(m);
	// }
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
	GenomePos QStart, QEnd, TStart, TEnd;

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

	GenomePos &qStart(int i) {return (*clusters)[ClusterIndex[i]].matches[chain[i]].first.pos;}
	GenomePos &tStart(int i) {return (*clusters)[ClusterIndex[i]].matches[chain[i]].second.pos;}
	GenomePos qEnd(int i) {return (*clusters)[ClusterIndex[i]].matches[chain[i]].first.pos + (*clusters)[ClusterIndex[i]].matchesLengths[chain[i]];}
	GenomePos tEnd(int i) {return (*clusters)[ClusterIndex[i]].matches[chain[i]].second.pos + (*clusters)[ClusterIndex[i]].matchesLengths[chain[i]];}
	// GenomePos trans_qStart(int i, int length) {
	// 	if (strand(i) == 0) {return qStart(i);}
	// 	else {return length - qEnd(i);} 
	// }	
	// GenomePos trans_qEnd(int i, int length) {
	// 	if (strand(i) == 0) {return qEnd(i);}
	// 	else {return length - qStart(i);}
	// }
	// GenomePos trans_tStart(int i, GenomePos offset) { return tStart(i) - offset;}
	// GenomePos trans_tEnd(int i, GenomePos offset) { return tEnd(i) - offset;}

	bool OverlapsOnT (GenomePos tS, GenomePos tE, float rate);
	void CleanSpuriousJumpingAnchors ();
	long diag(int i) {
		if (strand(i) == 1) { return (long) qEnd(i) + (long) tStart(i);}
		else { return (long) tStart(i) - (long) qStart(i);}
	}
	
	void Initialize(vector<unsigned int> &ch, vector<Fragment_Info> &Value, vector<bool> &lk) {
		chain.resize(ch.size());
		link = lk;
		ClusterIndex.resize(ch.size());
		for (int i = 0; i < ch.size(); i++) {
			ClusterIndex[i] = Value[ch[i]].clusterNum;
			assert((*clusters)[ClusterIndex[i]].matchStart != -1);
			chain[i] = ch[i] - (*clusters)[ClusterIndex[i]].matchStart;
		}
	}
};


bool UltimateChain::OverlapsOnT (GenomePos tS, GenomePos tE, float rate) {

	int ovp = 0;
	if (tS >= TStart and tS < TEnd) {
		ovp = min(tE, TEnd) - tS;
	}
	else if (tE > TStart and tE <= TEnd) {
		ovp = tE - max(tS, TStart);
	}
	else if (tS < TStart and tE > TEnd) {
		ovp = TEnd - TStart;
	}
	float denomA = TEnd - TStart;
	if (ovp/denomA <= rate) {return true;}
	else {return false;}
}

void UltimateChain::CleanSpuriousJumpingAnchors () {
	//
	// Clean spurious anchors jumping far
	//
	vector<bool> remove(chain.size(), 0);
	int im = 0, cur = 0, prev = 0;
	int jump = -1; GenomePos jump_tpos;
	while (im < chain.size() - 1) {
		cur = im + 1; prev = im;
		if (jump == -1) {
			if (strand(cur) == strand(prev)) {
				if (strand(cur) == 0) {
					if (tEnd(cur) > tStart(prev)) {jump = cur; jump_tpos = tStart(prev);}
				}
				else {
					if (tStart(cur) < tEnd(prev)) {jump = cur; jump_tpos = tEnd(prev);}
				}
			}
		}
		else {
			if (strand(cur) == 0) {
				if (tEnd(cur) <= jump_tpos and cur - jump <= 3) {
					for (int i = jump; i < cur; i++) {remove[i] = 1;}
					jump = -1;					
				}
			}
			else {
				if (tStart(cur) >= jump_tpos and cur - jump <= 3) {
					for (int i = jump; i < cur; i++) {remove[i] = 1;}
					jump = -1;					
				}
			}
		}
		im++;					
	}
	if (jump != -1 and cur - jump <= 3) {
		for (int i = jump; i <= cur; i++) {remove[i] = 1;}
	}
	int c = 0;
	for (int i = 0; i < remove.size(); i++) {
		if (!remove[i]) {
			chain[c] = chain[i]; 
			ClusterIndex[c] = ClusterIndex[i];
			if (c > 1) {link[c - 1] = link[i - 1];}
			c++;
		}
	}
	chain.resize(c);
	ClusterIndex.resize(c);
	link.resize(c-1);
	return;
}

class SplitChain {
public:
	vector<int> sptc;
	vector<bool> link;
	GenomePos QStart, QEnd, TStart, TEnd;
	int chromIndex;
	UltimateChain *chain;
	bool Strand;
	int clusterIndex; // for one refined_clusters
	vector<int> ClusterIndex; // one splitchain can correspond to multiple ext_clusters
	// int readlength;
	// int offset;
	SplitChain() {}
	~SplitChain() {}
	SplitChain (vector<int> &sp, vector<bool> &lk) {
		sptc = sp;
		link = lk;
	}
	SplitChain (vector<int> &sp, vector<bool> &lk, UltimateChain *c, bool str) : chain(c) {
		sptc = sp;
		link = lk;
		Strand = str;
	}
	// SplitChain (int c) {
	// 	clusterIndex = c;
	// }
	int size() const {
		return sptc.size();
	}

	int & operator[] (int i) {
		return sptc[i];
	}

	GenomePair & genomepair(int i) {
		return (*chain)[sptc[i]];
	}

	int length(int i) {
		return (*chain).length(sptc[i]);
	}

	int CHROMIndex(Genome & genome) {
		if (sptc.size() == 0) return 1;
		int firstChromIndex = genome.header.Find(TStart);
		int lastChromIndex;
		lastChromIndex = genome.header.Find(TEnd);
		if (firstChromIndex != lastChromIndex ) {return 1;}
		chromIndex = firstChromIndex;  
		return 0;
	}	
	bool strand (int i) { return chain->strand(sptc[i]);}
	GenomePos &qStart(int i) {return chain->qStart(sptc[i]);}
	GenomePos &tStart(int i) {return chain->tStart(sptc[i]);}
	GenomePos qEnd(int i) {return chain->qEnd(sptc[i]);}
	GenomePos tEnd(int i) {return chain->tEnd(sptc[i]);}

	// GenomePos trans_qStart(int i) {
	// 	if (strand(i) == 0) {return qStart(i);}
	// 	else {return readlength - qEnd(i);} 
	// }	
	// GenomePos trans_qEnd(int i) {
	// 	if (strand(i) == 0) {return qEnd(i);}
	// 	else {return readlength - qStart(i);}
	// }
	// GenomePos trans_tStart(int i) { return tStart(i) - offset;}
	// GenomePos trans_tEnd(int i) { return tEnd(i) - offset;}

	//
	// Sort first by tStart, then by qStart
	//
	int operator()(const int &a, const int &b) {
		if (tStart(a) != tStart(b)) { return tStart(a) < tStart(b);}
		else { return qStart(a) < qStart(b);}
	}

	int CartesianTargetLowerBound(int b, int e, int query) {
		return lower_bound(sptc.begin() - b, sptc.end() - e, query, *this) - sptc.begin();
	}

	int CartesianTargetUpperBound(int b, int e, int query) {
		return upper_bound(sptc.begin() - b, sptc.end() - e, query, *this) - sptc.begin();
	}
};

class Merge_SplitChain {
public:
	vector<Cluster *> *clusters;
	vector<int> merged_clusterIndex;
	GenomePos qStart, qEnd, tStart, tEnd;
	
	Merge_SplitChain (vector<int> &m, vector<Cluster *> * c) : merged_clusterIndex(m), clusters(c) {}
	void Update_boundary() {
		qStart = (*clusters)[merged_clusterIndex[0]]->qStart;
		tStart = (*clusters)[merged_clusterIndex[0]]->tStart;
		qEnd = (*clusters)[merged_clusterIndex[0]]->qEnd;
		tEnd = (*clusters)[merged_clusterIndex[0]]->tEnd;
		for (int t = 0; t < merged_clusterIndex.size(); t++) {
		qStart = min((*clusters)[merged_clusterIndex[0]]->qStart, qStart);
		tStart = min((*clusters)[merged_clusterIndex[0]]->tStart, tStart);
		qEnd = max((*clusters)[merged_clusterIndex[0]]->qEnd, qEnd);
		tEnd = max((*clusters)[merged_clusterIndex[0]]->tEnd, tEnd);			
		}
	} 
};


#endif