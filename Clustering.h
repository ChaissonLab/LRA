#ifndef CLUSTERING_H_
#define CLUSTERING_H_
#include "Options.h"
#include <vector>
#include <climits>
#include <list>
using namespace std;

template<typename Tup>
long DiagonalDifference(Tup &a, Tup &b, int strand=0) {
	if (strand == 0) { // forMatches
		long aDiag = (long)a.second.pos - (long)a.first.pos;
		long bDiag = (long)b.second.pos - (long)b.first.pos;
		return aDiag - bDiag;		
	}
	else { // revMatches
		long aDiag = a.first.pos + a.second.pos; 
		long bDiag= b.first.pos + b.second.pos;
		return aDiag - bDiag;				
	}
}

template<typename Tup>
long SkewDiagonalDifference(Tup &a, Tup &b, Options &opts, int strand=0) {
	if (strand == 0) { // forMatches
		long aDiag = (long)a.second.pos - (long)ceil(opts.slope*a.first.pos);
		long bDiag = (long)b.second.pos - (long)ceil(opts.slope*b.first.pos);
		return aDiag - bDiag;		
	}
	else { // revMatches
		long aDiag = ceil(opts.slope*a.first.pos) + a.second.pos; 
		long bDiag= ceil(opts.slope*b.first.pos) + b.second.pos;
		return aDiag - bDiag;				
	}
}

template<typename Tup>
long GapDifference(Tup &a, Tup &b) {
	long aDiff = abs((long)b.first.pos - (long)a.first.pos);
	long bDiff = abs((long)b.second.pos - (long)a.second.pos);
	return max(aDiff, bDiff);
}

template<typename Tup>
int DiagonalDrift(long curDiag, Tup &t, int strand=0) {
	long drift;
	if (strand == 0) drift= abs(curDiag - ((long)t.second.pos - (long)t.first.pos));
	else drift= abs(curDiag - ((long)t.first.pos + (long)t.second.pos));
	return drift;
}


template<typename Tup>
void CleanOffDiagonal(vector<pair<Tup, Tup> > &matches, Options &opts, int strand=0, int diagOrigin=-1, int diagDrift=-1) {
	if (matches.size() == 0) {
		return;
	}
	
	vector<bool> onDiag(matches.size(), false);
	int nOnDiag=0;
	if (matches.size() > 1 and abs(DiagonalDifference(matches[0], matches[1], strand)) < opts.cleanMaxDiag and 
				(diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[0], strand) < diagDrift )) { 
		onDiag[0] = true;
		nOnDiag++;
	}
	int m;
	
	for (int i = 1; i < matches.size() ; i++) {
		if (abs(DiagonalDifference(matches[i], matches[i-1], strand)) < opts.cleanMaxDiag and 
				(diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[i],strand) < diagDrift )) {	
			onDiag[i] = true;
			onDiag[i-1] = true;
			nOnDiag++;
		}
	}
	bool prevOnDiag = false;
	int  diagStart;
	int  Largest_ClusterNum = 0;
	//	ofstream diagSize("clustersize.tab");
	for (int i = 0; i < matches.size(); i++) {
		if (prevOnDiag == false and onDiag[i] == true) {
			diagStart = i;
		}
		if (prevOnDiag == true and onDiag[i] == false) {
			Largest_ClusterNum = max(Largest_ClusterNum, i - diagStart);
			//	diagSize << i - diagStart << "\t" << i << endl;
		}
		prevOnDiag = onDiag[i];
	}

	// Set the parameter minDiagCluster according to the value of Largest_ClusterNum
	// In this way, we won't lose small inversion.
	//cerr << "Largest_ClusterNum: " << Largest_ClusterNum << endl;
	/*
	if (Largest_ClusterNum < 10) {
		minDiagCluster = 1;
	}
	else if (Largest_ClusterNum < 20) {
		minDiagCluster = 2;
	} 
	else if (Largest_ClusterNum < 50) {
		minDiagCluster = 4;
	}
	else if (Largest_ClusterNum < 100) {
		minDiagCluster = 6;
	}
	else if (Largest_ClusterNum < 250) {
		minDiagCluster = 10;
	}
	else { // Largest_clusterNum >= 250 show obvious clusters
		minDiagCluster = 20;
	}
	*/
	int minDiagCluster=0.05*Largest_ClusterNum; // used to be 5, but if use smaller kmer, should be larger


	//minDiagCluster = (int) floor(Largest_ClusterNum/10);
	//cerr << "Largest_ClusterNum: " << Largest_ClusterNum << " minDiagCluster: " << minDiagCluster << endl;

	for (int i = 0; i < matches.size(); i++) {
		if (prevOnDiag == false and onDiag[i] == true) {
			diagStart = i;
		}
		if (prevOnDiag == true and onDiag[i] == false) {
			if (i - diagStart < minDiagCluster) {
				for (int j = diagStart; j < i; j++) {
					onDiag[j] = false;
				}
			}
		}
		prevOnDiag = onDiag[i];
	}

	int c   = 0;
	int pre = matches.size();
	for (int i=0; i < matches.size(); i++) {
		if (onDiag[i]) {
			matches[c] = matches[i]; c++;
		}
	}

	matches.resize(c);
}

class ClusterCoordinates {
 public:
	int start;
	int end;
	int strand;
	char *seq;
	GenomePos qStart, qEnd, tStart, tEnd;
	int chromIndex;	
	int coarseSubCluster;
	int clusterIndex;
	ClusterCoordinates() {
		clusterIndex=-1;
		qStart=-1;
		qEnd=0;
		tStart=-1;
		tEnd=0;
		seq=NULL;
		chromIndex=0;
		start=0;
		end=0;
		strand=-1;
		coarseSubCluster = -1;
	}
	bool EncompassesInRectangle(const ClusterCoordinates &b, const float frac) {
		int qovp=0;
		if (b.qStart >= qStart and b.qStart < qEnd) {
			qovp=min(qEnd, b.qEnd)-b.qStart;
		}
		else if (b.qEnd > qStart and b.qEnd <= qEnd) {
			qovp=b.qEnd-max(qStart, b.qStart);
		}
		else if (b.qStart <= qStart and b.qEnd > qEnd) {
			qovp=qEnd-qStart;
		}
		int tovp=0;
		if (b.tStart >= tStart and b.tStart < tEnd) {
			tovp=min(tEnd, b.tEnd)-b.tStart;
		}
		else if (b.tEnd > tStart and b.tEnd <= tEnd) {
			tovp=b.tEnd-max(tStart, b.tStart);
		}
		else if (b.tStart <= tStart and b.tEnd > tEnd) {
			tovp=tEnd-tStart;
		}
		//cerr << "encompass: " << float(qovp)/(b.qEnd-b.qStart) << "\t" << float(tovp)/(b.tEnd-b.tStart) << endl;
		return (float(qovp)/(b.qEnd-b.qStart) > frac) and (float(tovp)/(b.tEnd-b.tStart) > frac);
	}
	bool Encompasses(const ClusterCoordinates &b, const float frac) const {
		int qovp=0;
		if (b.qStart >= qStart and b.qStart < qEnd) {
			qovp=min(qEnd, b.qEnd)-b.qStart;
		}
		else if (b.qEnd > qStart and b.qEnd <= qEnd) {
			qovp=b.qEnd-max(qStart, b.qStart);
		}
		else if (b.qStart <= qStart and b.qEnd > qEnd) {
			qovp=qEnd-qStart;
		}
		int tovp=0;
		if (b.tStart >= tStart and b.tStart < tEnd) {
			tovp=min(tEnd, b.tEnd)-b.tStart;
		}
		else if (b.tEnd > tStart and b.tEnd <= tEnd) {
			tovp=b.tEnd-max(tStart, b.tStart);
		}
		else if (b.tStart <= tStart and b.tEnd > tEnd) {
			tovp=tEnd-tStart;
		}
		return ((float)qovp/(qEnd-qStart) > frac) or ((float)tovp/(tEnd-tStart) > frac);
	}

	bool OverlapsFracB(const ClusterCoordinates &b, float frac) const {
		int ovp=0;
		if (b.qStart >= qStart and b.qStart < qEnd) {
			ovp=min(qEnd, b.qEnd)-b.qStart;
		}
		else if (b.qEnd > qStart and b.qEnd < qEnd) {
			ovp=b.qEnd-max(qStart, b.qStart);
		}
		else if (b.qStart <= qStart and b.qEnd > qEnd) {
			ovp=qEnd-qStart;
		}
		float denomA=qEnd-qStart;
		float denomB=b.qEnd-b.qStart;
		//		cerr << "ovp: " << denomA << "\t" << denomB << "\t" << ovp << "\t" << ovp/denomB << endl;
		if ( ovp/denomB > frac) { return true; }
		else { return false; }
	}

	bool Overlaps(const ClusterCoordinates &b, float frac) const {
		int ovp=0;

		if (b.qStart >= qStart and b.qStart < qEnd) {
			ovp=min(qEnd, b.qEnd)-b.qStart;
		}
		else if (b.qEnd > qStart and b.qEnd < qEnd) {
			ovp=b.qEnd-max(qStart, b.qStart);
		}
		else if (b.qStart <= qStart and b.qEnd > qEnd) {
			ovp=qEnd-qStart;
		}
		float denomA=qEnd-qStart;
		float denomB=b.qEnd-b.qStart;
		if ( max(ovp/denomA, ovp/denomB) > frac) { return true; }
		else { return false; }
	}

	int Overlaps(const ClusterCoordinates &b) const {
		int ovp=0;

		if (b.qStart >= qStart and b.qStart < qEnd) {
			ovp=min(qEnd, b.qEnd)-b.qStart;
		}
		else if (b.qEnd > qStart and b.qEnd < qEnd) {
			ovp=b.qEnd-max(qStart, b.qStart);
		}
		else if (b.qStart <= qStart and b.qEnd > qEnd) {
			ovp=qEnd-qStart;
		}
		return ovp;
	}

	float OverlapsRate(const ClusterCoordinates &b) const {
		int ovp=0;

		if (b.qStart >= qStart and b.qStart < qEnd) {
			ovp=min(qEnd, b.qEnd)-b.qStart;
		}
		else if (b.qEnd > qStart and b.qEnd < qEnd) {
			ovp=b.qEnd-max(qStart, b.qStart);
		}
		else if (b.qStart <= qStart and b.qEnd > qEnd) {
			ovp=qEnd-qStart;
		}
		float denomA=qEnd-qStart;
		float denomB=b.qEnd-b.qStart;
		return max(ovp/denomA, ovp/denomB);
	}

	bool OverlapsOnRead(int & ReadLength, float frac) const {

		if (((float)(qEnd - qStart))/ReadLength > frac) {return true;}
		else { return false;}
	}

 ClusterCoordinates(int s,int e) : start(s), end(e) {
		qStart=qEnd=tStart=tEnd=strand=0;
		seq=NULL;
		chromIndex=0;
	}
 ClusterCoordinates(int s,int e, int st) : start(s), end(e), strand(st) {
		qStart=qEnd=tStart=tEnd=0;
		seq=NULL;
		chromIndex=0;
	}
  ClusterCoordinates(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, 
					int st) : start(s), end(e), strand(st), qStart(qs), qEnd(qe), tStart(ts), tEnd(te) {
		chromIndex=-1;
		seq=NULL;
	}
  ClusterCoordinates(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, 
					int st, int coarseSC) : start(s), end(e), strand(st), qStart(qs), qEnd(qe), tStart(ts), tEnd(te), coarseSubCluster(coarseSC) {
		chromIndex=-1;
		seq=NULL;
	}
};

class Cluster : public ClusterCoordinates {
 public:
	GenomePairs matches;
	vector<int> strands; // stores the strand of every GenomePair in matches
	vector<int> coarseSubCluster; // coarseSubCluster[i] means GenomePair i is from  the SubCluster[coarseSubCluster] 
	int64_t maxDiagNum;
	int64_t minDiagNum; // maxDiagNum and minDiagNum defines diagonal band boundary of the current cluster
	int coarse; 
	int Val; // Val stores the value of each Cluster;(using in SDP)
	int NumofAnchors0; // NumofAnchors0 stores the # of anchors of each splitcluster 
	vector<int> matchesLengths; // store the length of each anchor 
	bool refined; // refined == 0 means this Cluster has not been refined yet
	bool refinespace; // refinespace == 0 means this Cluster has not been add anchors in the step of RefineBtwnSpace;
	int outerCluster;
	int rank;
	Cluster() { refined=0; coarse=-1; rank=-1;}
 Cluster(int s, int e) : ClusterCoordinates(s,e) { coarse=-1; refined=0;}

 Cluster(int s, int e, int st) : ClusterCoordinates(s,e,st) { coarse=-1; refined=0;}

  Cluster(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, 
					int st) : ClusterCoordinates(s,e,qs,qe,ts,te,st) { coarse=-1; refined=0;} 
  Cluster(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, 
					int st, int cs) : ClusterCoordinates(s,e,qs,qe,ts,te,st) { coarse=cs; refined=0;} 
	
  Cluster(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, int st,
					GenomePairs::iterator gpBegin, GenomePairs::iterator gpEnd) : ClusterCoordinates(s,e,qs,qe,ts,te,st) {
		copy(gpBegin, gpEnd, back_inserter(matches));
		coarse=-1;
		maxDiagNum=0;
		minDiagNum=0;
		refined = 0;
		refinespace = 0;
		Val = 0;
		NumofAnchors0 = 0;
	}
  Cluster(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, int st,
					GenomePairs::iterator gpBegin, GenomePairs::iterator gpEnd, vector<int>::iterator stBegin, vector<int>::iterator stEnd) : ClusterCoordinates(s,e,qs,qe,ts,te,st) {
		copy(gpBegin, gpEnd, back_inserter(matches));
		copy(stBegin, stEnd, back_inserter(strands));
		coarse=-1;
		maxDiagNum=0;
		minDiagNum=0;
		Val = 0;
		NumofAnchors0 = 0;
	}
  Cluster(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te, int st, int coa) {
		qStart = qs;
		qEnd = qe;
		tStart = ts;
		tEnd = te;
		strand = st;
		coarse = coa;
		Val = 0;
		NumofAnchors0 = 0;
	}
  Cluster(int st, GenomePairs::iterator gpBegin, GenomePairs::iterator gpEnd, vector<int>::iterator stBegin, vector<int>::iterator stEnd) {
  		strand = st;
		copy(gpBegin, gpEnd, back_inserter(matches));
  }
	
  bool OverlapsPrevious(const Cluster &prev) {
		//
		// Assume clusters are sorted by target.
		//
		if (prev.strand == strand and 
				prev.tEnd >= tStart and
				prev.tStart < tStart and 
				prev.qEnd >= qStart and
				prev.qStart < qStart) {
			return true;
		}
		else {
			return false;
		}
	}
	void UpdateBoundaries(const Cluster &rhs) {
		tEnd   = max(tEnd, rhs.tEnd);
		tStart = min(tStart, rhs.tStart);
		qEnd   = max(qEnd, rhs.qEnd);
		qStart = min(qStart, rhs.qStart);
	}
	int size() const {
		return matches.size();
	}

	int operator<(const Cluster &rhs) const {
		if (strand != rhs.strand) {
			return strand != rhs.strand;
		}
		else if (tStart != rhs.tStart) {
			return tStart < rhs.tStart;
		}
		else {
			return qStart < rhs.qStart;
		}
	}

	void SetClusterBoundariesFromMatches (Options &opts) {
		qStart = matches[0].first.pos;
		qEnd = qStart + opts.globalK;
		tStart = matches[0].second.pos;
		tEnd = tStart + opts.globalK;
		for (int i = 1; i < matches.size(); i++) {
			tEnd = max(tEnd, matches[i].second.pos + opts.globalK);
			tStart = min(tStart, matches[i].second.pos);
			qEnd = max(qEnd, matches[i].first.pos + opts.globalK);
			qStart = min(qStart, matches[i].first.pos);
		}
	}		
};

class OrderClusterBySize {
 public:
	int operator()(const Cluster &a, const Cluster &b) {
		return a.size() > b.size();
	}
};

class ClusterOrder {
 public:
	vector<Cluster> *clusters;
	vector<int> index;
	int orderType;
 ClusterOrder(vector<Cluster> *c, int t=0) : clusters(c), orderType(t) {
		index.resize(clusters->size());
		for (int i=0;i<index.size();i++) { index[i]=i;}
		Sort();
	}
	
	//
	// Cartesian sort of clusters.
	//
	int operator()(const int i, const int j) {
		if (orderType == 1) { 
			return (*clusters)[i].size() > (*clusters)[j].size();
		}
		else {
			assert((*clusters)[i].strand == 0 or (*clusters)[i].strand == 1);
			assert((*clusters)[j].strand == 0 or (*clusters)[j].strand == 1);
			
			if ((*clusters)[i].tStart != (*clusters)[j].tStart) {
				return (*clusters)[i].tStart < (*clusters)[j].tStart;
			}
			else {
				return (*clusters)[i].qStart < (*clusters)[j].qStart;
			}
		}
	}
	template<typename t>
		void Sort() {
		sort(index.begin(), index.end(), t());
	}

	void Sort() {
		sort(index.begin(), index.end(), *this);
	}
	Cluster & operator[](int i) {
		return (*clusters)[index[i]];
	}
	int size() {
		return index.size();
	}
};


class Clusters_valueOrder {
 public:
	vector<float> *clusters_value;
	vector<int> index;
	
 	Clusters_valueOrder(vector<float> *c) : clusters_value(c) {
		index.resize((*clusters_value).size());
		for (int i=0;i<index.size();i++) { index[i]=i;}
		Sort();
	}

	int operator()(const int i, const int j) {
		assert(i < clusters_value->size());
		assert(j < clusters_value->size());
		return (*clusters_value)[i] > (*clusters_value)[j];			
	}

	void Sort() {
		sort(index.begin(), index.end(), *this);
	}

	float & operator[](int i) {
		return (*clusters_value)[index[i]];
	}
	
	int size() {
		return index.size();
	}
};


class LogCluster {
 public:
 	vector<Cluster> SubCluster;
 	Cluster * Hp;
 	bool ISsecondary; // ISsecondary == 1 means this is a secondary chain. Otherwise it's a primary chain
 	int primary; // When ISsecondary == 1, primary stores the index of the primary chain in vector<LogCluster>
 	vector<int> secondary; // When ISsecondary == 0, secondary stores the indices of the secondary chains
 	int coarse; // coarse means this LogCluster stores information about refinedCluster[coarse]
 	bool direction; // direction means how the read is mapped to the reference
 	bool split; // split == 1 means this chain has been split
 	int main; // If a read is splitted, "main" stores the index of the main alignment part

 	LogCluster () {
 		ISsecondary = 0;
 		primary = -1;
 		coarse = -1;
 		direction = 0;
 		split = 0;
 	};
 	~LogCluster() {};
 	
 	void setHp(Cluster & H) {
 		Hp = & H;
 	}
 	void SetCoarse () {
 		coarse = SubCluster[0].coarse;
 	}

 	void SetSubClusterBoundariesFromMatches (Options &opts) {
		// set the boundaries for SubCluster[i] -- the current last one in SubCluster
 		int i = SubCluster.size() - 1; // the last one in SubCluster
		for (int is = SubCluster[i].start; is < SubCluster[i].end; ++is) {

			if (is == SubCluster[i].start) {
		 		SubCluster[i].tStart = Hp->matches[is].second.pos;
		 		SubCluster[i].qStart = Hp->matches[is].first.pos;					
			}
			SubCluster[i].tEnd   = max(SubCluster[i].tEnd, Hp->matches[is].second.pos + opts.globalK);
			SubCluster[i].tStart = min(SubCluster[i].tStart, Hp->matches[is].second.pos);
			SubCluster[i].qEnd   = max(SubCluster[i].qEnd, Hp->matches[is].first.pos + opts.globalK);
			SubCluster[i].qStart = min(SubCluster[i].qStart, Hp->matches[is].first.pos); 	
		}
 	}
};


template<typename Tup>
void PrintDiagonal(vector<pair<Tup, Tup> > &matches, int strand=0) {
	for (int m=1; m < matches.size(); m++) {
		long d=DiagonalDifference(matches[m], matches[m-1], strand);
		cerr << matches[m-1].first.pos << "\t" << matches[m].first.pos << "\t" << matches[m-1].second.pos << "\t" << matches[m].second.pos << "\t" << d << endl;
	}
}
 
template<typename Tup>
long GetDiag(pair<Tup, Tup> &match, int strand, Options &opts) {
	if (strand == 0) {
		return (long) match.second.pos - (long) ceil(opts.slope*match.first.pos);
	}
	else {
		return (long) ceil(opts.slope*match.first.pos) + (long) match.second.pos;
	}
}
template<typename Tup>
void EstimateDiagonalSlope(vector<pair<Tup, Tup> > &matches, int s, int e, Options &opts, Genome &genome, int ri, int outerIteration, int strand=0) {
	if (e-s == 1) return;
	int cs = s;
	int ce = e;
	long driftOrigin = GetDiag(matches[cs], strand, opts); // Change to a one with rate as 1
	while (cs < e) {
		ce = cs+1;

		while (ce < e and (abs(DiagonalDifference(matches[ce], matches[ce-1], strand)) < opts.maxDiag ) 
					  and (GapDifference(matches[ce], matches[ce-1]) < opts.maxGap) 
					  and DiagonalDrift(driftOrigin, matches[ce], strand) < opts.maxDrift) { // quite large maxGap
			ce++;
		}	
		if (ce-cs >= opts.minTightCluster) {
			if (strand == 0) {
				opts.slope = (float)(((double)matches[cs].second.pos - (double)matches[ce-1].second.pos) / ((double)matches[cs].first.pos - (double)matches[ce-1].first.pos));
			}
			else {
				opts.slope = (float)(((double)matches[cs].second.pos - (double)matches[ce-1].second.pos) / ((double)matches[ce-1].first.pos - (double)matches[cs].first.pos));
			}
			// cerr << "outerIteration: " << outerIteration << endl; 
			// cerr << "matches[" << cs << "]: " << matches[cs].first.pos << ", " << matches[cs].second.pos << endl;
			// cerr << "matches[" << ce << "]: " << matches[ce-1].first.pos << ", " << matches[ce-1].second.pos << endl;
			// cerr << "opts.slope: " << opts.slope << endl;
			break;
		}
		cs = ce;
		driftOrigin = GetDiag(matches[cs], strand, opts);
	}
}

template<typename Tup>
void StoreFineClusters(int ri, vector<pair<Tup, Tup> > &matches, vector<Cluster> &clusters, Options &opts, int s, int e, 
											 Genome &genome,
											 GenomePos readLength, 
											// interval_map<GenomePos, int> &xIntv, 
											// interval_map<GenomePos, int> &yIntv,
											 int strand=0, int outerIteration=0) {
	//
	// StoreFineClusters based on unique part
	// The first step: store the number of matches corresponding to each minimizer from the read; -- a linear pass
	// matches are sorted in cartesian order;
	//
	if (e-s == 1) return;
	vector<int> match_num;
	vector<int> pos_start;
	int oc=1;
	int us=s;
	for (int i=s+1; i<e; i++) {
		if (matches[i].first.pos == matches[i-1].first.pos) {
			oc++;
		}
		else {
			match_num.push_back(oc);
			pos_start.push_back(us);
			us=i;
			oc = 1;
		}
		if (i == e-1) {
			match_num.push_back(oc);
			pos_start.push_back(us);				
		}
	}
	//
	// Find stretches of "1" on a diagonal in match_num; -- a linear pass
	//
	int u_start=0, u_end=0, u_maxstart=0, u_maxend=0;
	int max_pos = 0;
	bool reset=0; // reset == 0 means a new stretch of "1";
	vector<int> Start;
	vector<int> End;
	if (match_num.size() == 1) {
		assert(pos_start[0] == s);
		if (match_num[0] == 1) {
			u_maxstart = 0; 
			u_maxend = 1;	
			Start.push_back(u_maxstart);
			End.push_back(u_maxend);		
		}
	}
	else {
		int k=0;
		while (k<match_num.size()-1) {
			if (match_num[k] != 1) {
				k++;
				continue;
			}
			else if (reset == 0) {
				u_start = k;
				u_end = k+1;					
			}
			if (match_num[k+1] == match_num[k] 
					and abs(DiagonalDifference(matches[pos_start[k+1]], matches[pos_start[k]], strand)) < opts.maxDiag/2) { 
				u_end = k+2;
				k++;
				reset = 1;
			}
			else {
				if ((u_maxstart == 0 and u_maxend == 0) or (u_maxend-u_maxstart < u_end-u_start)) {
					u_maxstart = u_start;
					u_maxend = u_end;		
					max_pos=Start.size(); 
				}
				if (Start.size() > 0) assert(u_start != Start.back());
				Start.push_back(u_start);
				End.push_back(u_end);		
				assert(u_maxstart == Start[max_pos]);
				assert(u_maxend == End[max_pos]);
				//if ((u_maxstart == 0 and u_maxend == 0) or (u_maxend-u_maxstart < u_end-u_start)) max_pos=Start.size() - 1; 
				k++;
				reset = 0;
			}

		}	
		if (reset == 1) {

			if (k == match_num.size()-1) {
				if (u_maxstart == 0 and u_maxend == 0) { // the whole rough cluster is a unique linear part
					u_start = 0;
					u_end = match_num.size();
					max_pos=Start.size(); 
					u_maxstart = u_start;
					u_maxend = u_end;
				}	
				else if (u_maxend-u_maxstart < match_num.size()-u_start) {
					u_end = match_num.size();	
					max_pos=Start.size(); 
					u_maxstart = u_start;
					u_maxend = u_end;	
				}
				if (Start.size() > 0) assert(u_start != Start.back());
				Start.push_back(u_start);
				End.push_back(u_end);				
				assert(u_maxstart == Start[max_pos]);
				assert(u_maxend == End[max_pos]);
			}

			// if (u_maxstart == 0 and u_maxend == 0) { // the whole rough cluster is a unique linear part
			// 	u_start = 0;
			// 	u_end = match_num.size();
			// 	max_pos=Start.size(); 
			// }	
			// else if (k == match_num.size()-1 and u_maxend-u_maxstart < match_num.size()-u_start) {
			// 	u_end = match_num.size();	
			// 	max_pos=Start.size(); 	
			// }
			// u_maxstart = u_start;
			// u_maxend = u_end;
			// if (Start.size() > 0) assert(u_start != Start.back());
			// Start.push_back(u_start);
			// End.push_back(u_end);	
			// assert(u_maxstart == Start[max_pos]);
			// assert(u_maxend == End[max_pos]);	
			//if ((k == match_num.size()-1 and u_maxend-u_maxstart < match_num.size()-u_start) or (u_maxstart == 0 and u_maxend == 0)) max_pos=Start.size() - 1; 			
		}
	}
	//
	// Decide the u_minDiag and u_maxDiag of the unique stretches of anchors
	//
	int c_s = pos_start[u_maxstart], c_e = pos_start[u_maxend-1]+1;

	if (opts.dotPlot) {
		ofstream uniquematch("UniqueMatch.tab", std::ofstream::app);
		for (int h = 0; h < Start.size(); h++) {
			for (int he = Start[h]; he < End[h]; he++) {
				if (strand == 0) {
					uniquematch << matches[pos_start[he]].first.pos << "\t"
						  << matches[pos_start[he]].second.pos << "\t"
						  << matches[pos_start[he]].first.pos + opts.globalK << "\t"
						  << matches[pos_start[he]].second.pos + opts.globalK << "\t"
						  << outerIteration << "\t"
						  << opts.slope << "\t"
						  << h << "\t"
						  << End[h] - Start[h] << "\t"
						  << strand << endl;
				}
				else {
					uniquematch << matches[pos_start[he]].first.pos << "\t"
						  << matches[pos_start[he]].second.pos + opts.globalK << "\t"
						  << matches[pos_start[he]].first.pos + opts.globalK << "\t"
						  << matches[pos_start[he]].second.pos<< "\t"
						  << outerIteration << "\t"
						  << opts.slope << "\t"	
						  << h << "\t"		
						  << End[h] - Start[h] << "\t"		  
						  << strand << endl;					
				}					
			}
		}
		uniquematch.close();
	}	

	if (c_e - c_s >= opts.minUniqueStretchNum and 
			matches[c_e-1].first.pos + opts.globalK - matches[c_s].first.pos >= opts.minUniqueStretchDist) { 
		clusters.push_back(Cluster(0, 0, strand)); 

		if (c_e - c_s == e - s) { // no need to extend - the whole rough cluster is a unique linear part
			for (int i = c_s; i < c_e; i++) {
				clusters.back().matches.push_back(matches[i]);				
			}	
		}
		else {
			//
			// extend the unique stretches from the left; -- a linear pass in vector<int> Start and End;
			//
			list<int> StretchOfOne;
			int prev_anchor = c_s;
			int prev_stretch = max_pos;
			int addup = 0;
			if (max_pos >= 0) {
				StretchOfOne.push_back(max_pos);
				int i = max_pos - 1;
				while (i >= 0) {
					int i_m = pos_start[End[i]-1];
					
					// cerr << "Diag: " << prev_stretch << " and " << i << " " << abs(DiagonalDifference(matches[i_m], matches[prev_anchor], strand)) << endl;
					// cerr << "matches[i_m].first.pos: " << matches[i_m].first.pos << " matches[i_m].second.pos: " << matches[i_m].second.pos << endl;
					// cerr << "matches[prev_anchor].first.pos: " << matches[prev_anchor].first.pos << " matches[prev_anchor].second.pos: " << matches[prev_anchor].second.pos << endl;
					// cerr << "prev_anchor: " << prev_anchor << endl;
				
					if (/*End[i] - Start[i] >= (c_e-c_s) / 3 and */ End[i] - Start[i] >= 15) { // stretches of above 15 unqiue anchors are pretty trustworthy
						addup = 1000000; //50000
					}
					else addup = 0;
					if ( (End[i] - Start[i] >= (c_e-c_s) / 20 or End[i] - Start[i] >= 15) 
							and abs(DiagonalDifference(matches[i_m], matches[prev_anchor], strand)) < opts.maxDiag+addup) {
						assert(i < StretchOfOne.back());
						StretchOfOne.push_back(i);
						prev_anchor = pos_start[Start[i]];
						prev_stretch = i;
					}
					i--;
				}
			}
			//
			// extend the unique stretches from the right; -- a linear pass in vector<int> Start and End;
			//
			prev_anchor = c_e - 1;
			prev_stretch = max_pos;
			if (max_pos < Start.size()) {
				int i = max_pos + 1;
				while (i < Start.size()) {
					int i_m = pos_start[Start[i]];
					// cerr << "Diag: " << prev_stretch << " and " << i << " " << abs(DiagonalDifference(matches[i_m], matches[prev_anchor], strand)) << endl;
					// cerr << "matches[i_m].first.pos: " << matches[i_m].first.pos << " matches[i_m].second.pos: " << matches[i_m].second.pos << endl;
					// cerr << "matches[prev_anchor].first.pos: " << matches[prev_anchor].first.pos << " matches[prev_anchor].second.pos: " << matches[prev_anchor].second.pos << endl;
					// cerr << "prev_anchor: " << prev_anchor << endl;
					
					if (/*End[i] - Start[i] >= (c_e-c_s) / 3 and */ End[i] - Start[i] >= 15) {
						addup = 1000000;
					}
					else addup = 0;
					if ((End[i] - Start[i] >= (c_e-c_s) / 20 or End[i] - Start[i] >= 15) 
							and abs(DiagonalDifference(matches[i_m], matches[prev_anchor], strand)) < opts.maxDiag+addup) {
						assert(i > StretchOfOne.front());
						StretchOfOne.push_front(i);
						prev_anchor = pos_start[End[i]-1];
						prev_stretch = i;
					}
					i++;
				}
			}
			//
			// extend anchors between every two stretch of unqiue matches; Do not forget the ends;
			//
			assert(StretchOfOne.size() > 0);
			prev_stretch = -1;
			int p_s = 0, p_e = 0;
			for (list<int>::reverse_iterator it=StretchOfOne.rbegin(); it!=StretchOfOne.rend(); ++it) {
				
				vector<int> Cluster_index;
				c_s = pos_start[Start[*it]]; c_e = pos_start[End[*it]-1]+1;

				if (it == StretchOfOne.rbegin()) {
					p_s = s;
					p_e = pos_start[Start[*it]];
				}
				else {
					//p_s = pos_start[End[prev_stretch]-1]; ///??? -1
					p_s = pos_start[End[prev_stretch]]; ///??? -1
					p_e = pos_start[Start[*it]];
				}
				prev_stretch = *it;
				assert(p_s >= s and p_e <= e and p_s <= p_e and p_e == c_s);

				int prev_match = c_s;
				for (int si = p_e - 1; si >= p_s; si--) {
					if (abs(DiagonalDifference(matches[si], matches[prev_match], strand)) < opts.maxDiag) { //opts.maxDiag-200
						Cluster_index.push_back(si);
						prev_match = si;
					}
				}
				for (vector<int>::reverse_iterator ci = Cluster_index.rbegin(); ci != Cluster_index.rend(); ++ci) {
					clusters.back().matches.push_back(matches[*ci]);
					//
					// for debug
					//
					if (clusters.back().matches.size() >= 2) {
						int mb = clusters.back().matches.size() - 1;
						assert(clusters.back().matches[mb].first.pos >= clusters.back().matches[mb-1].first.pos);
						if (clusters.back().matches[mb].first.pos == clusters.back().matches[mb-1].first.pos) {
							assert(clusters.back().matches[mb].second.pos >= clusters.back().matches[mb-1].second.pos);
						}
					}
				}
				// insert the unique anchors
				for (int si = c_s; si < c_e; si++) {
					//Cluster_index.push_front(si);
					clusters.back().matches.push_back(matches[si]);
				}
				if (opts.dotPlot) {
					ofstream uniquepart("UniquePart.tab", std::ofstream::app);
					for (int si = c_s; si < c_e; si++) {
						if (strand == 0) {
							uniquepart << matches[si].first.pos << "\t"
								  << matches[si].second.pos << "\t"
								  << matches[si].first.pos + opts.globalK << "\t"
								  << matches[si].second.pos + opts.globalK << "\t"
								  << outerIteration << "\t"
								  << opts.slope << "\t"
								  << strand << endl;
						}
						else {
							uniquepart << matches[si].first.pos << "\t"
								  << matches[si].second.pos + opts.globalK << "\t"
								  << matches[si].first.pos + opts.globalK << "\t"
								  << matches[si].second.pos<< "\t"
								  << outerIteration << "\t"
								  << opts.slope << "\t"					  
								  << strand << endl;					
						}							
					}
					uniquepart.close();
				}						
				//
				// extend the right end;
				//
				if ((it != StretchOfOne.rend()) and (next(it) == StretchOfOne.rend())) {
					p_s = pos_start[End[*it]-1]+1;
					p_e = e;
					assert(p_s >= s and p_e <= e and p_s <= p_e);
					prev_match = c_e - 1;
					for (int si = p_s; si < p_e; si++) {
						if (abs(DiagonalDifference(matches[si], matches[prev_match], strand)) < opts.maxDiag) {//opts.maxDiag-200
							clusters.back().matches.push_back(matches[si]);
							prev_match = si;
							if (clusters.back().matches.size() >= 2) {
								int mb = clusters.back().matches.size() - 1;
								assert(clusters.back().matches[mb].first.pos >= clusters.back().matches[mb-1].first.pos);
								if (clusters.back().matches[mb].first.pos == clusters.back().matches[mb-1].first.pos) {
									assert(clusters.back().matches[mb].second.pos >= clusters.back().matches[mb-1].second.pos);
								}
							}
						}
					}						
				}
			}
		}
		//
		// Decide the cooridinates of the cluster;
		//
		clusters.back().SetClusterBoundariesFromMatches(opts);
		clusters.back().chromIndex = ri;

		if (opts.dotPlot) {
			ofstream fineclusters("fineclusters_byunique.tab", std::ofstream::app);
			if (clusters.size() > 0) {
				for (int h = 0; h < clusters.back().matches.size(); h++) {
					if (strand == 0) {
					   fineclusters << clusters.back().matches[h].first.pos << "\t"
							  << clusters.back().matches[h].second.pos << "\t"
							  << clusters.back().matches[h].first.pos + opts.globalK << "\t"
							  << clusters.back().matches[h].second.pos + opts.globalK << "\t"
							  << outerIteration << "\t"
							  << clusters.size() - 1 << "\t"
							  << strand << endl;
					}
					else {
						fineclusters << clusters.back().matches[h].first.pos << "\t"
							  << clusters.back().matches[h].second.pos + opts.globalK << "\t"
							  << clusters.back().matches[h].first.pos + opts.globalK << "\t"
							  << clusters.back().matches[h].second.pos<< "\t"
							  << outerIteration << "\t"
							  << clusters.size() - 1 << "\t"
							  << strand << endl;					
					}
				}				
			}
			fineclusters.close();
		}	
		
		if (clusters.back().matches.size() <= opts.minClusterSize) {
			clusters.pop_back();
		}
	}
	else {
		return;
	}
}

void 
SplitRoughClustersWithGaps(GenomePairs &matches, Cluster &OriginalClusters, vector<Cluster> &split, Options &opts, int strand=0) {

	if (OriginalClusters.end - OriginalClusters.start == 0) {
		return;
	}
	int split_cs = OriginalClusters.start;
	GenomePos split_qStart = matches[split_cs].first.pos;
	GenomePos split_tStart = matches[split_cs].second.pos;
	GenomePos split_qEnd = split_qStart + opts.globalK; 
	GenomePos split_tEnd = split_tStart + opts.globalK;

	for (int m=OriginalClusters.start+1; m < OriginalClusters.end; m++) {
		int gap=GapDifference(matches[m], matches[m-1]);
		//int diag_gap=DiagonalDifference(matches[m], matches[m-1], strand);

		if (gap > opts.maxGap /*or diag_gap > opts.maxGap*/) {
			// cerr << "GAP: " << gap << "\t" << m << "\t" << clusters[c].matches.size() << "\t" << clusters[c].matches[m].second.pos - clusters[c].matches[m-1].second.pos << "\t" << clusters[c].matches[m].first.pos - clusters[c].matches[m-1].first.pos << endl;
			if ( (m-split_cs) > opts.minClusterSize) split.push_back(Cluster(split_cs, m, split_qStart, split_qEnd, split_tStart, split_tEnd, OriginalClusters.strand));
			split_qStart = matches[m].first.pos;
			split_tStart = matches[m].second.pos;
			split_qEnd = split_qStart + opts.globalK; 
			split_tEnd = split_tStart + opts.globalK;
			split_cs = m;
		}
		else {
			split_qStart = min(split_qStart, matches[m].first.pos);
			split_tStart = min(split_tStart, matches[m].second.pos);
			split_qEnd = max(split_qEnd, matches[m].first.pos + opts.globalK);
			split_tEnd = max(split_tEnd, matches[m].second.pos + opts.globalK);
		}
	}
	int last=OriginalClusters.end;
	split.push_back(Cluster(split_cs, last, split_qStart, split_qEnd, split_tStart, split_tEnd, OriginalClusters.strand));
}

void SplitClustersWithGaps(vector<Cluster> &clusters, vector<Cluster> &split, Options &opts) {
	int curSplit=-1;

	for (int c=0; c < clusters.size(); c++) {
		split.push_back(Cluster());
		curSplit++;
		split[curSplit].tStart = clusters[c].tStart;
		split[curSplit].qStart = clusters[c].qStart;

		if (clusters[c].matches.size() == 0) {
			continue;
		}
		split[curSplit].matches.push_back(clusters[c].matches[0]);
		for (int m=1; m < clusters[c].matches.size(); m++) {
			int gap=GapDifference(clusters[c].matches[m], clusters[c].matches[m-1]);

			if (gap > 500) {
				// s				cerr << "GAP: " << gap << "\t" << m << "\t" << clusters[c].matches.size() << "\t" << clusters[c].matches[m].second.pos - clusters[c].matches[m-1].second.pos << "\t" << clusters[c].matches[m].first.pos - clusters[c].matches[m-1].first.pos << endl;
				split[curSplit].qEnd = clusters[c].matches[m-1].first.pos + opts.globalK;
				split[curSplit].tEnd = clusters[c].matches[m-1].second.pos + opts.globalK;

				split.push_back(Cluster());
				curSplit++;
				split[curSplit].qStart = clusters[c].matches[m].first.pos;
				split[curSplit].tStart = clusters[c].matches[m].second.pos;
			}
		}
		int last=clusters[c].matches.size();
		split[curSplit].qEnd = clusters[c].matches[last-1].first.pos + opts.globalK;
		split[curSplit].tEnd = clusters[c].matches[last-1].second.pos + opts.globalK;
	}
}


template<typename Tup>
void StoreDiagonalClusters(vector<pair<Tup, Tup> > &matches, vector<Cluster> &clusters, Options &opts, int s, int e, int strand=0) {

	// int maxGap = -1, maxDiag = -1;
	// if (rough == false) { maxGap = opts.maxGapBtwnAnchors;}
	// else { maxDiag = opts.maxDiag;} 

	int i;
	int cs = s, ce = e;
	int64_t dd,absdd;
	
	while (cs < e) {
		ce = cs+1;
		GenomePos qStart=matches[cs].first.pos, 
				  qEnd=matches[cs].first.pos + opts.globalK, 
			      tStart=matches[cs].second.pos, 
			      tEnd=matches[cs].second.pos + opts.globalK;

		// (TODO)Jingwen: Delete (opts.maxGap == -1 or GapDifference(matches[ce], matches[ce-1]) < opts.maxGap) in the below
		while (ce < e and abs(DiagonalDifference(matches[ce], matches[ce-1], strand)) < opts.maxDiag
						//and (GapDifference(matches[ce], matches[ce-1]) < opts.maxGap)
					  //and (maxGap == -1 or GapDifference(matches[ce], matches[ce-1]) < maxGap)
					  ) {

			qStart = min(qStart, matches[ce].first.pos);
			qEnd   = max(qEnd, matches[ce].first.pos + opts.globalK);
			tStart = min(tStart, matches[ce].second.pos);
			tEnd   = max(tEnd, matches[ce].second.pos + opts.globalK);
			/*		cerr << "rc:\t" << clusters.size() << "\t" << matches[ce].first.pos << "\t" << matches[ce].second.pos << "\t" 
					 << GetDiag(matches[ce], strand) << "\t" <<
					 abs(DiagonalDifference(matches[ce], matches[ce-1], strand))
					 << endl;
			*/
			ce++;
		}	
		
		if (ce - cs >= opts.minClusterSize and qEnd - qStart >= opts.minClusterLength and tEnd - tStart >= opts.minClusterLength) {
			//if (rough == true) {
				clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand));
			//}
			// else {
			// 	if (lite == false ) {
			// 		clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand, matches.begin() + cs, matches.begin()+ce));
			// 	}
			// 	else {
			// 		clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand));
			// 	}
			// }
		}
		cs=ce;
	}	
}


bool sign(int val) {
	if (val >= 0) return true;
	return false;
}

// TODO(Jingwen): delete this later
template<typename Tup>
void RemovePairedIndels(vector<pair<Tup, Tup>> & matches, vector<Cluster> & clusters, Options & opts) {
	if (clusters.size() < 3) { return;}
	vector<bool> remove(clusters.size(), false);

	for (int c = 1; c < clusters.size() - 1; c++) {
		GenomePos prevQEnd = clusters[c-1].qEnd;
		GenomePos prevTEnd = clusters[c-1].tEnd;
		GenomePos qStart = clusters[c].qStart;
		GenomePos tStart = clusters[c].tStart;
		GenomePos qEnd = clusters[c].qEnd;
		GenomePos tEnd = clusters[c].tEnd;

		GenomePos nextQStart = clusters[c+1].qStart;
		GenomePos nextTStart = clusters[c+1].tStart;

		int prevGap = (int)(qStart - prevQEnd) - (int)(tStart - prevTEnd);
		int nextGap = (int)(nextQStart - qEnd) - (int)(nextTStart - tEnd);

		if (sign(prevGap) != sign(nextGap) and
				abs(prevGap) + abs(nextGap) > abs(prevGap + nextGap)  ) { //(clusters[c].end +opts.k - clusters[c].start)) {
			remove[c] = true;
			for (int ci=clusters[c].start; ci < clusters[c].end; ci++) {
				matches[ci].first.pos = -1;
			}
		} 
	}	
	int m=0;

	for (int i=0; i < matches.size(); i++) {
		if (matches[i].first.pos != -1) {
			matches[m] = matches[i];
			m++;
		}
	}
	matches.resize(m);
}


//
// This function removes paired indels for chain after 2nd SDP
void RemovePairedIndels(GenomePos qAlnStart, GenomePos tAlnStart, GenomePos qAlnEnd, GenomePos tAlnEnd, vector<unsigned int> & matches, GenomePairs & Pairs, Options &opts) {
	unsigned int nMatches = matches.size();
	if ( nMatches < 3)   { return;}
	vector<bool> remove(nMatches, false);
	GenomePos prevQEnd, prevTEnd, qStart, tStart, qEnd, tEnd;
	GenomePos nextQStart;
	GenomePos nextTStart;

	for (unsigned int c = 0; c < nMatches ; c++) {

		if (c == 0) {
			prevQEnd = qAlnStart;
			prevTEnd = qAlnEnd;
		}		
		else {
			prevQEnd = Pairs[matches[c-1]].first.pos + opts.globalK;
			prevTEnd = Pairs[matches[c-1]].second.pos + opts.globalK;
		}
		qStart   = Pairs[matches[c]].first.pos;
		tStart   = Pairs[matches[c]].second.pos;
		qEnd     = Pairs[matches[c]].first.pos + opts.globalK;
		tEnd     = Pairs[matches[c]].second.pos + opts.globalK;

		if (c < nMatches-1) {
			nextQStart = Pairs[matches[c+1]].first.pos;
			nextTStart = Pairs[matches[c+1]].second.pos;
		}
		else {
			nextQStart = qAlnEnd;
			nextTStart = tAlnEnd;
		}


		int prevGap = (int)(qStart-prevQEnd) - (int)(tStart-prevTEnd);
		int nextGap = (int)(nextQStart - qEnd) - (int)(nextTStart - tEnd);

		if (sign(prevGap) != sign(nextGap) and
				abs(prevGap) + abs(nextGap) >  abs(prevGap + nextGap)  ) {
			remove[c] = true;
		} 
	}	
	int m=0;

	for (int i=0; i < nMatches; i++) {
		if (remove[i] == false) {
			matches[m] = matches[i];
			m++;
		}
	}
	matches.resize(m);
}


void SetClusterBoundariesFromSubCluster(Cluster &cluster, Options &opts, LogCluster &logCluster) {
	for (int i = 0; i < logCluster.SubCluster.size(); ++i) {
		cluster.tEnd = max(cluster.tEnd, logCluster.SubCluster[i].tEnd);
		cluster.tStart = min(cluster.tStart, logCluster.SubCluster[i].tStart);
		cluster.qEnd = max(cluster.qEnd, logCluster.SubCluster[i].qEnd);
		cluster.qStart = min(cluster.qStart, logCluster.SubCluster[i].qStart);			
	}
}

void SetCoarseFromSubClusters (Cluster & cluster, const LogCluster &logCluster) {
	for (int i = 0; i < logCluster.SubCluster.size(); ++i) {
		for (int j = logCluster.SubCluster[i].start; j < logCluster.SubCluster[i].end; ++j) {
			cluster.coarseSubCluster[j] = i;
		}
	}
}

#endif
