#ifndef CLUSTERING_H_
#define CLUSTERING_H_
#include "Options.h"
#include "seqan/seeds.h"

template<typename Tup>
int64_t DiagonalDifference(Tup &a, Tup &b) {
		int64_t aDiag = a.first.pos - a.second.pos, 
			bDiag= b.first.pos - b.second.pos;
		return aDiag - bDiag;
}
template<typename Tup>
int64_t GapDifference(Tup &a, Tup &b) {
	int64_t aDiff = abs((int)b.first.pos - (int)a.first.pos);
	int64_t bDiff = abs((int)b.second.pos - (int)a.second.pos);
	return max(aDiff, bDiff);
}

template<typename Tup>
int DiagonalDrift(int curDiag, Tup &t) {
	int drift= abs(curDiag - ((int)t.first.pos - (int)t.second.pos));
	return drift;
}

template<typename Tup>
void CleanOffDiagonal(vector<pair<Tup, Tup> > &matches, Options &opts, int diagOrigin=-1, int diagDrift=-1) {
	if (matches.size() == 0) {
		return;
	}
	
	vector<bool> onDiag(matches.size(), false);
	
	if (matches.size() > 1 and abs(DiagonalDifference(matches[0], matches[1])) < opts.maxDiag and 
				(diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[0]) < diagDrift )) {
		onDiag[0] = true;
	}
	int m;
	for (int i = 1; i < matches.size() ; i++) {
		if (abs(DiagonalDifference(matches[i], matches[i-1])) < opts.maxDiag and 
				(diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[i]) < diagDrift ) ) {	
			onDiag[i] = true;
		}
	}
	bool prevOnDiag = false;
	int  diagStart;
	for (int i = 0; i < matches.size(); i++) {
		if (prevOnDiag == false and onDiag[i] == true) {
			diagStart = i;
		}
		if (prevOnDiag == true and onDiag[i] == false) {
			if (i - diagStart < opts.minDiagCluster) {
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
	ClusterCoordinates() {
		qStart=-1;
		qEnd=0;
		tStart=-1;
		tEnd=0;
		seq=NULL;
		chromIndex=0;
		start=0;
		end=0;
		strand=-1;
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
		float denomB=b.qEnd- b.qStart;
		if ( max(ovp/denomA, ovp/denomB)  > frac) { return true; }
		else { return false; }
	}

 ClusterCoordinates(int s,int e) : start(s), end(e) {
		qStart=qEnd=tStart=tEnd=strand=0;
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
};

class Cluster : public ClusterCoordinates {
 public:
	GenomePairs matches;
	int coarse;
	Cluster() {}
  Cluster(int s, int e) : ClusterCoordinates(s,e) { coarse=0;}


  Cluster(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, 
					int st) : ClusterCoordinates(s,e,qs,qe,ts,te,st) { coarse=0;} 
	
  Cluster(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, int st,
					GenomePairs::iterator gpBegin, GenomePairs::iterator gpEnd) : ClusterCoordinates(s,e,qs,qe,ts,te,st) {
		copy(gpBegin, gpEnd, back_inserter(matches));
		coarse=0;
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
		return end - start;
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
	void SetClusterBoundariesFromMatches(Options &opts) {
		for (int i=0; i < matches.size(); i++) {
			tEnd   = max(tEnd, matches[i].second.pos+opts.globalK);
			tStart = min(tStart, matches[i].second.pos );
			qEnd   = max(qEnd, matches[i].first.pos+opts.globalK);
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
	
  ClusterOrder(vector<Cluster> *c) : clusters(c) {
		index.resize(clusters->size());
		for (int i=0;i<index.size();i++) { index[i]=i;}
		Sort();
	}
		
	int operator()(const int i, const int j) {
			assert((*clusters)[i].strand == 0 or (*clusters)[i].strand == 1);
			assert((*clusters)[j].strand == 0 or (*clusters)[j].strand == 1);

		if ((*clusters)[i].strand != (*clusters)[j].strand) {
			return (*clusters)[i].strand < (*clusters)[j].strand;
		}
		else if ((*clusters)[i].tStart != (*clusters)[j].tStart) {
			return (*clusters)[i].tStart < (*clusters)[j].tStart;
		}
		else {
			return (*clusters)[i].qStart < (*clusters)[j].qStart;
		}
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

template<typename Tup>
void PrintDiagonal(vector<pair<Tup, Tup> > &matches) {
	for (int m=1; m < matches.size(); m++) {
		int64_t d=DiagonalDifference(matches[m], matches[m-1]);
		cout << matches[m-1].first.pos << "\t" << matches[m].first.pos << "\t" << matches[m-1].second.pos << "\t" << matches[m].second.pos << "\t" << d << endl;
	}
}
		

template<typename Tup>
void StoreDiagonalClusters(vector<pair<Tup, Tup> > &matches, vector<Cluster > &clusters, Options &opts, bool lite, int strand=0) {
	int i;
	int cs = 0, ce =0;
	cs = 0;
	int64_t dd,absdd;
	
	while (cs < matches.size()) {
		ce = cs+1;
		GenomePos qStart=matches[cs].first.pos, 
			qEnd=matches[cs].first.pos+opts.globalK, 
			tStart=matches[cs].second.pos, 
			tEnd=matches[cs].second.pos+opts.globalK;
		int diff=0;
		while (ce < matches.size() and 
					 abs(DiagonalDifference(matches[ce], matches[ce-1])) < opts.maxDiag  and 
					 (opts.maxGap == -1 or GapDifference(matches[ce], matches[ce-1]) < opts.maxGap) ) {
			qStart = min(qStart, matches[ce].first.pos);
			qEnd   = max(qEnd, matches[ce].first.pos+opts.globalK);
			tStart = min(tStart, matches[ce].second.pos);
			tEnd   = max(tEnd, matches[ce].second.pos+opts.globalK);
			ce++;
			if (ce < matches.size()) {
				diff = GapDifference(matches[ce], matches[ce-1]);
			}

		}			
		if (ce - cs >= opts.minClusterSize) {
			if (lite == false ) {
				clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand, matches.begin() + cs, matches.begin()+ce));
			}
			else {
				clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand));
			}
		}
		cs=ce;
	}
}

bool sign(int val) {
	if (val >= 0) return true;
	return false;
}

void RemovePairedIndels(GenomePos qAlnStart, GenomePos tAlnStart,
												GenomePos qAlnEnd,   GenomePos tAlnEnd,
												seqan::String<seqan::Seed<seqan::Simple> > &matches,
												 Options &opts) {
	int nMatches=seqan::length(matches);
	if ( nMatches < 3)   { return;}
	vector<bool> remove(nMatches, false);
	GenomePos prevQEnd, prevTEnd, qStart, tStart, qEnd, tEnd ;
	GenomePos nextQStart;
	GenomePos nextTStart;

	for (int c = 0; c < nMatches ; c++) {

		if (c == 0) {
			prevQEnd = qAlnStart;
			prevTEnd = qAlnEnd;
		}		
		else {
			prevQEnd = seqan::beginPositionV(matches[c-1]) + opts.globalK;
			prevTEnd = seqan::beginPositionH(matches[c-1])+opts.globalK;
		}
		qStart   = seqan::beginPositionV(matches[c]);
		tStart   = seqan::beginPositionH(matches[c]);
		qEnd     = seqan::beginPositionV(matches[c])+opts.globalK;
		tEnd     = seqan::beginPositionH(matches[c])+opts.globalK;

		if (c < nMatches-1) {
			nextQStart = seqan::beginPositionV(matches[c+1]);
			nextTStart = seqan::beginPositionH(matches[c+1]);
		}
		else {
			nextQStart = qAlnEnd;
			nextTStart = tAlnEnd;
		}


		int prevGap = (int)(qStart-prevQEnd) - (int)(tStart-prevTEnd);
		int nextGap = (int)(nextQStart - qEnd) - (int)(nextTStart - tEnd);

		if (sign(prevGap) != sign(nextGap) and
				abs(prevGap) + abs(nextGap) >  abs(prevGap + nextGap)  ) { //(matches[c].end +opts.k - matches[c].start)) {
			remove[c] = true;
		} 
	}	
	int m=0;

	for (int i=0; i < seqan::length(matches); i++) {
		if (remove[i] == false) {
			matches[m] = matches[i];
			m++;
		}
	}
	seqan::resize(matches, m);
}


template<typename Tup>
void RemovePairedIndels(vector<pair<Tup, Tup> > &matches, 
												vector<Cluster > &clusters, 
												Options &opts) {
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

		int prevGap = (int)(qStart-prevQEnd) - (int)(tStart-prevTEnd);
		int nextGap = (int)(nextQStart - qEnd) - (int)(nextTStart - tEnd);

		if (sign(prevGap) != sign(nextGap) and
				abs(prevGap) + abs(nextGap) >  abs(prevGap + nextGap)  ) { //(clusters[c].end +opts.k - clusters[c].start)) {
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


	
#endif
