#ifndef CLUSTERING_H_
#define CLUSTERING_H_
#include "Options.h"
#include <vector>
#include <climits>
#include <list>
#include <regex>
#include <string>
#include <unordered_map>
#include "Timing.h"

using namespace std;

class ClusterCoordinates {
 public:
	int start;
	int end;
	bool strand;
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
	vector<int> splitmatchindex; // stores the index for the split clusters
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
	bool split; // whether to enter the splitting step
	float anchorfreq;
	float refineEffiency;
	int matchStart;
	vector<bool> overlap;
	bool flip = 0;

	Cluster() { refined=0; coarse=-1; rank=-1; flip = 0;}
 Cluster(int s, int e) : ClusterCoordinates(s,e) { coarse=-1; refined=0; flip = 0;}

 Cluster(int s, int e, int st) : ClusterCoordinates(s,e,st) { coarse=-1; refined=0;flip = 0;}

  Cluster(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, 
					int st) : ClusterCoordinates(s,e,qs,qe,ts,te,st) { coarse=-1; refined=0;flip = 0;} 
  Cluster(int s, int e, 
					GenomePos qs, GenomePos qe,
					GenomePos ts, GenomePos te, 
					int st, int cs) : ClusterCoordinates(s,e,qs,qe,ts,te,st) { coarse=cs; refined=0;flip = 0;} 
	
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
		flip = 0;
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
		flip = 0;
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
		split = 0;
		flip = 0;
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

	void SetClusterBoundariesFromMatches (Options &opts, bool append_prev_cluster=0) {
		if (!append_prev_cluster) {
			qStart = matches[0].first.pos;
			qEnd = qStart + opts.globalK;
			tStart = matches[0].second.pos;
			tEnd = tStart + opts.globalK;
		}	
		for (int i = 1; i < matches.size(); i++) {
			tEnd = max(tEnd, matches[i].second.pos + opts.globalK);
			tStart = min(tStart, matches[i].second.pos);
			qEnd = max(qEnd, matches[i].first.pos + opts.globalK);
			qStart = min(qStart, matches[i].first.pos);
		}
		refineEffiency = ((float) matches.size()) / min(qEnd - qStart, tEnd - tStart);
	}	
	//
	// This function decide the chromIndex
	//
	bool CHROMIndex(Genome & genome) {
		if (matches.size() == 0) return 1;
		int firstChromIndex = genome.header.Find(tStart);
		int lastChromIndex;
		lastChromIndex = genome.header.Find(tEnd);
		if (firstChromIndex != lastChromIndex ) {
			return 1;
		}
		chromIndex = firstChromIndex;  
		return 0;
	}	
	GenomePos GetqStart(int i) {
		return matches[i].first.pos;
	}
	GenomePos GetqEnd(int i) {
		return matches[i].first.pos + matchesLengths[i];
	}
	GenomePos GettStart(int i) {
		return matches[i].second.pos;			
	}
	GenomePos GettEnd(int i) {
		return matches[i].second.pos + matchesLengths[i];
	}
};

class Cluster_SameDiag {
public:
	vector<int> start;
	vector<int> end;
	Cluster * cluster;
	GenomePos tStart, tEnd, qStart, qEnd;
	bool strand;
	int matchStart;
	int coarse;
	float anchorfreq;
	Cluster_SameDiag() {}
	Cluster_SameDiag(Cluster * c) : cluster(c) {
		tStart = cluster->tStart;
		tEnd = cluster->tEnd;
		qStart = cluster->qStart;
		qEnd = cluster->qEnd;
		strand = cluster->strand;
		matchStart = -1;
		coarse = -1;
		anchorfreq = cluster->anchorfreq;
	}
	int length(int i) {
		if (cluster->matches[end[i] - 1].first.pos + cluster->matchesLengths[end[i] - 1] >= cluster->matches[start[i]].first.pos) {
			return (cluster->matches[end[i] - 1].first.pos + cluster->matchesLengths[end[i] - 1] - cluster->matches[start[i]].first.pos);
		}
		else {return 0;}
	}
	GenomePos GetqStart(int i) {
		return cluster->matches[start[i]].first.pos;
	}
	GenomePos GetqEnd(int i) {
		return cluster->matches[end[i] - 1].first.pos + length(i);
	}
	GenomePos GettStart(int i) {
		if (strand == 0) return cluster->matches[start[i]].second.pos;			
		else return cluster->matches[end[i] - 1].second.pos;
	}
	GenomePos GettEnd(int i) {
		if (strand == 0) return cluster->matches[start[i]].second.pos + length(i);
		return cluster->matches[end[i] - 1].second.pos + length(i);
	}
	int size() {
		return start.size();
	}
	float OverlaprateOnGenome(const Cluster_SameDiag *b) const {
		int ovp=0;
		if (tEnd <= b->tStart or b->tEnd <= tStart) {
			return 0;
		}
		ovp = min(tEnd, b->tEnd) - max(tStart, b->tStart);
		float denomA = tEnd - tStart;
		float denomB = b->tEnd - b->tStart;
		return max(ovp/denomA, ovp/denomB);
	}
	int OverlapOnGenome(const Cluster_SameDiag *b) const {
		int ovp=0;
		if (tEnd <= b->tStart or b->tEnd <= tStart) {
			return 0;
		}
		ovp = min(tEnd, b->tEnd) - max(tStart, b->tStart);
		return ovp;
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

template<typename Tup>
long DiagonalDifference(Tup &a, Tup &b, int strand=0) {
	if (strand == 0) { // Matches
		long aDiag = (long)a.second.pos - (long)a.first.pos;
		long bDiag = (long)b.second.pos - (long)b.first.pos;
		return aDiag - bDiag;		
	}
	else { // revMathches
		long aDiag = a.first.pos + a.second.pos; 
		long bDiag= b.first.pos + b.second.pos;
		return aDiag - bDiag;				
	}
}

template<typename Tup>
int DiagonalDrift(long curDiag, Tup &t, int strand=0) {
	long drift;
	if (strand == 0) drift= abs(curDiag - ((long)t.second.pos - (long)t.first.pos));
	else drift= abs(curDiag - ((long)t.first.pos + (long)t.second.pos));
	return drift;
}

template<typename Tup>
long maxGapDifference(Tup &a, Tup &b) {
	long aDiff = abs((long)b.first.pos - (long)a.first.pos);
	long bDiff = abs((long)b.second.pos - (long)a.second.pos);
	return max(aDiff, bDiff);
}

template<typename Tup>
long minGapDifference (Tup &a, Tup &b) {
	long aDiff = abs((long)b.first.pos - (long)a.first.pos);
	long bDiff = abs((long)b.second.pos - (long)a.second.pos);
	return min(aDiff, bDiff);
}

template<typename Tup>
long GapDifference (Tup &a, Tup &b, int &len) {
	long Diff = abs((long)b.first.pos - (long)a.first.pos + len);
	return Diff;
}

bool sign(int val) {
	if (val >= 0) return true;
	return false;
}

template<typename Tup>
void AVGfreq(int &as, int &ae, vector<pair<Tup, Tup> > &matches, float &avgfreq) {
	unordered_map<Tuple, int> miniCount;
	for (int r = as; r < ae; r++) {
		unordered_map<Tuple, int>::const_iterator got = miniCount.find(matches[r].first.t);
		if (got == miniCount.end()) {
			miniCount[matches[r].first.t] = 0;
		}
	}
	// // cerr << "pr.second - pr.first: " << pr.second - pr.first << " miniCount.size(): " << miniCount.size() << endl;
	avgfreq = (float)(ae - as) / miniCount.size();
}

// template<typename Tup>
// void CleanOffDiagonal(Genome &genome, vector<float> &matches_freq, vector<pair<Tup, Tup> > &matches, Options &opts, Read &read, int strand=0, int diagOrigin=-1, int diagDrift=-1) {
// 	if (matches.size() == 0) {
// 		return;
// 	}

// 	vector<bool> onDiag(matches.size(), false);
// 	int nOnDiag=0;
// 	if (matches.size() > 1 and abs(DiagonalDifference(matches[0], matches[1], strand)) < opts.cleanMaxDiag 
// 						   and (diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[0], strand) < diagDrift)) { 
// 		onDiag[0] = true;
// 		nOnDiag++;
// 	}
	
// 	for (int i = 1; i < matches.size(); i++) {
// 		if (abs(DiagonalDifference(matches[i], matches[i-1], strand)) < opts.cleanMaxDiag 
// 			and (diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[i], strand) < diagDrift)) {	
// 			onDiag[i] = true;
// 			onDiag[i-1] = true;
// 			nOnDiag++;
// 		}
// 	}
// 	bool prevOnDiag = false;
// 	int  diagStart;
// 	int  Largest_ClusterNum = 0;
// 	int ct = 0;
// 	for (int i = 0; i < matches.size(); i++) {
// 		if (prevOnDiag == false and onDiag[i] == true) {
// 			diagStart = i;
// 		}
// 		if (prevOnDiag == true and onDiag[i] == false) {
// 			Largest_ClusterNum = max(Largest_ClusterNum, i - diagStart);
// 			// if (opts.dotPlot and strand == 1){
// 			// 	ofstream rclust("rev-matches_cleanoffdiagonal.dots", ofstream::app);
// 			// 	for (int j = diagStart; j < i; j++) {			
// 			// 		rclust << matches[j].first.pos << "\t" << matches[j].second.pos + opts.globalK << "\t" << opts.globalK + matches[j].first.pos << "\t"
// 			// 				 << matches[j].second.pos << "\t" << ct << "\t" << 0 << endl;
// 			// 	}
// 			// 	rclust.close();
// 			// }
// 		}
// 		prevOnDiag = onDiag[i];
// 		ct++;
// 	}
// 	int sz = matches.size();
// 	Largest_ClusterNum = max(Largest_ClusterNum, sz - diagStart);

// 	int minDiagCluster = (int) floor(Largest_ClusterNum/10);
// 	if (minDiagCluster >= opts.minDiagCluster) minDiagCluster = opts.minDiagCluster;
// 	// cerr << "Largest_ClusterNum: " << Largest_ClusterNum << " minDiagCluster: " << minDiagCluster << endl;
// 	// if (opts.readType == Options::contig) {minDiagCluster = 30;}

// 	// 
// 	// Remove bins with number of anchors <= minDiagCluster
// 	//
// 	int counter = 0;
// 	if (minDiagCluster > 0) {
// 		for (int i = 0; i < matches.size(); i++) {
// 			if (prevOnDiag == false and onDiag[i] == true) {
// 				diagStart = i;
// 			}
// 			if (prevOnDiag == true and onDiag[i] == false) {
// 				if (i - diagStart < minDiagCluster) { // used to be minDiagCluster not opts.minDiagCluster
// 					for (int j = diagStart; j < i; j++) {
// 						onDiag[j] = false;
// 					}					
// 				}
// 				else {
// 					float avgfreq;
// 					AVGfreq(diagStart, i, matches, avgfreq);

// 					if (opts.dotPlot and opts.readname == read.name and strand == 0) {
// 						ofstream fclust("for-matches_cleanoffdiagonal.dots", ofstream::app);
// 						for (int j = diagStart; j < i; j++) {
// 							fclust << matches[j].first.pos << "\t" << matches[j].second.pos << "\t" << opts.globalK + matches[j].first.pos << "\t"
// 									<< matches[j].second.pos + opts.globalK << "\t" << counter << "\t" << avgfreq << endl;
// 						}
// 						fclust.close();
// 					}
// 					if (opts.dotPlot and opts.readname == read.name and strand == 1){
// 						ofstream rclust("rev-matches_cleanoffdiagonal.dots", ofstream::app);
// 						for (int j = diagStart; j < i; j++) {			
// 							rclust << matches[j].first.pos << "\t" << matches[j].second.pos + opts.globalK << "\t" << opts.globalK + matches[j].first.pos << "\t"
// 									 << matches[j].second.pos << "\t" << counter << "\t" << avgfreq << endl;
// 						}
// 						rclust.close();
// 					}
// 					for (int j = diagStart; j < i; j++) {matches_freq[j] = avgfreq;}
// 					int MinDiagCluster = 0;	
// 					if (avgfreq >= 2.0f and i - diagStart >= opts.cleanClustersize) {
// 						MinDiagCluster = opts.SecondCleanMinDiagCluster + floor((avgfreq - 1.5f)/ 1.0f) * opts.punish_anchorfreq + 
// 								floor((i - diagStart - opts.cleanClustersize)/opts.cleanClustersize) * opts.anchorPerlength;
// 						SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);						
// 					}
// 					else if (avgfreq >= 1.5f and i - diagStart >= opts.cleanClustersize) {
// 						MinDiagCluster = opts.SecondCleanMinDiagCluster + floor((avgfreq - 1.5f)/ 1.5f) * opts.punish_anchorfreq + 
// 								floor((i - diagStart - opts.cleanClustersize)/opts.cleanClustersize) * opts.anchorPerlength;
// 						SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);						
// 					}
// 					// cerr << "avgfreq: " << avgfreq << " MinDiagCluster: " << MinDiagCluster << " size: " << i - diagStart << endl;

// 					// else if (avgfreq > 1.0f) {
// 					// 	MinDiagCluster = opts.SecondCleanMinDiagCluster - (5 - floor((avgfreq - 1.0f) / 0.1f)) * (opts.punish_anchorfreq/2) - 
// 					// 			floor((opts.cleanClustersize - i + diagStart)/15) * (opts.anchorPerlength / 2);
// 					// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
// 					// }
// 					counter++;					
// 				}
// 			}
// 			prevOnDiag = onDiag[i];
// 		}		
// 	}
// 	else {
// 		for (int i = 0; i < matches.size(); i++) {onDiag[i] = false;}
// 	}

// 	int c = 0;
// 	for (int i = 0; i < matches.size(); i++) {
// 		if (onDiag[i]) {
// 			matches[c] = matches[i]; 
// 			matches_freq[c] = matches_freq[i];
// 			c++;
// 		}
// 	}
// 	matches.resize(c);
// 	matches_freq.resize(c);

// 	// prevOnDiag = false;
// 	// for (int i = 0; i < matches.size(); i++) {
// 	// 	if (prevOnDiag == false and onDiag[i] == true) {diagStart = i;}
// 	// 	if (prevOnDiag == true and onDiag[i] == false) {
// 	// 		clusters.push_back(Cluster(0, 0, strand)); 
// 	// 		for (int j = diagStart; j < i; j++) {
// 	// 			clusters.back().matches.push_back(matches[j]);
// 	// 		}		
// 	// 		clusters.back().SetClusterBoundariesFromMatches(opts);
// 	// 		int rci = genome.header.Find(clusters.back().tStart);
// 	// 		clusters.back().chromIndex = rci;
// 	// 	}
// 	// 	prevOnDiag = onDiag[i];
// 	// }
// }

template<typename Tup>
void CleanOffDiagonal(vector<pair<Tup, Tup> > &matches, vector<float> &matches_freq, Options &opts, Read &read, int strand=0, int diagOrigin=-1, int diagDrift=-1) {
	if (matches.size() == 0) {
		return;
	}

	vector<bool> onDiag(matches.size(), false);
	int nOnDiag=0;
	if (matches.size() > 1 and abs(DiagonalDifference(matches[0], matches[1], strand)) < opts.cleanMaxDiag 
						   and (diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[0], strand) < diagDrift)) { 
		onDiag[0] = true;
		nOnDiag++;
	}
	
	for (int i = 1; i < matches.size(); i++) {
		if (abs(DiagonalDifference(matches[i], matches[i-1], strand)) < opts.cleanMaxDiag 
			and (diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[i],strand) < diagDrift)) {	
			onDiag[i] = true;
			onDiag[i-1] = true;
			nOnDiag++;
		}
	}
	bool prevOnDiag = false;
	int  diagStart;
	int  Largest_ClusterNum = 0;
	for (int i = 0; i < matches.size(); i++) {
		if (prevOnDiag == false and onDiag[i] == true) {
			diagStart = i;
		}
		if (prevOnDiag == true and onDiag[i] == false) {
			Largest_ClusterNum = max(Largest_ClusterNum, i - diagStart);
		}
		prevOnDiag = onDiag[i];
	}
	int sz = matches.size();
	Largest_ClusterNum = max(Largest_ClusterNum, sz - diagStart);
	int minDiagCluster = (int) floor(Largest_ClusterNum/10);
	if (minDiagCluster >= opts.minDiagCluster) minDiagCluster = opts.minDiagCluster;
	// cerr << "Largest_ClusterNum: " << Largest_ClusterNum << " minDiagCluster: " << minDiagCluster << endl;
	// if (opts.readType == Options::contig) {minDiagCluster = 30;}
	// 
	// Remove bins with number of anchors <= minDiagCluster
	//
	int counter = 0;
	if (minDiagCluster > 0) {
		for (int i = 0; i < matches.size(); i++) {
			if (prevOnDiag == false and onDiag[i] == true) {
				diagStart = i;
			}
			if (prevOnDiag == true and onDiag[i] == false) {
				if (i - diagStart < minDiagCluster) { // used to be minDiagCluster not opts.minDiagCluster
					for (int j = diagStart; j < i; j++) {
						onDiag[j] = false;
					}
					// if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname and strand == 0) {
					// 	ofstream fclust("for-matches_1.dots", ofstream::app);
					// 	for (int j = diagStart; j < i; j++) {
					// 		fclust << matches[j].first.pos << "\t" << matches[j].second.pos << "\t" << opts.globalK + matches[j].first.pos << "\t"
					// 				<< matches[j].second.pos + opts.globalK << "\t" << counter << "\t" << 0 << endl;
					// 	}
					// 	fclust.close();
					// }
					// if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname and strand == 1){
					// 	ofstream rclust("rev-matches_1.dots", ofstream::app);
					// 	for (int j = diagStart; j < i; j++) {			
					// 		rclust << matches[j].first.pos << "\t" << matches[j].second.pos + opts.globalK << "\t" << opts.globalK + matches[j].first.pos << "\t"
					// 				 << matches[j].second.pos <<"\t" << counter << "\t" << 0 << endl;
					// 	}
					// 	rclust.close();
					// }
				}
				else {
					float avgfreq;
					AVGfreq(diagStart, i, matches, avgfreq);
					if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname and strand == 0) {
						ofstream fclust("for-matches_cleanoffdiagonal.dots", ofstream::app);
						for (int j = diagStart; j < i; j++) {
							fclust << matches[j].first.pos << "\t" << matches[j].second.pos << "\t" << opts.globalK + matches[j].first.pos << "\t"
									<< matches[j].second.pos + opts.globalK << "\t" << counter << "\t" << avgfreq << endl;
						}
						fclust.close();
					}
					if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname and strand == 1){
						ofstream rclust("rev-matches_cleanoffdiagonal.dots", ofstream::app);
						for (int j = diagStart; j < i; j++) {			
							rclust << matches[j].first.pos << "\t" << matches[j].second.pos + opts.globalK << "\t" << opts.globalK + matches[j].first.pos << "\t"
									 << matches[j].second.pos << "\t" << counter << "\t" << avgfreq << endl;
						}
						rclust.close();
					}
					for (int j = diagStart; j < i; j++) {
						matches_freq[j] = avgfreq;
					}
					int MinDiagCluster = 0;

					if (opts.bypassClustering) {
						if (avgfreq >= 2.0f and i - diagStart >= opts.cleanClustersize) {
							MinDiagCluster = opts.SecondCleanMinDiagCluster + floor((avgfreq - 1.5f)/ 1.0f) * opts.punish_anchorfreq + 
									floor((i - diagStart - opts.cleanClustersize)/opts.cleanClustersize) * opts.anchorPerlength;
							SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);						
						}
						else if (avgfreq >= 1.5f and i - diagStart >= opts.cleanClustersize) {
							MinDiagCluster = opts.SecondCleanMinDiagCluster + floor((avgfreq - 1.5f)/ 1.5f) * opts.punish_anchorfreq + 
									floor((i - diagStart - opts.cleanClustersize)/opts.cleanClustersize) * opts.anchorPerlength;
							SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);	
						}	
					}
					else {
						// if (avgfreq >= 40.0f) {
						// 	MinDiagCluster = 1000;
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);						
						// }
						// else if (avgfreq >= 30.0f) {
						// 	MinDiagCluster = 500;
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);						
						// }		
						if (avgfreq >= 4.0f and i - diagStart >= opts.cleanClustersize) {
							MinDiagCluster = opts.SecondCleanMinDiagCluster + floor((avgfreq - 1.5f)/ 1.0f) * opts.punish_anchorfreq + 
									floor((i - diagStart - opts.cleanClustersize)/opts.cleanClustersize) * opts.anchorPerlength;
							SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);						
						}
						else if (avgfreq >= 1.5f and i - diagStart >= opts.cleanClustersize) {
							MinDiagCluster = opts.SecondCleanMinDiagCluster + floor((avgfreq - 1.5f)/ 1.5f) * opts.punish_anchorfreq + 
									floor((i - diagStart - opts.cleanClustersize)/opts.cleanClustersize) * opts.anchorPerlength;
							SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);						
						}
						// else if (avgfreq >= 15.0f) {
						// 	int MinDiagCluster = 100;
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);						
						// }
						// else if (avgfreq >= 5.0f) {
						// 	int MinDiagCluster = opts.SecondCleanMinDiagCluster + 30;
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						// }
						// else if (avgfreq >= 4.0f) {
						// 	int MinDiagCluster = opts.SecondCleanMinDiagCluster + 20;
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						// }	
						// else if (avgfreq >= 3.0f) {
						// 	int MinDiagCluster = opts.SecondCleanMinDiagCluster + 10; //30
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						// }
						// else if (avgfreq >= 1.5f and i - diagStart >= 100) {
						// 	int MinDiagCluster = opts.SecondCleanMinDiagCluster; //30
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						// }
						// else if (avgfreq >= 1.4f and i - diagStart >= 100) {
						// 	int MinDiagCluster = opts.SecondCleanMinDiagCluster - 5; //30
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						// }
						else if (avgfreq > 1.0f and i - diagStart >= opts.cleanClustersize) {
							MinDiagCluster = opts.SecondCleanMinDiagCluster - (5 - floor((avgfreq - 1.0f) / 0.1f)) * (opts.punish_anchorfreq/2) + 
									floor((i - diagStart - opts.cleanClustersize)/opts.cleanClustersize) * (opts.anchorPerlength / 2);
							SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						}
						else if (avgfreq > 1.0f) {
							MinDiagCluster = opts.SecondCleanMinDiagCluster - (5 - floor((avgfreq - 1.0f) / 0.1f)) * (opts.punish_anchorfreq/2) - 
									floor((opts.cleanClustersize - i + diagStart)/15) * (opts.anchorPerlength / 2);
							SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						}
						// else if (avgfreq >= 1.3f and i - diagStart >= 100) {
						// 	int MinDiagCluster = opts.SecondCleanMinDiagCluster - 5; //30
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						// }
						// else if (avgfreq >= 1.2f and i - diagStart >= 100) {
						// 	int MinDiagCluster = opts.SecondCleanMinDiagCluster - 15; //30
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						// }
						// else if (avgfreq >= 1.1f and i - diagStart >= 100) {
						// 	int MinDiagCluster = opts.SecondCleanMinDiagCluster - 10; //30
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						// }
						// else if (avgfreq > 1.0f) {
						// 	int MinDiagCluster = 10; //30
						// 	SecondRoundCleanOffDiagonal(matches, MinDiagCluster, opts.SecondCleanMaxDiag, onDiag, diagStart, i, strand, diagOrigin, diagDrift);
						// }						
					}

					if (read.name == opts.readname) cerr << "avgfreq: " << avgfreq << " MinDiagCluster: " << MinDiagCluster << endl;
					counter++;
				}
			}
			prevOnDiag = onDiag[i];
		}		
	}
	else {
		for (int i = 0; i < matches.size(); i++) {onDiag[i] = false;}
	}

	int c = 0;
	for (int i = 0; i < matches.size(); i++) {
		if (onDiag[i]) {
			matches[c] = matches[i]; 
			matches_freq[c] = matches_freq[i];
			c++;
		}
	}
	matches.resize(c);
	matches_freq.resize(c);

	// if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname and strand == 0) {
	// 	ofstream fclust("for-matches_clean.dots");
	// 	for (int m = 0; m < matches.size(); m++) {
	// 		fclust << matches[m].first.pos << "\t" << matches[m].second.pos << "\t" << opts.globalK + matches[m].first.pos << "\t"
	// 				<<  matches[m].second.pos + opts.globalK << "\t" << m << "\t" <<endl;
	// 	}
	// 	fclust.close();
	// }
	// if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname and strand == 1){
	// 	ofstream rclust("rev-matches_clean.dots");
	// 	for (int m=0; m < matches.size(); m++) {			
	// 		rclust <<  matches[m].first.pos << "\t" <<  matches[m].second.pos + opts.globalK << "\t" << opts.globalK + matches[m].first.pos  << "\t"
	// 				 <<  matches[m].second.pos << m << "\t" << endl;
	// 	}
	// 	rclust.close();
	// }
}

template<typename Tup>
void SecondRoundCleanOffDiagonal(vector<pair<Tup, Tup> > &matches, int &MinDiagCluster,int &CleanMaxDiag, vector<bool> &OriginalOnDiag, int &os, int &oe, 
									int strand=0, int diagOrigin=-1, int diagDrift=-1) {
	if (MinDiagCluster < 0 or MinDiagCluster >= oe - os) return;
	int nOnDiag2=0;
	for (int i = os; i < oe; i++) {
		OriginalOnDiag[i] = false;
	}	
	if ((oe - os) > 1 and abs(DiagonalDifference(matches[os], matches[os + 1], strand)) < CleanMaxDiag 
						   and (diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[os], strand) < diagDrift)) { 
		OriginalOnDiag[os] = true;
		nOnDiag2++;
	}
	
	for (int i = os + 1; i < oe; i++) {
		if (abs(DiagonalDifference(matches[i], matches[i - 1], strand)) < CleanMaxDiag 
			and (diagOrigin == -1 or DiagonalDrift(diagOrigin, matches[i],strand) < diagDrift)) {	
			OriginalOnDiag[i] = true;
			OriginalOnDiag[i- 1] = true;
			nOnDiag2++;
		}
	}

	bool prevOnDiag = false;
	int  diagStart;
	for (int i = os; i < oe; i++) {
		if (prevOnDiag == false and OriginalOnDiag[i] == true) {
			diagStart = i;
		}
		if (prevOnDiag == true and OriginalOnDiag[i] == false) {
			if (i - diagStart < MinDiagCluster) { // used to be minDiagCluster not opts.minDiagCluster
				for (int j = diagStart; j < i; j++) {
					OriginalOnDiag[j] = false;
				}
			}
		}
		prevOnDiag = OriginalOnDiag[i];
	}
}

// class LogCluster {
//  public:
//  	vector<Cluster> SubCluster;
//  	Cluster * Hp;
//  	bool ISsecondary; // ISsecondary == 1 means this is a secondary chain. Otherwise it's a primary chain
//  	int primary; // When ISsecondary == 1, primary stores the index of the primary chain in vector<LogCluster>
//  	vector<int> secondary; // When ISsecondary == 0, secondary stores the indices of the secondary chains
//  	int coarse; // coarse means this LogCluster stores information about refinedCluster[coarse]
//  	bool direction; // direction means how the read is mapped to the reference
//  	bool split; // split == 1 means this chain has been split
//  	int main; // If a read is splitted, "main" stores the index of the main alignment part

//  	LogCluster () {
//  		ISsecondary = 0;
//  		primary = -1;
//  		coarse = -1;
//  		direction = 0;
//  		split = 0;
//  	};
//  	~LogCluster() {};
 	
//  	void setHp(Cluster & H) {
//  		Hp = & H;
//  	}
//  	void SetCoarse () {
//  		coarse = SubCluster[0].coarse;
//  	}

//  	void SetSubClusterBoundariesFromMatches (Options &opts) {
// 		// set the boundaries for SubCluster[i] -- the current last one in SubCluster
//  		int i = SubCluster.size() - 1; // the last one in SubCluster
// 		for (int is = SubCluster[i].start; is < SubCluster[i].end; ++is) {

// 			if (is == SubCluster[i].start) {
// 		 		SubCluster[i].tStart = Hp->matches[is].second.pos;
// 		 		SubCluster[i].qStart = Hp->matches[is].first.pos;					
// 			}
// 			SubCluster[i].tEnd   = max(SubCluster[i].tEnd, Hp->matches[is].second.pos + opts.globalK);
// 			SubCluster[i].tStart = min(SubCluster[i].tStart, Hp->matches[is].second.pos);
// 			SubCluster[i].qEnd   = max(SubCluster[i].qEnd, Hp->matches[is].first.pos + opts.globalK);
// 			SubCluster[i].qStart = min(SubCluster[i].qStart, Hp->matches[is].first.pos); 	
// 		}
//  	}
// };


template<typename Tup>
void PrintDiagonal(vector<pair<Tup, Tup> > &matches, int strand=0) {
	for (int m=1; m < matches.size(); m++) {
		long d=DiagonalDifference(matches[m], matches[m-1], strand);
		cerr << matches[m-1].first.pos << "\t" << matches[m].first.pos << "\t" << matches[m-1].second.pos << "\t" << matches[m].second.pos << "\t" << d << endl;
	}
}
 
template<typename Tup>
long GetDiag(pair<Tup, Tup> &match, int strand, Options &opts) {
	if (strand == 0) return (long) match.second.pos - (long) ceil(opts.slope*match.first.pos);
	else return (long) ceil(opts.slope*match.first.pos) + (long) match.second.pos;
}

template<typename Tup>
long GetDiag(pair<Tup, Tup> &match, int &len, bool strand) {
	if (strand == 0) return (long) match.second.pos - (long) match.first.pos;
	else return (long) match.first.pos + (long) match.second.pos + len;
}

template<typename Tup>
void StoreFineClusters(int ri, vector<pair<Tup, Tup> > &matches, float anchorfreq, vector<Cluster> &clusters, Options &opts, vector<int> &splitmatchindex,
											 Genome &genome, Read &read, GenomePos readLength, int strand=0, int outerIteration=0) {
	//
	// StoreFineClusters based on unique part
	// The first step: store the number of matches corresponding to each minimizer from the read; -- a linear pass
	// matches are sorted in cartesian order; (first, second)
	//
	if (splitmatchindex.size() == 1) return;
	if (anchorfreq == 1.0f) {
		clusters.push_back(Cluster(0, 0, strand)); 
		for (int i = 0; i < splitmatchindex.size(); i++) {
			clusters.back().matches.push_back(matches[splitmatchindex[i]]);
		}		
		clusters.back().SetClusterBoundariesFromMatches(opts);
		clusters.back().chromIndex = ri;
		clusters.back().anchorfreq = 1.0f;
		// cerr << "clusters.back().anchorfreq: " << clusters.back().anchorfreq << endl;
		if (clusters.back().CHROMIndex(genome) == 1) {
			clusters.pop_back();
		}		
		if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
		return;
	}

	vector<int> match_num;
	vector<int> pos_start;
	int oc = 1; 
	int us = 0;
	for (int i = 1; i < splitmatchindex.size(); i++) {
		if (matches[splitmatchindex[i]].first.pos ==  matches[splitmatchindex[i - 1]].first.pos) {
			oc++;
		}
		else {
			match_num.push_back(oc);
			pos_start.push_back(us);
			us = i;
			oc = 1;
		}
		if (i == splitmatchindex.size() - 1) {
			match_num.push_back(oc);
			pos_start.push_back(us);				
		}
	}
	//
	// Find stretches of "1" on a diagonal in match_num; -- a linear pass
	//
	int u_start = 0, u_end = 0, u_maxstart = 0, u_maxend = 0;
	int max_pos = 0;
	bool reset = 0; // reset == 0 means a new stretch of "1";
	vector<int> Start;
	vector<int> End;
	if (match_num.size() == 1) {
		assert(pos_start[0] == 0);
		u_maxstart = 0; 
		u_maxend = 1;	
		Start.push_back(u_maxstart);
		End.push_back(u_maxend);		
	}
	else {
		int k = 0;
		while (k < match_num.size() - 1) {
			// Find a start of a stretch
			while (k < match_num.size() - 1 and match_num[k] != 1) {
				// (Check) if reset == 1, push u_start and u_end. reset = 0
				k++;
			}
			assert(reset == 0);
			u_start = k;
			u_end = k + 1;	
			reset = 1;				
			while (k < match_num.size() - 1 and match_num[k+1] == match_num[k] 
				and abs(DiagonalDifference(matches[splitmatchindex[pos_start[k+1]]], matches[splitmatchindex[pos_start[k]]], strand)) < opts.maxDiag
				and minGapDifference(matches[splitmatchindex[pos_start[k+1]]], matches[splitmatchindex[pos_start[k]]]) <= opts.maxGap) { 
				u_end = k + 2;
				k++;
				// reset = 1;
			}
			// insert a new stretch of unique minimizers
			Start.push_back(u_start);
			End.push_back(u_end);
			reset = 0;
			k++;
			if ((u_maxstart == 0 and u_maxend == 0) or (u_maxend - u_maxstart < u_end - u_start)) {
				u_maxstart = u_start;
				u_maxend = u_end;		
				max_pos = Start.size() - 1; 
			}		
		}	
	}
	//
	// Decide the u_minDiag and u_maxDiag of the unique stretches of anchors
	//
	if (u_maxstart == 0 and u_maxend == 0) { // no unique stretch
		return;
	}
	int c_s = pos_start[u_maxstart], c_e = pos_start[u_maxend - 1] + 1;
	assert(c_s >= 0 and c_e <= splitmatchindex.size() and c_s <= c_e);

	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
		ofstream uniquematch("UniqueMatch.tab", ofstream::app);
		for (int h = 0; h < Start.size(); h++) {
			for (int he = Start[h]; he < End[h]; he++) {
				if (strand == 0) {
					uniquematch << matches[splitmatchindex[pos_start[he]]].first.pos << "\t"
						  << matches[splitmatchindex[pos_start[he]]].second.pos << "\t"
						  << matches[splitmatchindex[pos_start[he]]].first.pos  + opts.globalK << "\t"
						  << matches[splitmatchindex[pos_start[he]]].second.pos + opts.globalK << "\t"
						  << outerIteration << "\t"
						  << h << "\t"
						  << End[h] - Start[h] << "\t"
						  << strand << endl;
				}
				else {
					uniquematch << matches[splitmatchindex[pos_start[he]]].first.pos << "\t"
						  << matches[splitmatchindex[pos_start[he]]].second.pos + opts.globalK << "\t"
						  << matches[splitmatchindex[pos_start[he]]].first.pos + opts.globalK << "\t"
						  << matches[splitmatchindex[pos_start[he]]].second.pos << "\t"
						  << outerIteration << "\t"
						  << h << "\t"		
						  << End[h] - Start[h] << "\t"		  
						  << strand << endl;					
				}					
			}
		}
		uniquematch.close();
	}	

	if (c_e - c_s >= opts.minUniqueStretchNum 
		and matches[splitmatchindex[c_e-1]].first.pos + opts.globalK - matches[splitmatchindex[c_s]].first.pos >= opts.minUniqueStretchDist) { 
		clusters.push_back(Cluster(0, 0, strand)); 
		vector<bool> AddOrNot(Start.size(), 0);
		list<int> StretchOfOne;
		bool append_prev_cluster = 0;
		int clast = clusters.size() - 2;
		if (c_e - c_s == splitmatchindex.size()) { // no need to extend - the whole rough cluster is a unique linear part
			for (int i = c_s; i < c_e; i++) {
				clusters.back().matches.push_back(matches[splitmatchindex[i]]);
			}	
			clusters.back().anchorfreq = anchorfreq;
			// cerr << "clusters.back().anchorfreq: " << clusters.back().anchorfreq << endl;
			AddOrNot[0] = 1;
		}
		else {
			//
			// extend the unique stretches from the left; -- a linear pass in vector<int> Start and End;
			//
			int prev_anchor = c_s;
			int prev_stretch = max_pos;
			int addup = 0;
			if (max_pos >= 0) {
				StretchOfOne.push_back(max_pos);
				AddOrNot[max_pos] = 1; // Add the stretch
				int i = max_pos - 1;
				while (i >= 0) {
					int i_m = pos_start[End[i]-1];
					// cerr << "Diag: " << prev_stretch << " and " << i << " " << abs(DiagonalDifference(matches[i_m], matches[prev_anchor], strand)) << endl;
					// cerr << "matches[i_m].first.pos: " << matches[i_m].first.pos << " matches[i_m].second.pos: " << matches[i_m].second.pos << endl;
					// cerr << "matches[prev_anchor].first.pos: " << matches[prev_anchor].first.pos << " matches[prev_anchor].second.pos: " << matches[prev_anchor].second.pos << endl;
					// cerr << "prev_anchor: " << prev_anchor << endl;
					if ((abs(DiagonalDifference(matches[splitmatchindex[i_m]], matches[splitmatchindex[prev_anchor]], strand)) <= opts.maxDiag 
							and minGapDifference(matches[splitmatchindex[i_m]], matches[splitmatchindex[prev_anchor]]) <= opts.maxGap)
									or minGapDifference(matches[splitmatchindex[i_m]], matches[splitmatchindex[prev_anchor]]) <= opts.maxGap/2) {
						assert(i < StretchOfOne.back());
						StretchOfOne.push_back(i);
						AddOrNot[i] = 1; // Add the stretch
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
					if ((abs(DiagonalDifference(matches[splitmatchindex[i_m]], matches[splitmatchindex[prev_anchor]], strand)) <= opts.maxDiag 
							and minGapDifference(matches[splitmatchindex[i_m]], matches[splitmatchindex[prev_anchor]]) <= opts.maxGap)
										or minGapDifference(matches[splitmatchindex[i_m]], matches[splitmatchindex[prev_anchor]]) <= opts.maxGap/2) {
						assert(i > StretchOfOne.front());
						StretchOfOne.push_front(i);
						AddOrNot[i] = 1; // Add the stretch
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
			for (list<int>::reverse_iterator it = StretchOfOne.rbegin(); it != StretchOfOne.rend(); ++it) {
				vector<int> Cluster_index;
				c_s = pos_start[Start[*it]]; c_e = pos_start[End[*it]-1]+1;

				if (it == StretchOfOne.rbegin()) {
					if (*it == 0) p_s = 0;
					else p_s = pos_start[End[*it - 1]];
					p_e = pos_start[Start[*it]];
				}
				else {
					//p_s = pos_start[End[prev_stretch]-1]; ///??? -1
					p_s = pos_start[End[prev_stretch]]; ///??? -1
					p_e = pos_start[Start[*it]];
				}
				prev_stretch = *it;
				assert(p_s >= 0 and p_e <= splitmatchindex.size() and p_s <= p_e and p_e == c_s);

				int prev_match = c_s;
				for (int si = p_e - 1; si >= p_s; si--) {
					if (abs(DiagonalDifference(matches[splitmatchindex[si]], matches[splitmatchindex[prev_match]], strand)) < opts.maxDiag) { //opts.maxDiag-200
						Cluster_index.push_back(si);
						prev_match = si;
					}
				}
				for (vector<int>::reverse_iterator ci = Cluster_index.rbegin(); ci != Cluster_index.rend(); ++ci) {
					clusters.back().matches.push_back(matches[splitmatchindex[*ci]]);
					// for debug
					if (opts.debug and clusters.back().matches.size() >= 2) {
						int mb = clusters.back().matches.size() - 1;
						assert(clusters.back().matches[mb].first.pos >= clusters.back().matches[mb-1].first.pos);
						if (clusters.back().matches[mb].first.pos == clusters.back().matches[mb-1].first.pos) {
							assert(clusters.back().matches[mb].second.pos >= clusters.back().matches[mb-1].second.pos);
						}
					}
					
				}
				// insert the unique anchors
				for (int si = c_s; si < c_e; si++) {
					clusters.back().matches.push_back(matches[splitmatchindex[si]]);
				}
				if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
					ofstream uniquepart("UniquePart.tab", std::ofstream::app);
					for (int si = c_s; si < c_e; si++) {
						if (strand == 0) {
							uniquepart << matches[splitmatchindex[si]].first.pos << "\t"
								  << matches[splitmatchindex[si]].second.pos << "\t"
								  << matches[splitmatchindex[si]].first.pos+ opts.globalK << "\t"
								  << matches[splitmatchindex[si]].second.pos + opts.globalK << "\t"
								  << outerIteration << "\t"
								  << strand << endl;
						}
						else {
							uniquepart << matches[splitmatchindex[si]].first.pos << "\t"
								  << matches[splitmatchindex[si]].second.pos + opts.globalK << "\t"
								  << matches[splitmatchindex[si]].first.pos + opts.globalK << "\t"
								  << matches[splitmatchindex[si]].second.pos << "\t"
								  << outerIteration << "\t"
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
					if (*it == AddOrNot.size() - 1) {p_e = splitmatchindex.size();}
					else { p_e = pos_start[Start[*it + 1]];}
					assert(p_s >= 0 and p_e <= splitmatchindex.size() and p_s <= p_e);
					prev_match = c_e - 1;
					for (int si = p_s; si < p_e; si++) {
						if (abs(DiagonalDifference(matches[splitmatchindex[si]], matches[splitmatchindex[prev_match]], strand)) < opts.maxDiag) {//opts.maxDiag-200
							clusters.back().matches.push_back(matches[splitmatchindex[si]]);
							prev_match = si;
							// debug
							if (opts.debug and clusters.back().matches.size() >= 2) {
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
			clusters.back().anchorfreq = anchorfreq;
			// cerr << "clusters.back().anchorfreq: " << clusters.back().anchorfreq << endl;
		}
		//
		// Decide the cooridinates of the cluster;
		//
		clusters.back().SetClusterBoundariesFromMatches(opts, append_prev_cluster);
		clusters.back().chromIndex = ri;
		// // debug
		// if (clusters.size() >= 2) {
		// 	assert(!(clusters.back().qStart == clusters[clusters.size() - 2].qStart and clusters.back().tStart == clusters[clusters.size() - 2].tStart
		// 		and clusters.back().qEnd == clusters[clusters.size() - 2].qEnd and clusters.back().tEnd == clusters[clusters.size() - 2].tEnd)
		// 		and clusters.back().matches[0].first.t != clusters[clusters.size() - 2].matches[0].first.t);
		// }

		if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
		
		//
		// remove Cluster that firstChromIndex != lastChromIndex; or remove Cluster of only a few matches;
		//
		if (clusters.back().matches.size() <= opts.minClusterSize and clusters.back().CHROMIndex(genome) == 1) {
			clusters.pop_back();
		}
		//
		// remove Cluster that has too high slope
		//
		if (((long)clusters.back().tEnd - (long)clusters.back().tStart) / ((long) clusters.back().qEnd - (long) clusters.back().qStart) >= 5) {
			clusters.pop_back();
		}
		
		//Check the rest unadded stretches. Add them to clusters
		
		for (int ar = 0; ar < AddOrNot.size(); ar++) { 
			if (!AddOrNot[ar] and End[ar] - Start[ar] >= 15) {
				clusters.push_back(Cluster(0, 0, strand)); 
				for (int i = pos_start[Start[ar]]; i < pos_start[End[ar] - 1] + 1; i++) {
					clusters.back().matches.push_back(matches[splitmatchindex[i]]);
				}		
				clusters.back().SetClusterBoundariesFromMatches(opts);
				clusters.back().chromIndex = ri;
				clusters.back().anchorfreq = anchorfreq;
				// cerr << "clusters.back().anchorfreq: " << clusters.back().anchorfreq << endl;
				if (clusters.back().CHROMIndex(genome) == 1) {
					clusters.pop_back();
				}
				if (((long)clusters.back().tEnd - (long)clusters.back().tStart) / ((long) clusters.back().qEnd - (long) clusters.back().qStart) >= 5) {
					clusters.pop_back();
				}
				else {
					// debug
					if (opts.debug and clusters.size() >= 2) {
						assert(!(clusters.back().qStart == clusters[clusters.size() - 2].qStart and clusters.back().tStart == clusters[clusters.size() - 2].tStart
							and clusters.back().qEnd == clusters[clusters.size() - 2].qEnd and clusters.back().tEnd == clusters[clusters.size() - 2].tEnd));
					}	

					if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
				}		
			} 
		}
	}
	else {
		return;
	}
}

bool CloseToPreviousCluster(Cluster &a, GenomePos &qS, GenomePos &tS, GenomePos &tE, Options &opts) {
	long aDiff = abs((long)qS - (long)a.qEnd);
	long bDiff;
	if (a.strand == 0) bDiff = abs((long)tS - (long)a.tEnd);
	else bDiff = abs((long)a.tStart - (long)tE);	
	long aDiag = 0, bDiag = 0;	
	if (a.strand == 0) {aDiag = (long)a.tEnd - (long)a.qEnd; bDiag = (long)tS - (long)qS;}
	else {aDiag = (long) a.qEnd + (long) a.tStart; bDiag = (long) qS + (long) tE;}
	if (min(aDiff, bDiff) <= opts.RoughClustermaxGap and abs(aDiag - bDiag) < opts.maxDiag) return true;
	return false;
}

void UpdateCluster(Cluster &a, GenomePos &qS, GenomePos &qE, GenomePos &tS, GenomePos &tE) {
	a.qStart = min(a.qStart, qS);
	a.qEnd = max(a.qEnd, qE);
	a.tStart = min(a.tStart, tS);
	a.tEnd = max(a.tEnd, tE);
}

void MergeTwoClusters(Cluster &a, GenomePos &qS, GenomePos &qE, GenomePos &tS, GenomePos &tE, int &start, int &end) {
	UpdateCluster(a, qS, qE, tS, tE);
	for (int q = start; q < end; q++) { a.splitmatchindex.push_back(q); }
	a.end = end;	
}

void 
SplitRoughClustersWithGaps(vector<pair<GenomeTuple, GenomeTuple> > &matches, Cluster &OriginalClusters, vector<Cluster> &split, 
				Options &opts, int &outIter, Read &read, int strand=0) {

	if (OriginalClusters.end - OriginalClusters.start == 0) {
		return;
	}
	if (OriginalClusters.anchorfreq >= 10.0f) {
		split.push_back(Cluster(OriginalClusters.start, OriginalClusters.end, OriginalClusters.qStart, OriginalClusters.qEnd, 
								OriginalClusters.tStart, OriginalClusters.tEnd, OriginalClusters.strand, outIter));
		for (int q = OriginalClusters.start; q < OriginalClusters.end; q++) { split.back().splitmatchindex.push_back(q);}	
		split.back().anchorfreq = OriginalClusters.anchorfreq;
		return;		
	}
	int cur_s = split.size();
	int split_cs = OriginalClusters.start;
	GenomePos split_qStart = matches[split_cs].first.pos,
			  split_tStart = matches[split_cs].second.pos,
	          split_qEnd = split_qStart + opts.globalK, 
	          split_tEnd = split_tStart + opts.globalK;
	// if (read.name == "S3_12348!22!16246770!16262110!-") cerr << "opts.RoughClustermaxGap: " << opts.RoughClustermaxGap << " opts.minClusterSize: " << opts.minClusterSize << endl; 
	for (int m = OriginalClusters.start + 1; m < OriginalClusters.end; m++) {
		int gap = maxGapDifference(matches[m], matches[m-1]);
		// int diag_gap=DiagonalDifference(matches[m], matches[m-1], strand);

		if (gap > opts.RoughClustermaxGap /*or diag_gap > opts.cleanMaxDiag*/) {
			// cerr << "GAP: " << gap << "\t" << m << "\t" << clusters[c].matches.size() << "\t" << clusters[c].matches[m].second.pos - clusters[c].matches[m-1].second.pos << "\t" << clusters[c].matches[m].first.pos - clusters[c].matches[m-1].first.pos << endl;
			if ((m - split_cs) > opts.minClusterSize) {
				if (split.size() > cur_s and CloseToPreviousCluster(split.back(), split_qStart, split_tStart, split_tEnd, opts)) {
					MergeTwoClusters(split.back(), split_qStart, split_qEnd, split_tStart, split_tEnd, split_cs, m);
				}
				else {
					split.push_back(Cluster(split_cs, m, split_qStart, split_qEnd, split_tStart, split_tEnd, OriginalClusters.strand, outIter));
					split.back().anchorfreq = OriginalClusters.anchorfreq;
					for (int q = split_cs; q < m; q++) { split.back().splitmatchindex.push_back(q);}
				}
			}
			else if (opts.debug and opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
				ofstream fclust("SplitRoughClusters_removed.dots", ofstream::app);
				for (int t = split_cs; t < m; t++) {
					fclust << matches[t].first.pos << "\t" << matches[t].second.pos << "\t" << opts.globalK + matches[t].first.pos << "\t"
							<< matches[t].second.pos + opts.globalK << "\t" << outIter << "\t" << strand << endl;
				}
				fclust.close();
			}
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
	int last = OriginalClusters.end;
	if ((last - split_cs) > opts.minClusterSize) {
		if (split.size() > cur_s and CloseToPreviousCluster(split.back(), split_qStart, split_tStart, split_tEnd, opts)) {
			MergeTwoClusters(split.back(), split_qStart, split_qEnd, split_tStart, split_tEnd, split_cs, last);
		}
		else {
			split.push_back(Cluster(split_cs, last, split_qStart, split_qEnd, split_tStart, split_tEnd, OriginalClusters.strand, outIter));
			split.back().anchorfreq = OriginalClusters.anchorfreq;
			for (int q = split_cs; q < last; q++) { split.back().splitmatchindex.push_back(q);}			
		}
	}
}

// template<typename Tup>
// void StoreDiagonalClusters(vector<pair<Tup, Tup> > &matches, vector<float> & matches_freq, vector<Cluster> &clusters, Options &opts, int s, int e, int strand=0) {
// 	int cs = s, ce = e;
// 	float totalfreq = 0.0f;
// 	while (cs < e) {
// 		ce = cs+1;
// 		GenomePos qStart = matches[cs].first.pos, 
// 				  qEnd = matches[cs].first.pos + opts.globalK, 
// 			      tStart = matches[cs].second.pos, 
// 			      tEnd = matches[cs].second.pos + opts.globalK;
// 		totalfreq += matches_freq[cs];
// 		while (ce < e and abs(DiagonalDifference(matches[ce], matches[ce-1], strand)) < opts.maxDiag) {
// 			qStart = min(qStart, matches[ce].first.pos);
// 			qEnd   = max(qEnd, matches[ce].first.pos + opts.globalK);
// 			tStart = min(tStart, matches[ce].second.pos);
// 			tEnd   = max(tEnd, matches[ce].second.pos + opts.globalK);
// 			totalfreq += matches_freq[ce];
// 			ce++;
// 		}	
// 		if (ce - cs >= opts.minClusterSize and qEnd - qStart >= opts.minClusterLength and tEnd - tStart >= opts.minClusterLength) {
// 			clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand));
// 			clusters.back().anchorfreq = totalfreq / (ce - cs);
// 			totalfreq = 0.0f;
// 		}
// 		cs=ce;
// 	}	
// }

template<typename Tup>
int RemoveSuperRepetitiveClusters (int s, int e, vector<pair<Tup, Tup> > &matches) {
	Tuple t = matches[s].first.t;
	for (int m = s + 1; m < e; m++) {
		if (matches[m].first.t != t) return false;
	}
	return true;
}

template<typename Tup>
void StoreDiagonalClusters(Genome &genome, vector<float> &matches_freq, vector<pair<Tup, Tup> > &matches, vector<Cluster> &clusters, Options &opts, int s, int e, int strand=0) {
	int cs = s, ce = e;
	float totalfreq = 0.0f;
	while (cs < e) {
		ce = cs+1;
		GenomePos qStart = matches[cs].first.pos, 
				  qEnd = matches[cs].first.pos + opts.globalK, 
			      tStart = matches[cs].second.pos, 
			      tEnd = matches[cs].second.pos + opts.globalK;
		totalfreq += matches_freq[cs];

		while (ce < e and abs(DiagonalDifference(matches[ce], matches[ce-1], strand)) < opts.maxDiag) {
			qStart = min(qStart, matches[ce].first.pos);
			qEnd   = max(qEnd, matches[ce].first.pos + opts.globalK);
			tStart = min(tStart, matches[ce].second.pos);
			tEnd   = max(tEnd, matches[ce].second.pos + opts.globalK);
			totalfreq += matches_freq[ce];
			ce++;
		}	
		if (ce - cs >= opts.minClusterSize and qEnd - qStart >= opts.minClusterLength and tEnd - tStart >= opts.minClusterLength
			and !RemoveSuperRepetitiveClusters(cs, ce, matches)) {
			if (opts.bypassClustering) { // needs to insert matches
				clusters.push_back(Cluster(cs, ce, qStart, qEnd, tStart, tEnd, strand));
				for (int b = cs; b < ce; b++) {
					clusters.back().matches.push_back(matches[b]);
				}
				clusters.back().anchorfreq = totalfreq / (ce - cs);
				clusters.back().SetClusterBoundariesFromMatches(opts);
				int rci = genome.header.Find(clusters.back().tStart);
				clusters.back().chromIndex = rci;
				totalfreq = 0.0f;				
			}
			else {
				clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand));
				clusters.back().anchorfreq = totalfreq / (ce - cs);
				totalfreq = 0.0f;				
			}
		}
		cs=ce;
	}	
}


void MatchesToFineClusters (vector<GenomePair> &Matches, vector<Cluster> &clusters, Genome &genome, Read &read, Options &opts, Timing &timing, bool ma_strand = 0) {
	if (read.unaligned) return;
	if (ma_strand == 0) {
		//
		// Guess that 500 is where the setup/takedown overhead of an array index is equal to sort by value.
		// sort fragments in allMatches by forward diagonal, then by first.pos(read)
		//
		DiagonalSort<GenomeTuple>(Matches, 500);
		vector<float> Matches_freq(Matches.size(), 1);
		CleanOffDiagonal(Matches, Matches_freq, opts, read);
		timing.Tick("CleanOffDiagonal");

		//
		// Matches --> roughclusters --> split_roughClusters --> fineclusters 
		//
		vector<Cluster> roughClusters;
		vector<Cluster> split_roughClusters;
		StoreDiagonalClusters(genome, Matches_freq, Matches, roughClusters, opts, 0, Matches.size()); // rough == true means only storing "start and end" in every clusters[i]
		for (int c = 0; c < roughClusters.size(); c++) {
			CartesianSort(Matches, roughClusters[c].start, roughClusters[c].end);
			SplitRoughClustersWithGaps(Matches, roughClusters[c], split_roughClusters, opts, c, read, 0);
		}
		timing.Tick("roughclusters");

		//cerr << "roughClusters.size(): " << roughClusters.size() << " split_roughClusters.size(): " << split_roughClusters.size()<< endl;
		if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
			ofstream fclust("for-matches.dots");
			for (int m = 0; m < Matches.size(); m++) {
				fclust << Matches[m].first.pos << "\t" << Matches[m].second.pos << "\t" << opts.globalK + Matches[m].first.pos << "\t"
						<< Matches[m].second.pos + opts.globalK << "\t" << m << endl;
			}
			fclust.close();

			ofstream rclust("roughClusters.dots");
			for (int m = 0; m < roughClusters.size(); m++) {
				int rci = genome.header.Find(roughClusters[m].tStart);
				for (int c = roughClusters[m].start; c < roughClusters[m].end; ++c) {
					rclust << Matches[c].first.pos << "\t" << Matches[c].second.pos << "\t" << opts.globalK + Matches[c].first.pos << "\t"
						   << Matches[c].second.pos  + opts.globalK << "\t" << m << "\t" << genome.header.names[rci] << endl;				
				}
			}
			rclust.close();

			ofstream wclust("split_roughClusters.dots");
			for (int m=0; m < split_roughClusters.size(); m++) {
				int rci = genome.header.Find(split_roughClusters[m].tStart);
				for (int t = 0; t < split_roughClusters[m].splitmatchindex.size(); t++) {
					int c = split_roughClusters[m].splitmatchindex[t];
					wclust << Matches[c].first.pos << "\t" << Matches[c].second.pos << "\t" << opts.globalK + Matches[c].first.pos << "\t"
						   << Matches[c].second.pos + opts.globalK << "\t" << m << "\t" << split_roughClusters[m].coarse << "\t"<< genome.header.names[rci] << endl;				
				}
			}
			wclust.close();
		}	
		for (int c = 0; c < split_roughClusters.size(); c++) {
			int rci = genome.header.Find(split_roughClusters[c].tStart);
			StoreFineClusters(rci, Matches, split_roughClusters[c].anchorfreq, clusters, opts, split_roughClusters[c].splitmatchindex, genome, read, read.length, 0, c);
		}
		Matches_freq.clear();
		timing.Tick("fineclusters");
	}
	else {
		AntiDiagonalSort<GenomeTuple>(Matches, 500);
		vector<float> Matches_freq(Matches.size(), 1);
		CleanOffDiagonal(Matches, Matches_freq, opts, read, 1);
		timing.Tick("CleanOffDiagonal");

		vector<Cluster> revroughClusters;
		vector<Cluster> split_revroughClusters;
		StoreDiagonalClusters(genome, Matches_freq, Matches, revroughClusters, opts, 0, Matches.size(), 1);
		for (int c = 0; c < revroughClusters.size(); c++) {
			CartesianSort(Matches, revroughClusters[c].start, revroughClusters[c].end);
			SplitRoughClustersWithGaps(Matches, revroughClusters[c], split_revroughClusters, opts, c, read, 1);
		}
		timing.Tick("roughclusters");

		// cerr << "revroughClusters.size(): " << revroughClusters.size() << " split_revroughClusters.dots: " << split_revroughClusters.size()<< endl;
		if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
			ofstream rclust("rev-matches.dots");
			for (int m=0; m < Matches.size(); m++) {			
				rclust << Matches[m].first.pos << "\t" << Matches[m].second.pos + opts.globalK << "\t" << opts.globalK + Matches[m].first.pos  << "\t"
						 << Matches[m].second.pos << "\t" << m << endl;
			}
			rclust.close();

			ofstream revsclust("revroughClusters.dots");
			for (int m=0; m < revroughClusters.size(); m++) {
				int rci = genome.header.Find(revroughClusters[m].tStart);
				for (int c = revroughClusters[m].start; c < revroughClusters[m].end; ++c) {
					revsclust << Matches[c].first.pos << "\t" << Matches[c].second.pos + opts.globalK << "\t" << opts.globalK + Matches[c].first.pos << "\t"
							 << Matches[c].second.pos << "\t" << m << "\t" << genome.header.names[rci] << endl;				
				}
			}
			revsclust.close();

			ofstream wclust("split_revroughClusters.dots");
			for (int m=0; m < split_revroughClusters.size(); m++) {
				int rci = genome.header.Find(split_revroughClusters[m].tStart);
				for (int t = 0; t < split_revroughClusters[m].splitmatchindex.size(); t++) {
					int c = split_revroughClusters[m].splitmatchindex[t];
					wclust << Matches[c].first.pos << "\t" << Matches[c].second.pos << "\t" << opts.globalK + Matches[c].first.pos << "\t"
						   << Matches[c].second.pos + opts.globalK << "\t" << m << "\t" << split_revroughClusters[m].coarse << "\t" << genome.header.names[rci] << endl;				
				}
			}
			wclust.close();
		}
		for (int c = 0; c < split_revroughClusters.size(); c++) {
			int rci = genome.header.Find(split_revroughClusters[c].tStart);
			StoreFineClusters(rci, Matches, split_revroughClusters[c].anchorfreq, clusters, opts, split_revroughClusters[c].splitmatchindex, genome, read, read.length, 1, c);
		}	
		Matches_freq.clear();
		timing.Tick("fineclusters");
	}
}


void CleanMatches (vector<GenomePair> &Matches, vector<Cluster> &clusters, Genome &genome, Read &read, Options &opts, Timing &timing, bool ma_strand = 0) {
	if (read.unaligned) return;
	if (ma_strand == 0) {
		//
		// Guess that 500 is where the setup/takedown overhead of an array index is equal to sort by value.
		// sort fragments in allMatches by forward diagonal, then by first.pos(read)
		//
		DiagonalSort<GenomeTuple>(Matches, 500);
		if (opts.dotPlot and read.name == opts.readname) {
			ofstream fclust("for-matches.dots");
			for (int m = 0; m < Matches.size(); m++) {
				fclust << Matches[m].first.pos << "\t" << Matches[m].second.pos << "\t" << opts.globalK + Matches[m].first.pos << "\t"
						<< Matches[m].second.pos + opts.globalK << "\t" << m << endl;
			}
			fclust.close();
		}
		vector<float> Matches_freq(Matches.size(), 1);
		CleanOffDiagonal(Matches, Matches_freq, opts, read);
		StoreDiagonalClusters(genome, Matches_freq, Matches, clusters, opts, 0, Matches.size(), 0); 
		timing.Tick("CleanOffDiagonal");

		if (opts.dotPlot and read.name == opts.readname) {
			ofstream fclust("for-matches_clean.dots");
			for (int m = 0; m < Matches.size(); m++) {
				fclust << Matches[m].first.pos << "\t" << Matches[m].second.pos << "\t" << opts.globalK + Matches[m].first.pos << "\t"
						<< Matches[m].second.pos + opts.globalK << "\t" << m << endl;
			}
			fclust.close();
		}	
	}
	else {
		AntiDiagonalSort<GenomeTuple>(Matches, 500);
		if (opts.dotPlot and read.name == opts.readname) {
			ofstream rclust("rev-matches.dots");
			for (int m=0; m < Matches.size(); m++) {			
				rclust << Matches[m].first.pos << "\t" << Matches[m].second.pos + opts.globalK << "\t" << opts.globalK + Matches[m].first.pos  << "\t"
						 << Matches[m].second.pos << "\t" << m << endl;
			}
			rclust.close();
		}
		vector<float> revMatches_freq(Matches.size(), 1);
		CleanOffDiagonal(Matches, revMatches_freq, opts, read, 1);
		StoreDiagonalClusters(genome, revMatches_freq, Matches, clusters, opts, 0, Matches.size(), 1); 
		timing.Tick("CleanOffDiagonal");
		// cerr << "revroughClusters.size(): " << revroughClusters.size() << " split_revroughClusters.dots: " << split_revroughClusters.size()<< endl;
		if (opts.dotPlot and read.name == opts.readname) {
			ofstream rclust("rev-matches_clean.dots");
			for (int m=0; m < Matches.size(); m++) {			
				rclust << Matches[m].first.pos << "\t" << Matches[m].second.pos + opts.globalK << "\t" << opts.globalK + Matches[m].first.pos  << "\t"
						 << Matches[m].second.pos << "\t" << m << endl;
			}
			rclust.close();
		}
	
	}
}

// void CheckTrueIntervalInFineCluster(vector<Cluster> &clusters, string &rname, Genome &genome, Read &read) {
// 	if (read.unaligned) return;
// 	string chrom;
// 	GenomePos i_s, i_e;
// 	string bname = rname;
// 	regex re("[!]");
// 	sregex_token_iterator first{bname.begin(), bname.end(), re, -1}, last;//the '-1' is what makes the regex split (-1 := what was not matched)
// 		vector<string> tokens{first, last};
// 	chrom = tokens[1];
// 	string::size_type sz;
// 	i_s = stoi(tokens[2], &sz);
// 	i_e = stoi(tokens[3], &sz);

// 	ofstream cclust("readsCheckTrueIntervalInFineCluster.txt", ofstream::app);
// 	bool check = false;
// 	for (int i = 0; i < clusters.size(); i++) {
// 		if (chrom == genome.header.names[clusters[i].chromIndex]) {
// 			GenomePos offset = genome.header.pos[clusters[i].chromIndex];
// 			if ((i_e + offset >= clusters[i].tStart and i_s + offset <= clusters[i].tEnd)) {
// 				int o = min(i_e + offset, clusters[i].tEnd) - max(i_s + offset, clusters[i].tStart);
// 				if ((float) o / (i_e - i_s) >= 0.2) {
// 					check = true;
// 				}
// 			}
// 		}		
// 	}
// 	if (check) {
// 		cclust << rname << endl;
// 	}
// 	cclust.close();
// }
#endif
