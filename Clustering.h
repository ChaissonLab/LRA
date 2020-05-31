#ifndef CLUSTERING_H_
#define CLUSTERING_H_
#include "Options.h"
#include <vector>
#include <climits>
#include <boost/icl/interval_map.hpp>
using namespace boost::icl;
using namespace std;
template<typename Tup>
int64_t DiagonalDifference(Tup &a, Tup &b, int strand=0) {
	if (strand == 0) { // forMatches
		int64_t aDiag = (int64_t)a.first.pos - (int64_t)a.second.pos, 
			bDiag = (int64_t)b.first.pos - (int64_t)b.second.pos;
		return aDiag - bDiag;		
	}
	else { // revMathches
		int64_t aDiag = a.first.pos + a.second.pos, 
			bDiag= b.first.pos + b.second.pos;
		return aDiag - bDiag;				
	}
}


template<typename Tup>
int64_t GapDifference(Tup &a, Tup &b) {
	int64_t aDiff = abs((int)b.first.pos - (int)a.first.pos);
	int64_t bDiff = abs((int)b.second.pos - (int)a.second.pos);
	return max(aDiff, bDiff);
}

template<typename Tup>
int DiagonalDrift(int curDiag, Tup &t, int strand=0) {
	int drift;
	if (strand == 0) drift= abs(curDiag - ((int)t.first.pos - (int)t.second.pos));
	else drift= abs(curDiag - ((int)t.first.pos + (int)t.second.pos));
	return drift;
}


template<typename Tup>
void CleanOffDiagonal(vector<pair<Tup, Tup> > &matches, Options &opts, int &minDiagCluster, int strand=0, int diagOrigin=-1, int diagDrift=-1) {
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
	minDiagCluster=0.05*Largest_ClusterNum; // used to be 5, but if use smaller kmer, should be larger


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
		//cout << "encompass: " << float(qovp)/(b.qEnd-b.qStart) << "\t" << float(tovp)/(b.tEnd-b.tStart) << endl;
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
		//		cout << "ovp: " << denomA << "\t" << denomB << "\t" << ovp << "\t" << ovp/denomB << endl;
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
	long long int maxDiagNum;
	long long int minDiagNum; // maxDiagNum and minDiagNum defines diagonal band boundary of the current cluster
	int coarse; 
	int Val; // Val stores the value of each Cluster;(using in SDP)
	int NumofAnchors; // NumofAnchors stores the # of anchors of each splitcluster 
	vector<int> matchesLengths; // store the length of each anchor 
	bool refined; // refined == 0 means this Cluster has not been refined yet
	bool refinespace; // refinespace == 0 means this Cluster has not been add anchors in the step of RefineBtwnSpace;
	int outerCluster;
	bool used;
	int rank;
	Cluster() { refined=0; coarse=-1; used=0; rank=-1;}
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
		NumofAnchors = 0;
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
		NumofAnchors = 0;
	}
  Cluster(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te, int st, int coa) {
		qStart = qs;
		qEnd = qe;
		tStart = ts;
		tEnd = te;
		strand = st;
		coarse = coa;
		Val = 0;
		used = 0;
		NumofAnchors = 0;
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
		int64_t d=DiagonalDifference(matches[m], matches[m-1], strand);
		cout << matches[m-1].first.pos << "\t" << matches[m].first.pos << "\t" << matches[m-1].second.pos << "\t" << matches[m].second.pos << "\t" << d << endl;
	}
}
 
template<typename Tup>
long GetDiag(pair<Tup, Tup> &match, int strand) {
	if (strand == 0) {
		return (long) match.first.pos - (long) match.second.pos;
	}
	else {
		return (long) match.first.pos + (long) match.second.pos;
	}
}

template<typename Tup>
void StoreFineClusters(vector<pair<Tup, Tup> > &matches, vector<Cluster> &clusters, Options &opts, int s, int e, 
											 Genome &genome,
											 GenomePos readLength, 
											 interval_map<GenomePos, int> &xIntv, 
											 interval_map<GenomePos, int> &yIntv,
											 int strand=0, int outerIteration=0) {

	int localMinClusterSize=3;

	int startClusterIndex=clusters.size();
	if (e==s) {
		return;
	}
	long minDiag=GetDiag(matches[s], strand), maxDiag=GetDiag(matches[s], strand);
	for (int i=s+1; i < e; i++) {
		long diag=GetDiag(matches[i], strand);

		if (diag < minDiag) { minDiag = diag;}
		if (diag > maxDiag) { maxDiag = diag;}
	}
	int binSize=50;
	long span=maxDiag-minDiag;
	//	cout << "max diag "<< maxDiag << " min: " << minDiag << endl;
	//cerr << "span: " << span << " binSize: " << binSize << " minClusterSize: " << opts.minClusterSize << endl;
	vector<int> bins(span/binSize + 1,0);
	int maxBin=0;
	int binTotal=0;
	for (int i=s; i < e; i++) {
		long diag  = GetDiag(matches[i], strand);
		long index = (diag - minDiag) / binSize;
		bins[index] += 1;
		if (bins[index] > maxBin) {
			maxBin=bins[index];
		}
	}

	vector<int> sortedBins=bins;
	sort(sortedBins.begin(), sortedBins.end());
	int cumulativeSum=0;
	int cutoff=2;
	for (int i=sortedBins.size(); i >0; i--) {
		cumulativeSum+= sortedBins[i-1];
		if (float(cumulativeSum)/(e-s) > 0.95) {
			cutoff=sortedBins[i-1];
			break;
		}
	}
	/*
	cout << "pre cutoff " << cutoff << "\t" << e-s << endl;
	for (int i=0; i < bins.size(); i++) {
		cout << "bin " << i << "\t" << bins[i] << endl;
	}
	*/
	for (int i=0; i < bins.size(); i++) {
		if (bins[i] < cutoff) { // TODO:Jingwen: minClusterSize is too small?
			bins[i] = 0;
		}
	}
	if (maxBin < 0.0005 * readLength) {
		//		cout << "Switching to wide bin " << maxBin << "\t" << 0.0005*readLength << endl;
		binSize=200;
		localMinClusterSize=2;
		bins.resize(span/binSize+1);
		std::fill(bins.begin(), bins.end(), 0);
		
		for (int i=s; i < e; i++) {
			long diag  = GetDiag(matches[i], strand);
			long index = (diag - minDiag) / binSize;
			bins[index] += 1;
			if (bins[index] > maxBin) {
				maxBin=bins[index];
			}
		}
	}
	for (int i=0; i < bins.size(); i++) {
		if (bins[i] < localMinClusterSize) { // TODO:Jingwen: minClusterSize is too small?
			bins[i] = 0;
		}
	}
	/*	
	cout << "bin size: " << binSize << endl;
	for (int i=0; i<bins.size(); i++) {
		cout << "bin " << i << "\t" << bins[i] << endl;
	}
	*/
	vector<int> diagStart, diagEnd, diagSize;
	std::map<int, int> diagToCluster;
	long curDiagIndex=-1;
	int curCluster=-1;
	for (int i=s; i < e; i++) {
		long diag  = GetDiag(matches[i], strand);
		long index = (diag - minDiag) / binSize;
		//
		// Low count bins were previously filtered out to prevent spurious clusters from forming.
		//
		if ( bins[index] > 0 ) {

			//
			// Need to store this match. curClusterIndex points to the diagonal index
			// where matches are being stored. 
			//
			
			int minDistance = -1;
			int bestDiag = -1;
			int minGap=-1;
			int minDiag=-1;
			int auxMinDist=-1;
			//			cout << "Search " << i << endl;
			for (int clSearch=startClusterIndex; clSearch < clusters.size(); clSearch++) {
				if ( clusters[clSearch].matches.size() > 0 ) {
					int gap  = GapDifference(matches[i], clusters[clSearch].matches[clusters[clSearch].matches.size()-1]);
					int diag = DiagonalDifference(matches[i], clusters[clSearch].matches[clusters[clSearch].matches.size()-1], strand);
					float diagRatio = fabs(diag/float(gap));
					//					cout << "\t" << clSearch << "\t" << gap << "\t" << diag << "\t" << diagRatio << endl;

					if ( gap < 1000 and (diagRatio < 0.2 or gap < 200))  {
						if (minDistance == -1 or gap < minDistance ) {
							minDistance = gap;
							bestDiag    = clSearch;
						}
					}
					if (minGap == -1 or minGap < abs(gap)) {
						minGap = gap;
					}
					if (minDiag == -1 or minDiag < abs(diag)) {
						minDiag = diag;
					}
				}
			}
			if (minDistance != -1 and minDistance < 4000) {
				assert(bestDiag != -1);
				curDiagIndex         = bestDiag;
				diagToCluster[index] = curDiagIndex;
				curCluster           = curDiagIndex;
			}
			else {
				/*
					cout << "Creating a new cluster because of diagonal drift or no cluster found " 
						 << "\t" << minGap << "\t" << minDiag << "\t" << clusters.size() << endl;
				*/
				diagToCluster[index] = clusters.size();
				curDiagIndex = index;
				
				curCluster   = clusters.size();								
				clusters.push_back(Cluster(0, 0, matches[i].first.pos, 
																	 matches[i].first.pos + opts.globalK, 
																	 matches[i].second.pos, 
																	 matches[i].second.pos + opts.globalK, strand));
				clusters[clusters.size()-1].outerCluster = outerIteration;
				clusters[clusters.size()-1].clusterIndex=clusters.size()-1;
			}
			/*					
			if ( curDiagIndex != index) {
				//
				// The diagonal at index is not the one currently pointed to
				// by curDiagIndex, need to find where to store matches.  
				//
				if (diagToCluster.find(index) != diagToCluster.end() ) {
					//
					// This diagonal points to a cluster, reset curDiagIndex
					// to use this index to point to the current cluster.
					//
					curDiagIndex = index;
					curCluster   = diagToCluster[index];
				}
				else {
					//
					// Need to find out which diagonal to point a cluster to.
					//
					if (curDiagIndex >=0 and abs(index - curDiagIndex) <= 3) {
						assert(curCluster != -1);
						diagToCluster[index] = curCluster;
						curDiagIndex = index;
					}
					else {
						//
						// Too far off, need to start a new cluster.

					}
				}
			*/
				
			int clusterSize;
			assert(curCluster != -1);
			
			if (clusters[curCluster].matches.size() > 0) {
				int last=clusters[curCluster].size()-1;
				if (matches[i].first.pos != clusters[curCluster].matches[last].first.pos) {
					assert(matches[i].first.pos >= clusters[curCluster].matches[last].first.pos);
				}
				else {
					assert(matches[i].second.pos >= clusters[curCluster].matches[last].second.pos);
				}
			}
			clusters[curCluster].matches.push_back(matches[i]);

			/*
			cout << curCluster << "\t" << curDiagIndex << "\t" << index <<
			"\t" << clusters[curCluster].matches.size() << "\t" <<
			matches[i].first.pos << "\t" << matches[i].second.pos << endl;
			*/
			// Update endpoint of cluster.
			assert(curCluster < clusters.size());
			clusters[curCluster].qEnd = max(clusters[curCluster].qEnd, matches[i].first.pos + opts.globalK);
			clusters[curCluster].tEnd = max(clusters[curCluster].tEnd, matches[i].second.pos + opts.globalK);
			clusters[curCluster].qStart = min(clusters[curCluster].qStart, matches[i].first.pos);
			clusters[curCluster].tStart = min(clusters[curCluster].tStart, matches[i].second.pos);
			
			string chromName;
			long chromPos;
			genome.GlobalPosToChrom(clusters[curCluster].tStart, chromPos, chromName);
			/*
			cout << curCluster << "\t" << index << "\t"  
					 << s << "-" << i << "-" << e << "\t" << minDistance << "\t"
					 << e-s << "\tcoords\t"
           << matches[i].first.pos << "\t" << matches[i].second.pos << "\tcluster:\t"
					 << clusters[curCluster].qStart << "\t" 
					 << clusters[curCluster].qEnd << "\t" 
					 << clusters[curCluster].qEnd - clusters[curCluster].qStart << "\t" 
					 << clusters[curCluster].tStart << "\t" 
					 << clusters[curCluster].tEnd << "\t"
					 << clusters[curCluster].tEnd - clusters[curCluster].tEnd << "\t" 
					 << chromName << "\t" << chromPos << "\t"
					 << clusters[curCluster].matches.size() << "\t"
					 << outerIteration << endl;
			*/
			assert(clusters[curCluster].tEnd >= clusters[curCluster].tStart);
			assert(clusters[curCluster].qEnd >= clusters[curCluster].qStart);
		}
	}
	for (int c=0; c < clusters.size(); c++) {
		assert(clusters[c].tStart >= 0);
		assert(clusters[c].tEnd >= 0);
	}

	int cn=startClusterIndex;

	bool sufficientSize=false;
	for (int c=startClusterIndex; c < clusters.size(); c++) {
		if (clusters[c].matches.size() > localMinClusterSize) {
			//			cout << "adding intervals " << clusters[c].qStart << "\t" << clusters[c].qEnd << "\t" << clusters[c].tStart << "\t" << clusters[c].tEnd << endl;
			if (clusters[c].matches.size() >= localMinClusterSize*4) {
				sufficientSize=true;
			}
			xIntv.add(make_pair(interval<GenomePos>::right_open(clusters[c].qStart, 
																													clusters[c].qEnd), 
													clusters[c].matches.size()));
			yIntv.add(make_pair(interval<GenomePos>::right_open(clusters[c].tStart, 
																													clusters[c].tEnd), 
													clusters[c].matches.size()));
		}
	}
	int nContained=0;

	for (int c=startClusterIndex; c < clusters.size(); c++) {
		bool sizeThresh=clusters[c].matches.size() > 1;
		bool xIntvS = contains(xIntv, clusters[c].qStart);
		bool xIntvE = contains(xIntv, clusters[c].qEnd-1);
		bool yIntvS = contains(yIntv, clusters[c].tStart);
		bool yIntvE = contains(yIntv, clusters[c].tEnd-1);
		bool contained = xIntvS or xIntvE or yIntvS or yIntvE;
		if (sizeThresh == false and contained == true) { ++nContained;}
		//
		// If the anchor is large enough, or it is 
		if (sizeThresh or (contained == false and sufficientSize == true) ) {
			clusters[cn] = clusters[c];
			cn++;
		}
	}
	//	if (nContained > 20) { cerr << clusters.size() << "\t" << nContained << endl;}
	//  std::cout << outerIteration << "\t" << clusters.size() << "\t" << cn << endl;
	clusters.resize(cn);


			
			//
			// Add this point to a cluster.

	
//	while ( foundCluster == true ) {
//		int maxDiagSize=0;
//		int maxDiag=0;
//		bool clusterStarted = false;
//		for ( int i=0; i < bins.size(); i++ ) {
//			if ( bins[i] > maxDiagSize and bins[i] > opts.minClusterSize ) { 
//				maxDiagSize = bins[i];
//				maxDiag = i;
//			}
//		}
//		if ( maxDiagSize  == 0 ) {
//			break;
//		}
//		else {
//			int j=maxDiag;
//			while (j > 0 and maxDiag - j < 3 and bins[j-1] > opts.minClusterSize) { j--; } //3;
//			int k=maxDiag;
//			while (k < bins.size() and k-maxDiag < 3 and bins[k] > 0) { k++;} //3;
//			int totalSize=0;
//			for (int l=j; l < k; l++) {
//				totalSize+=bins[l];
//				bins[l] = 0;
//			}
//
//			diagSize.push_back(totalSize);
//			diagStart.push_back(j);
//			diagEnd.push_back(k);
//		}
//	}
//	
}

void SplitClustersWithGaps(vector<Cluster> &clusters, vector<Cluster> &split, Options &opts ) {
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
				// s				cout << "GAP: " << gap << "\t" << m << "\t" << clusters[c].matches.size() << "\t" << clusters[c].matches[m].second.pos - clusters[c].matches[m-1].second.pos << "\t" << clusters[c].matches[m].first.pos - clusters[c].matches[m-1].first.pos << endl;
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
void StoreDiagonalClusters(vector<pair<Tup, Tup> > &matches, vector<Cluster> &clusters, Options &opts, int s, int e, 
								bool rough, bool lite, int strand=0) {

	int maxGap = -1, maxDiag = -1;
	if (rough == false) { maxGap = opts.maxGapBtwnAnchors;}

	else { maxDiag = opts.maxDiag;} 
	int i;
	int cs = s, ce = e;
	int64_t dd,absdd;
	
	while (cs < e) {
		ce = cs+1;
		GenomePos qStart=matches[cs].first.pos, 
				  qEnd=matches[cs].first.pos + opts.globalK, 
			      tStart=matches[cs].second.pos, 
			      tEnd=matches[cs].second.pos + opts.globalK;
		int diff=0;

		// (TODO)Jingwen: Delete (opts.maxGap == -1 or GapDifference(matches[ce], matches[ce-1]) < opts.maxGap) in the below
		while (ce < e and (abs(DiagonalDifference(matches[ce], matches[ce-1], strand)) < opts.maxDiag or maxDiag == -1) 
					  and (maxGap == -1 or GapDifference(matches[ce], matches[ce-1]) < maxGap)) {

			qStart = min(qStart, matches[ce].first.pos);
			qEnd   = max(qEnd, matches[ce].first.pos + opts.globalK);
			tStart = min(tStart, matches[ce].second.pos);
			tEnd   = max(tEnd, matches[ce].second.pos + opts.globalK);
			/*		cout << "rc:\t" << clusters.size() << "\t" << matches[ce].first.pos << "\t" << matches[ce].second.pos << "\t" 
					 << GetDiag(matches[ce], strand) << "\t" <<
					 abs(DiagonalDifference(matches[ce], matches[ce-1], strand))
					 << endl;
			*/
			ce++;
		}	
		
		if (ce - cs >= opts.minClusterSize and qEnd - qStart >= opts.minClusterLength and tEnd - tStart >= opts.minClusterLength) {
			if (rough == true) {
				clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand));
			}
			else {
				if (lite == false ) {
					clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand, matches.begin() + cs, matches.begin()+ce));
				}
				else {
					clusters.push_back(Cluster(cs,ce, qStart, qEnd, tStart, tEnd, strand));
				}
			}
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
