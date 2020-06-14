#ifndef SPLIT_CLUSTERS_H_
#define SPLIT_CLUSTERS_H_

#include <iostream> //std::cout 
#include <fstream>   
#include <cstdlib>   // std::labs, std::EXIT FAILURE, std::EXIT SUCCESS
#include <vector>
#include <set>
#include <cstdint>
#include <utility>

using namespace std;

typedef int64_t Ngp;


class IntervalSet {
public:
	double slope;
	double intercept;
	bool strand;
	vector<pair<GenomePos, bool>> Set;

	IntervalSet (Cluster & cluster) {
		slope = (double)((Ngp)cluster.tEnd - (Ngp)cluster.tStart)/((Ngp)cluster.qEnd - (Ngp)cluster.qStart);
		if (cluster.strand == 0) {
			intercept = ((double)((Ngp)cluster.qEnd*cluster.tStart - (Ngp)cluster.qStart*cluster.tEnd))/
												((Ngp) cluster.qEnd - (Ngp) cluster.qStart);
		}
		else {
			slope = -1*slope;
			intercept = (double)((Ngp)cluster.qStart*cluster.tStart - (Ngp)cluster.qEnd*cluster.tEnd)/((Ngp)cluster.qStart - (Ngp)cluster.qEnd);
		}
		strand = cluster.strand;
	};

	~IntervalSet() {};

	int operator() (const pair<GenomePos, bool> &a, const pair<GenomePos, bool> &b) {
		if (a.second == b.second and a.second == 0) {
			return a.first < b.first;
		}
		else if (a.second == b.second and a.second == 1) {
			if (strand == 0) return a.first < b.first;
			else return a.first > b.first;
		}
		else if (a.second == 0 and b.second == 1) {
			if (strand == 0) return a.first*slope + intercept < (double) b.first;
			else return a.first*slope + intercept > (double) b.first;
		}
		else {
			if (strand == 0) return (double) a.first < b.first*slope + intercept;
			else return (double) a.first > b.first*slope + intercept;
		}		
	}

	void Sort () {
		sort(Set.begin(), Set.end(), *this);
	}
};


void SplitClusters(vector<Cluster> & clusters, vector<Cluster> & splitclusters) {
	set<GenomePos> qSet;
	set<GenomePos> tSet;

	//
	// insert q/t coordinates of each cluster into qSet/tSet;
	//
	for (int m = 0; m < clusters.size(); m++) {
		qSet.insert(clusters[m].qStart);
		qSet.insert(clusters[m].qEnd);
		tSet.insert(clusters[m].tStart);
		tSet.insert(clusters[m].tEnd);
	} 

	//
	// Find what coordinates appear in the interval of each cluster
	//
	for (int m = 0; m < clusters.size(); m++) {
	
		IntervalSet itlSet(clusters[m]);
		set<GenomePos>::iterator its, ite;
		its = qSet.upper_bound(clusters[m].qStart);
		ite = qSet.lower_bound(clusters[m].qEnd); // this points to clusters[m].tEnd

		for (set<GenomePos>::iterator it = its; it != ite; it++) {
			itlSet.Set.push_back(make_pair(*it, 0)); 
		}

		its = tSet.upper_bound(clusters[m].tStart);
		ite = tSet.lower_bound(clusters[m].tEnd);

		for (set<GenomePos>::iterator it = its; it != ite; it++) {
			itlSet.Set.push_back(make_pair(*it, 1)); 
		}		

		itlSet.Sort();

		//
		// Split clusters[m] 
		//
		pair<GenomePos, GenomePos> prev;
		if (clusters[m].strand == 0) prev = make_pair(clusters[m].qStart, clusters[m].tStart);
		else prev = make_pair(clusters[m].qStart, clusters[m].tEnd); 

		vector<pair<GenomePos, bool>>::iterator it = itlSet.Set.begin();
		for (; it < itlSet.Set.end(); it++) {

			if (it->second == 0) { // split on q coord.
				GenomePos t = (GenomePos) ceil(itlSet.slope * it->first + itlSet.intercept);

				if (prev.first < it->first) { 
					if (clusters[m].strand == 0) 
						splitclusters.push_back(Cluster(prev.first, it->first, prev.second, t, clusters[m].strand, m)); // initialize coarse to specify the index of original index
					else 
						splitclusters.push_back(Cluster(prev.first, it->first, t, prev.second, clusters[m].strand, m));					
				}
				else continue;

				prev = make_pair(it->first, t);					
			}
			else { // split on t coord.
				GenomePos q = (GenomePos) ceil((it->first - itlSet.intercept) / itlSet.slope);
				
				if (prev.first < q) {
					if (clusters[m].strand == 0) 
						splitclusters.push_back(Cluster(prev.first, q, prev.second, it->first, clusters[m].strand, m));
					else 
						splitclusters.push_back(Cluster(prev.first, q, it->first, prev.second, clusters[m].strand, m));

				}
				else continue;

				prev = make_pair(q, it->first);					
			}

		} 

		if (prev.first < clusters[m].qEnd) {
			if (clusters[m].strand == 0) splitclusters.push_back(Cluster(prev.first, clusters[m].qEnd, prev.second, clusters[m].tEnd, clusters[m].strand, m));
			else splitclusters.push_back(Cluster(prev.first, clusters[m].qEnd, clusters[m].tStart, prev.second, clusters[m].strand, m));			
		}

	}
}


void DecideSplitClustersValue (vector<Cluster> & clusters, vector<Cluster> & splitclusters, Options & opts) {
	if (splitclusters.size() == 0) return;
	//
	// Compute matches bases/the Cluster length for each Cluster in clusters;
	//
	for (int m = 0; m < clusters.size(); m++) {

		if (clusters[m].matches.size() == 0) continue;
		GenomePos cur_len = clusters[m].matches[0].first.pos;
		GenomePos MatNum = 0;

		for (int n = 0; n < clusters[m].matches.size(); n++) {

			if (cur_len > clusters[m].matches[n].first.pos) {
				assert(cur_len <= clusters[m].matches[n].first.pos + opts.globalK);
				MatNum += clusters[m].matches[n].first.pos + opts.globalK - cur_len;
			}
			else {
				MatNum += opts.globalK;
			}
			cur_len = clusters[m].matches[n].first.pos + opts.globalK;
		}
		clusters[m].Val = MatNum;
	} 
	//
	// Compute the value for each Cluster in splitclusters;
	//
	for (int m = 0; m < splitclusters.size(); m++) {
		int ic = splitclusters[m].coarse;
		float pika = (float)min(splitclusters[m].qEnd - splitclusters[m].qStart, splitclusters[m].tEnd - splitclusters[m].tStart) / 
										(float)min(clusters[ic].qEnd - clusters[ic].qStart, clusters[ic].tEnd - clusters[ic].tStart);
		splitclusters[m].Val = (int)clusters[ic].Val*pika;
		//cerr << "m: " << m << " ic: " << ic << " length/totallenth: " << pika << 
		//" clusters[ic].Val: " << clusters[ic].Val << endl;
	}
	//
	// Compute # of anchors in each splitcluster
	//
	int m = 0;
	int n = 1;
	int ic_m = splitclusters[m].coarse;

	int ic_n;
	if (splitclusters.size() > n) {
		ic_n = splitclusters[n].coarse;
	}
	
	int matchS = 0, matchE = 0;
	while (n < splitclusters.size()) {
		if (ic_m == ic_n) {	
			matchE = CartesianLowerBound<GenomeTuple>(clusters[ic_n].matches.begin(), 
														  clusters[ic_n].matches.end(), splitclusters[n].qStart);		
			assert(matchE >= matchS);
			splitclusters[m].NumofAnchors = matchE - matchS;		
			matchS = matchE;														  	
		}
		else {
			matchE = clusters[ic_m].matches.size();
			assert(matchE >= matchS);
			splitclusters[m].NumofAnchors = matchE - matchS;		
			matchS = 0;	
		}
		m = n;
		ic_m = ic_n;
		n++;
		if (n < splitclusters.size()) ic_n = splitclusters[n].coarse;
	}
	if (splitclusters.size() > 1) {	
		splitclusters[n-1].NumofAnchors = clusters[ic_m].matches.size() - matchS;	
	}
	
}


#endif

