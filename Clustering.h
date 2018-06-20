#ifndef CLUSTERING_H_
#define CLUSTERING_H_
#include "Options.h"

template<typename Tup>
int64_t DiagonalDifference(Tup &a, Tup &b) {
		int64_t aDiag = a.first.pos - a.second.pos, 
			bDiag= b.first.pos - b.second.pos;
		return aDiag - bDiag;
}


template<typename Tup>
void CleanOffDiagonal(vector<pair<Tup, Tup> > &matches, Options &opts) {
	if (matches.size() == 0) {
		return;
	}
	
	vector<bool> keep(matches.size(), true);

	int m;
	if (matches.size() > 1) {
		if (abs(DiagonalDifference(matches[0], matches[1])) > opts.maxDiag) {
			keep[0] = false;
		}
	}
	for (int i = 1; i < matches.size() - 1; i++) {
		if (abs(DiagonalDifference(matches[i], matches[i-1])) > opts.maxDiag &&
				abs(DiagonalDifference(matches[i], matches[i+1])) > opts.maxDiag) {
			keep[i] = false;
		}
	}
	if (matches.size() > 1) {
		int last = matches.size()-1;
		if (abs(DiagonalDifference(matches[last], matches[last-1])) > opts.maxDiag) {
			keep[last] = false;
		}
	}
	int c   = 0;
	int pre = matches.size();
	for (int i=0; i < matches.size(); i++) {
		if (keep[i]) {
			matches[c] = matches[i]; c++;
		}
	}
	matches.resize(c);
}


class Cluster {
 public:
	int start;
	int end;
	int strand;
	char *seq;
 Cluster(int s, int e) : start(s), end(e), strand(0) {}
	int size() {
		return end - start;
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
void StoreDiagonalClusters(vector<pair<Tup, Tup> > &matches, vector<Cluster > &clusters, Options &opts) {
	int i;
	int cs = 0, ce =0;
	cs = 0;
	int64_t dd,absdd;

	while (cs < matches.size()) {
		ce = cs+1;
		while (ce < matches.size() and 
					 abs(DiagonalDifference(matches[ce], matches[ce-1])) < opts.maxDiag) {
			ce++;
			dd=DiagonalDifference(matches[ce], matches[ce-1]);
			absdd=abs(DiagonalDifference(matches[ce], matches[ce-1]));			
		}			
		if (ce - cs >= opts.minClusterSize) {
			clusters.push_back(Cluster(cs,ce));
		}
		cs=ce;
	}
}
	

#endif
