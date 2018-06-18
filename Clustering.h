#ifndef CLUSTERING_H_
#define CLUSTERING_H_
#include "Options.h"

template<typename Tup>
int DiagonalDifference(Tup &a, Tup &b) {
		int aDiag = a.first.pos - a.second.pos, 
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
	int c = 0;
	for (int i=0; i < matches.size(); i++) {
		if (keep[i]) {
			matches[c] = matches[i]; c++;
		}
	}

	matches.resize(c);
}

template<typename Tup>
void OffDiagonalClusters(vector<pair<Tup, Tup> > &matches, vector< vector<pair<Tup, Tup> > > &clusters, Options &opts) {
	int i;
	int cs = 0, ce =0;
	cs = 0;
	while (cs < matches.size()) {
		ce = cs+1;
		while (ce < matches.size() and 
					 abs(DiagonalDifference(matches[ce], matches[ce-1])) < opts.maxDiag) {
			ce++;
		}
			
		if (ce - cs >= opts.minClusterSize) {
			clusters.push_back(vector<pair<Tup, Tup> >());
			int l=clusters.size()-1;
			for (int ci = cs; ci < ce; ci++ ){
				clusters[l].push_back(matches[ci]);
			}
		}
		cs=ce;
	}
}
	

#endif
