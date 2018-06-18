#ifndef SORTING_H_
#define SORTING_H_

#include <algorithm>

template<typename Tup> 
class DiagonalSortOp {
 public:
	int operator()(const pair<Tup, Tup> &a, const pair<Tup, Tup> &b) {
		int aDiag = a.first.pos - a.second.pos, 
			bDiag= b.first.pos - b.second.pos;

		if (aDiag != bDiag) {
			return aDiag < bDiag;
		}
		else {
			return a.first.pos < b.first.pos;
		}
	}
};


template<typename Tup>
void DiagonalSort(vector<pair<Tup, Tup> > &vals) {
	sort(vals.begin(), vals.end(), DiagonalSortOp<Tup>());
}

template<typename Tup> 
class CartesianSortOp {
 public:
	int operator()(const pair<Tup, Tup> &a, const pair<Tup, Tup> &b) {
		if (a.first.pos != b.first.pos) {
			return a.first.pos < b.first.pos;
		}
		else {
			return a.second.pos < b.second.pos;
		}
	}
};
template<typename Tup>
void CartesianSort(vector<pair<Tup, Tup> > &vals) {
	sort(vals.begin(), vals.end(), CartesianSortOp<Tup>());
}
		
template<typename Tup> 
class CartesianTargetSortOp {
 public:
	int operator()(const pair<Tup, Tup> &a, const pair<Tup, Tup> &b) {
		if (a.second.pos != b.second.pos) {
			return a.second.pos < b.second.pos;
		}
		else {
			return a.first.pos < b.first.pos;
		}
	}
};
template<typename Tup>
void CartesianTargetSort(vector<pair<Tup, Tup> > &vals) {
	sort(vals.begin(), vals.end(), CartesianTargetSortOp<Tup>());
}

#endif
