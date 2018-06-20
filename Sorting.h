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
void DiagonalSort(typename vector<pair<Tup, Tup> >::iterator  begin,
									typename vector<pair<Tup, Tup> >::iterator  end) {
	sort(begin, end, DiagonalSortOp<Tup>());
}
template<typename Tup>
void DiagonalSort(vector<pair<Tup, Tup> > &vals) {
	DiagonalSort<Tup>(vals.begin(), vals.end());
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
void CartesianSort(typename vector<pair<Tup, Tup> >::iterator  begin,
									 typename vector<pair<Tup, Tup> >::iterator  end) {

	sort(begin, end, CartesianSortOp<Tup>());
}

template<typename Tup>
void CartesianSort(vector<pair<Tup, Tup> > &vals) {
	CartesianSort<Tup>(vals.begin(), vals.end());
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
	CartesianTargetSort<Tup>(vals.begin(), vals.end());
}

template<typename Tup>
void CartesianTargetSort(typename vector<pair<Tup, Tup> >::iterator  begin,
												 typename vector<pair<Tup, Tup> >::iterator  end) {
	sort(begin, end, CartesianTargetSortOp<Tup>());
}


template<typename Tup>
int CartesianTargetLowerBound(typename vector<pair<Tup, Tup> >::iterator  begin,
															typename vector<pair<Tup, Tup> >::iterator  end, 
															int64_t query) {
	pair<Tup, Tup> queryTup;
	queryTup.second.pos = query;
	return lower_bound(begin, end, queryTup, CartesianTargetSortOp<Tup>()) - begin;
}

template<typename Tup>
int CartesianTargetUpperBound(typename vector<pair<Tup, Tup> >::iterator  begin,
															typename vector<pair<Tup, Tup> >::iterator  end, 
															int64_t query) {
	pair<Tup, Tup> queryTup;
	queryTup.second.pos = query;
	return upper_bound(begin, end, queryTup, CartesianTargetSortOp<Tup>()) - begin;
}

#endif
