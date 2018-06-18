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
			return aDiag > bDiag;
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

		
#endif
