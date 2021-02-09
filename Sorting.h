#ifndef SORTING_H_
#define SORTING_H_

#include <algorithm>
#include "Types.h"
#include <vector>
#include <iterator>
#include <numeric>


using std::sort;
using std::pair;
using std::tuple;

template<typename Tup> 
class DiagonalSortOp {
 public:
	int operator()(const pair<Tup, Tup> &a, const pair<Tup, Tup> &b) {
		int aDiag = (int)a.first.pos - (int)a.second.pos, 
			bDiag= (int)b.first.pos - (int)b.second.pos;
		if (aDiag != bDiag) {
			return aDiag < bDiag;
		}
		else {
			return a.first.pos < b.first.pos; 
		}
	}
};

// Sort anchors with frequence
template<typename Tup, typename T> 
class DiagonalSortOp_freq {
 public:
	int operator()(const tuple<Tup, Tup, T> &a, const tuple<Tup, Tup, T> &b) {
		int aDiag = (int)get<0>(a).pos - get<1>(b).pos, 
			bDiag= (int)get<0>(b).pos - (int)get<1>(b).pos;
		if (aDiag != bDiag) {
			return aDiag < bDiag;
		}
		else {
			return get<0>(a).pos < get<0>(b).pos; 
		}
	}
};

template<typename Tup, typename T>
class DiagonalIndexSort {
public:
	typename vector<tuple<Tup, Tup, T> >::iterator tuples;
  int operator()(const int &a, const int &b) {
		typename vector<tuple<Tup, Tup, T> >::iterator ap = tuples+a;
		typename vector<tuple<Tup, Tup, T> >::iterator bp = tuples+b;
		int aDiag = (int)get<0>(*ap).pos - (int)get<1>(*ap).pos, 
			bDiag = (int)get<0>(*bp).pos - (int)get<1>(*bp).pos;
		if (aDiag != bDiag) {
			return aDiag < bDiag;
		}
		else {
			return get<0>(*ap).pos < get<0>(*bp).pos; 
		}
	}
};
	
template<typename Tup, typename T>
void DiagonalSort(typename vector<tuple<Tup, Tup, T> >::iterator begin, typename vector<tuple<Tup, Tup, T> >::iterator end, int minRange=0) {
	if (minRange == 0 or end - begin < minRange) {
		sort(begin, end, DiagonalSortOp_freq<Tup, T>());
	}
	else {
		DiagonalIndexSort<Tup, T> sorter;
		sorter.tuples=begin;
		vector<int> index(end - begin);
		std::iota(index.begin(), index.end(), 0);
		sort(index.begin(), index.end(), sorter);
		vector<tuple<Tup, Tup, T> > temp(end - begin);
		copy(begin, end, temp.begin());
		
		for (int i=0; i < index.size(); i++) {
			temp[i] = *(begin + index[i]);
		}
		copy(temp.begin(), temp.end(), begin);
	}	
}

template<typename Tup, typename T>
void DiagonalSort(vector<tuple<Tup, Tup, T> > &vals, int minRange=0) {
	DiagonalSort<Tup, T>(vals.begin(), vals.end(), minRange);
}

template<typename Tup>
void DiagonalSort(typename vector<pair<Tup, Tup> >::iterator begin, typename vector<pair<Tup, Tup> >::iterator end) {
		sort(begin, end, DiagonalSortOp<Tup>());
}

template<typename Tup> 
class AntiDiagonalSortOp {
 public:
 AntiDiagonalSortOp(GenomePos l) : length(l) {}
	GenomePos length;
	int operator()(const pair<Tup, Tup> &a, const pair<Tup, Tup> &b) {
		int aDiag = (int)a.first.pos - (int)(length-a.second.pos),  
			bDiag= (int)b.first.pos - (int)(length-b.second.pos); 

		if (aDiag != bDiag) {
			return aDiag < bDiag;
		}
		else {
			return a.first.pos < b.first.pos; 
		}
	}
};

template<typename Tup, typename T> 
class AntiDiagonalSortOp_freq {
 public:
 AntiDiagonalSortOp_freq(GenomePos l) : length(l) {}
	GenomePos length;
	int operator()(const tuple<Tup, Tup, T> &a, const tuple<Tup, Tup, T> &b) {
		int aDiag = (int)get<0>(a).pos - (int)(length - get<1>(a).pos),  
			bDiag= (int)get<0>(b).pos - (int)(length - get<1>(b).pos); 

		if (aDiag != bDiag) {
			return aDiag < bDiag;
		}
		else {
			return get<0>(a).pos < get<0>(b).pos; 
		}
	}
};

template<typename Tup, typename T>
class AntiDiagonalIndexSort {
public:
	GenomePos length;
	typename vector<tuple<Tup, Tup, T> >::iterator tuples;
  int operator()(const int &a, const int &b) {
		typename vector<std::tuple<Tup, Tup, T> >::iterator ap=tuples+a;
		typename vector<std::tuple<Tup, Tup, T> >::iterator bp=tuples+b;
		int aDiag = (int)get<0>(*ap).pos - (int)(length - get<1>(*ap).pos); 
		int bDiag = (int)get<0>(*bp).pos - (int)(length - get<1>(*bp).pos);
		
		if (aDiag != bDiag) {
			return aDiag < bDiag;
		}
		else {
			return get<0>(*ap).pos < get<0>(*bp).pos; 
		}
	}
};

template<typename Tup, typename T>
void AntiDiagonalSort(typename vector<tuple<Tup, Tup, T> >::iterator  begin, typename vector<tuple<Tup, Tup, T> >::iterator  end, GenomePos genomeLength, int sortByIndex=0) {

	if (sortByIndex == 0 or end-begin < sortByIndex) {
		sort(begin, end, AntiDiagonalSortOp_freq<Tup, T>(genomeLength));
	}
	else {
		AntiDiagonalIndexSort<Tup, T> sorter;
		sorter.tuples=begin;
		sorter.length=genomeLength;
		vector<int> index(end-begin);
		std::iota(index.begin(), index.end(), 0);
		sort(index.begin(), index.end(), sorter);
		GenomePos pos;
		vector<tuple<Tup, Tup, T >> temp(end-begin);
		copy(begin,end, temp.begin());
		
		for (int i = 0; i < index.size(); i++) {
			temp[i] = *(begin+index[i]);
		}
		copy(temp.begin(), temp.end(), begin);
	}	
}

template<typename Tup, typename T>
void AntiDiagonalSort(vector<tuple<Tup, Tup, T> > &vals, GenomePos length, int sortByIndex=0) {
	AntiDiagonalSort<Tup, T>(vals.begin(), vals.end(), length, sortByIndex);
}

template<typename Tup>
void AntiDiagonalSort(typename vector<pair<Tup, Tup> >::iterator  begin, typename vector<pair<Tup, Tup> >::iterator  end, GenomePos genomeLength) {
	sort(begin, end, AntiDiagonalSortOp<Tup>(genomeLength));
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

template<typename Tup, typename T> 
class CartesianSortOp_freq {
 public:
	int operator()(const tuple<Tup, Tup, T> &a, const tuple<Tup, Tup, T> &b) {
		if (get<0>(a).pos != get<0>(b).pos) {
			return get<0>(a).pos < get<0>(b).pos;
		}
		else {
			return get<1>(a).pos < get<1>(b).pos;
		}
	}
};

template<typename Tup, typename T>
void CartesianSort(vector<tuple<Tup, Tup, T> > &vals, int s, int e) {
	sort(vals.begin() + s, vals.begin() + e, CartesianSortOp_freq<Tup, T>());
}

template<typename Tup>
void CartesianSort(typename vector<pair<Tup, Tup> >::iterator begin, typename vector<pair<Tup, Tup> >::iterator end) {
	sort(begin, end, CartesianSortOp<Tup>());
}

template<typename Tup>
void CartesianSort(vector<pair<Tup, Tup> > &vals) {
	CartesianSort<Tup>(vals.begin(), vals.end());
}
		
template<typename Tup>
int CartesianLowerBound(typename vector<pair<Tup, Tup> >::iterator  begin,
						typename vector<pair<Tup, Tup> >::iterator  end, int64_t query) {
	pair<Tup, Tup> queryTup;
	queryTup.first.pos = query;
	queryTup.second.pos = 0;
	return lower_bound(begin, end, queryTup, CartesianSortOp<Tup>()) - begin;
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
void CartesianTargetSort(typename vector<pair<Tup, Tup> >::iterator begin, typename vector<pair<Tup, Tup> >::iterator end) {
	sort(begin, end, CartesianTargetSortOp<Tup>());
}


template<typename Tup>
void CartesianTargetSort(vector<pair<Tup, Tup>> &matches, int s, int e) {
	sort(matches.begin() + s, matches.begin() + e, CartesianTargetSortOp<Tup>());
}


template<typename Tup>
int CartesianTargetLowerBound(typename vector<pair<Tup, Tup> >::iterator  begin,
							 typename vector<pair<Tup, Tup> >::iterator  end, int64_t query) {
	pair<Tup, Tup> queryTup;
	queryTup.second.pos = query;
	queryTup.first.pos = 0;
	return lower_bound(begin, end, queryTup, CartesianTargetSortOp<Tup>()) - begin;
}

template<typename Tup>
int CartesianTargetUpperBound(typename vector<pair<Tup, Tup> >::iterator  begin, typename vector<pair<Tup, Tup> >::iterator  end, 
															int64_t query) {
	pair<Tup, Tup> queryTup;
	queryTup.second.pos = query;
	queryTup.first.pos = 0;
	return upper_bound(begin, end, queryTup, CartesianTargetSortOp<Tup>()) - begin;
}


template<typename T> 
class SortByRowOp {
 public:
	int operator()(const T & a, const T & b) {
		if (a.se.first != b.se.first) {
			return a.se.first < b.se.first;
		}
		else if (a.se.second != b.se.second){
			return a.se.second < b.se.second; 
		}
		else {
			return a.ind < b.ind; 
		}
	}
};


template<typename T1, typename T2> 
class SortByColOp {
 public:

 	SortByColOp(std::vector<T1> & H);// constructor
 	std::vector<T1> * Hp;

	int operator()(const T2 & a, const T2 & b) {
		if ((*Hp)[a].se.second != (*Hp)[b].se.second) {
			return (*Hp)[a].se.second < (*Hp)[b].se.second;
		}
		else if ((*Hp)[a].se.first != (*Hp)[b].se.first){
			return (*Hp)[a].se.first < (*Hp)[b].se.first; 
		}
		else {
			return (*Hp)[a].ind < (*Hp)[b].ind; 
		}
	}
};

// constructor
template<typename T1, typename T2>
SortByColOp<T1, T2>::SortByColOp(std::vector<T1> & H) {
	Hp = & H;
}


template<typename T1, typename T2>
class SortByBackDiagOp
{
public:
 	SortByBackDiagOp(std::vector<T1> & H); // constructor && initialization list

 	std::vector<T1> * Hp;

 	int operator()(const T2 & a, const T2 & b) {
 		long int aBackDiag = (*Hp)[a].se.first + (*Hp)[a].se.second;
 		long int bBackDiag = (*Hp)[b].se.first + (*Hp)[b].se.second;
 		if (aBackDiag != bBackDiag) {
 			return aBackDiag < bBackDiag;
 		}
 		else if ((*Hp)[a].se.first != (*Hp)[b].se.first){
 			return (*Hp)[a].se.first < (*Hp)[b].se.first;
 		}
 		else {
 			return (*Hp)[a].ind < (*Hp)[b].ind;
 		}
 	}
};


// Constructor
template<typename T1, typename T2>
SortByBackDiagOp<T1, T2>::SortByBackDiagOp(std::vector<T1> & H) {
	Hp = & H;
}


// This Lower_bound function return the index of the element
// the first element in the range [first,last) which is greater than or equal to val.
template <typename T1, typename T2>
T1 Lower_Bound (T1 first, T1 last, long int val, std::vector<T2> & E_1) {
	
	T1 it;
	unsigned int count, step;
	count = std::distance(first, last);
	while (count > 0) {
		it = first; step = count/2; std::advance(it, step);
		if ( E_1[*it] < val) {
			first = ++it;
			count -= step + 1;
		}
		else count = step;
	}
	return first;
}


#endif
