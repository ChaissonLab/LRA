#ifndef SORTING_H_
#define SORTING_H_

#include <algorithm>
#include "Types.h"
#include <vector>
#include <iterator>
#include <numeric>


using std::sort;
using std::pair;


template<typename Tup>
class DiagonalIndexSort {
public:

	typename vector<pair<Tup, Tup> >::iterator tuples;
  int operator()(const int &a, const int &b) {
		typename vector<std::pair<Tup, Tup> >::iterator ap=tuples+a;
		typename vector<std::pair<Tup, Tup> >::iterator bp=tuples+b;
		int aDiag = (int)ap->first.pos - (int)ap->second.pos, 
			bDiag= (int)bp->first.pos - (int)bp->second.pos;
		
		if (aDiag != bDiag) {
			return aDiag < bDiag;
		}
		else {
			return ap->first.pos < bp->first.pos; 
		}
	}

};
	
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


template<typename Tup>
void DiagonalSort(typename vector<pair<Tup, Tup> >::iterator begin, typename vector<pair<Tup, Tup> >::iterator end, int minRange=0) {
	if (minRange == 0 or end-begin < minRange) {

		sort(begin, end, DiagonalSortOp<Tup>());
	}
	else {
		
		DiagonalIndexSort<Tup> sorter;
		sorter.tuples=begin;
		vector<int> index(end-begin);
		std::iota(index.begin(), index.end(), 0);
		sort(index.begin(), index.end(), sorter);
		GenomePos pos;
		vector<pair<Tup, Tup> > temp(end-begin);
		copy(begin,end, temp.begin());
		
		for (int i=0; i < index.size(); i++) {
			temp[i]=*(begin+index[i]);
		}
		copy(temp.begin(), temp.end(), begin);
	}	
}
template<typename Tup>
void DiagonalSort(vector<pair<Tup, Tup> > &vals, int minRange=0) {
	DiagonalSort<Tup>(vals.begin(), vals.end(), minRange);
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


template<typename Tup>
class AntiDiagonalIndexSort {
public:
	GenomePos length;
	typename vector<pair<Tup, Tup> >::iterator tuples;
  int operator()(const int &a, const int &b) {
		typename vector<std::pair<Tup, Tup> >::iterator ap=tuples+a;
		typename vector<std::pair<Tup, Tup> >::iterator bp=tuples+b;
		int aDiag = (int)ap->first.pos - (int)(length - ap->second.pos); 
		int bDiag = (int)bp->first.pos - (int)(length - bp->second.pos);
		
		if (aDiag != bDiag) {
			return aDiag < bDiag;
		}
		else {
			return ap->first.pos < bp->first.pos; 
		}
	}

};

template<typename Tup>
void AntiDiagonalSort(typename vector<pair<Tup, Tup> >::iterator  begin,
											typename vector<pair<Tup, Tup> >::iterator  end, GenomePos genomeLength, int sortByIndex=0) {

	if (sortByIndex == 0 or end-begin < sortByIndex) {
		sort(begin, end, AntiDiagonalSortOp<Tup>(genomeLength));
	}
	else {
		
		AntiDiagonalIndexSort<Tup> sorter;
		sorter.tuples=begin;
		sorter.length=genomeLength;
		vector<int> index(end-begin);
		std::iota(index.begin(), index.end(), 0);
		sort(index.begin(), index.end(), sorter);
		GenomePos pos;
		vector<pair<Tup, Tup> > temp(end-begin);
		copy(begin,end, temp.begin());
		
		for (int i=0; i < index.size(); i++) {
			temp[i]=*(begin+index[i]);
		}
		copy(temp.begin(), temp.end(), begin);
	}	

}


template<typename Tup>
void AntiDiagonalSort(vector<pair<Tup, Tup> > &vals, GenomePos length, int sortByIndex=0) {
	AntiDiagonalSort<Tup>(vals.begin(), vals.end(), length, sortByIndex);
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
void CartesianSort(typename vector<pair<Tup, Tup> >::iterator  begin, typename vector<pair<Tup, Tup> >::iterator  end) {
	sort(begin, end, CartesianSortOp<Tup>());
}

template<typename Tup>
void CartesianSort(vector<pair<Tup, Tup> > &vals, int s, int e) {
	sort(vals.begin() + s, vals.begin() + e, CartesianSortOp<Tup>());
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
void CartesianTargetSort(vector<pair<Tup, Tup>> &matches, int s, int e) {
	sort(matches.begin() + s, matches.begin() + e, CartesianTargetSortOp<Tup>());
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
