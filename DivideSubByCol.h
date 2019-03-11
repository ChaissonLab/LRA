#ifndef DIVIDE_SUB_BY_COL_H_
#define DIVIDE_SUB_BY_COL_H_


#include <iostream>
#include <utility>
#include <numeric> //std::floor
#include <cmath>
#include <iterator>
#include <map>

#include "SubProblem.h"
#include "Fragment_Info.h"
#include "Info.h"
#include "overload.h"
#include "Types.h"
#include "Point.h"

using std::cerr;
using std::cout;
using std::endl;
using std::iota;


// GetColInfo summarizes the col information in H2
void
GetColInfo (std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<info> & M) {

	unsigned int col = H1[H2[0]].se.second;
	unsigned int pstart = 0;
	unsigned int pend = 1;
	for (unsigned int i = 0; i < H2.size(); ++i) {
		//cerr << "i: " << i ;
		if (col == H1[H2[i]].se.second) {
			pend = i + 1;
		}
		else {
			//cerr << "create a info\n";
			info p(pstart, pend, col);
			M.push_back(p);
			pstart = i;
			pend = i + 1;
			col = H1[H2[i]].se.second;
		}

		if (i == H2.size() - 1) {
			//cerr << "create a info\n";
			info p(pstart, pend, col);
			M.push_back(p);
		}
	}
}

/*
// This function works for E array
void 
ScanPoints_Col(std::vector<info> & V, std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi, unsigned int & s, unsigned int & e) {


	std::map<long int, unsigned int> fmap; // <forward diagonal, number of points which have forward diagonal <= the current forward diagonal>
	for (unsigned int i = s; i < e; ++i) {
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {	
			long int l = static_cast<long int>(H1[H2[j]].second) - static_cast<long int>(H1[H2[j]].first);
			std::pair<std::map<long int, unsigned int>::iterator,bool> ret;
			ret = fmap.insert(std::pair<long int, unsigned int>(l, 1));
			if (ret.second == false) { // element with such forward diagonal l is already existed
				++(ret.first->second);
			}

		}
	}

	for (std::map<long int, unsigned int>::reverse_iterator it = fmap.rbegin(); it != fmap.rend(); ++it) {
		Bi.push_back(it->first);
	}	
}
*/


/*
// This function works for D array
void 
ScanPoints_Col(std::vector<info> & V, std::vector<Pair> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi, std::vector<unsigned int> & counter_D, unsigned int & s, unsigned int & e) {


	std::map<long int, unsigned int> fmap; // <forward diagonal, number of points which have forward diagonal <= the current forward diagonal>
	for (unsigned int i = s; i < e; ++i) {
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {	
			long int l = static_cast<long int>(H1[H2[j]].second) - static_cast<long int>(H1[H2[j]].first);
			std::pair<std::map<long int, unsigned int>::iterator,bool> ret;
			ret = fmap.insert(std::pair<long int, unsigned int>(l, 1));
			if (ret.second == false) { // element with such forward diagonal l is already existed
				++(ret.first->second);
			}

		}
	}

	for (std::map<long int, unsigned int>::reverse_iterator it = fmap.rbegin(); it != fmap.rend(); ++it) {
		Bi.push_back(it->first);
		counter_D.push_back(it->second);
	}	

	for (unsigned int p = 1; p < counter_D.size(); ++p) {
		counter_D[p] = counter_D[p - 1] + counter_D[p];
	}
}
*/


//ScanPoints_col find Di and Ei array
//Note this function didn't count the number of points which have forward diagonal <= the current forward diagonal>
void 
ScanPoints_Col(std::vector<info> & V, std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi, unsigned int & s, unsigned int & e,  bool & DE, unsigned int & n) {
	
	std::set<long int> ForwardIndex;
	for (unsigned int i = s; i < e; ++i) {
		unsigned int count = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {
			if (H1[H2[j]].ind == DE) {
				long int l = static_cast<long int>(H1[H2[j]].se.second) - static_cast<long int>(H1[H2[j]].se.first);
				ForwardIndex.insert(l);		
				++count;		
			}
		}

		if (count != 0) {
			if (DE == 1) V[i].SS_B.push_back(n);
			else V[i].SS_A.push_back(n);
		}
	}	

	for (std::set<long int>::reverse_iterator it = ForwardIndex.rbegin(); it != ForwardIndex.rend(); ++it) { // elements in D array are in the decreasing order
		Bi.push_back(*it);
	}
}



//ScanPoints_col find Di and Ei array
//Note this function didn't count the number of points which have forward diagonal <= the current forward diagonal>
void 
ScanPoints_Col(std::vector<info> & V, std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi,  std::vector<long int> & Ci, 
						unsigned int & s, unsigned int & e, unsigned int & n) {
	
	std::set<long int> ForwardIndex1;
	std::set<long int> ForwardIndex2;

	for (unsigned int i = s; i < e; ++i) {
		unsigned int count1 = 0;
		unsigned int count2 = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {
			if (H1[H2[j]].ind == 1) {
				long int l = static_cast<long int>(H1[H2[j]].se.second) - static_cast<long int>(H1[H2[j]].se.first);
				ForwardIndex1.insert(l);		
				++count1;		
			}
			else {
				long int r = static_cast<long int>(H1[H2[j]].se.second) - static_cast<long int>(H1[H2[j]].se.first);
				ForwardIndex2.insert(r);		
				++count2;		
			}
		}

		if (count1 != 0) {
			 V[i].SS_B.push_back(n);
		}
		if (count1 != 0 and count2 != 0) {
			++n;
			V[i].SS_B.push_back(n);
			V[i].SS_A.push_back(n);				
		}
	}	

	for (std::set<long int>::reverse_iterator it = ForwardIndex1.rbegin(); it != ForwardIndex1.rend(); ++it) { // elements in D array are in the decreasing order
		Bi.push_back(*it);
	}
	for (std::set<long int>::reverse_iterator it = ForwardIndex2.rbegin(); it != ForwardIndex2.rend(); ++it) { // elements in D array are in the decreasing order
		Ci.push_back(*it);
	}
}


void
Decide_Eb_Db_C (std::vector<long int> & Di, std::vector<long int> & Ei, std::vector<long int> & Db, std::vector<long int> & Eb, std::vector<unsigned int> & E) {

	for (unsigned int s = 0; s < Di.size(); ++s) {
		std::vector<unsigned int>::reverse_iterator t = Lower_Bound<std::vector<unsigned int>::reverse_iterator,long int>(E.rbegin(), E.rend(), Di[s], Ei);
		if (t == E.rbegin()) {
			break;
		}
		else{
			//std::prev(t);
			--t; 
			Db[s] = *t;
			Eb[*t] = s;
		}
	}

	unsigned int cur = -1;
	for (unsigned int s = 0; s < Eb.size(); ++s) {
		if (Eb[s] == -1 and cur == -1) {
			continue;
		}
		else if (Eb[s] != -1) {
			cur = Eb[s];
		}
		else {
			Eb[s] = cur;
		}
	}
}



void
DivideSubProbByCol (std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<info> & V, unsigned int start, unsigned int end, 
							unsigned int & n, std::vector<Subproblem> & Sub) { // [start, end) is a half open interval

	if (end == start + 1) { // subproblem A is empty, while B contains only one row. This is a leaf case.

		Subproblem s = Subproblem(n); // s is leaf subproblem
		Sub.push_back(s); 

		Subproblem ss = Subproblem(n + 1);
		Sub.push_back(ss); // ss is a subproblem which Di and Ei coming from one row
		unsigned int last = Sub.size();

		ScanPoints_Col(V, H1, H2, Sub[last - 2].Ei, Sub[last - 1].Di, start, end, n);	


		if (!Sub[last - 2].Ei.empty() and !Sub[last - 1].Di.empty()) { 
	
			// initialize Sub[last - 1]
			Sub[last - 1].Ei = Sub[last - 2].Ei;

			unsigned int l = Sub[last - 1].Di.size();
			unsigned int h = Sub[last - 1].Ei.size();

			std::vector<long int> p(h, -1);
			std::vector<long int> z(l, -1);
			std::vector<unsigned int> t(h, 0);
			Sub[last - 1].E = t;
			std::iota(Sub[last - 1].E.begin(), Sub[last - 1].E.end(), 0);
			Sub[last - 1].Eb = p;
			Sub[last - 1].Db = z;
			Decide_Eb_Db_C(Sub[last - 1].Di, Sub[last - 1].Ei, Sub[last - 1].Db, Sub[last - 1].Eb, Sub[last - 1].E);

			// initialize other attributes of this subproblem
			std::vector<float> v(l, 0);
			std::vector<unsigned int> w(l, 0);
			Sub[last - 1].Dv = v; 
			Sub[last - 1].Dp = w;
			Sub[last - 1].D = w;
			std::iota(Sub[last - 1].D.begin(), Sub[last - 1].D.end(), 0);

			std::vector<float> q(h, 0);
			Sub[last - 1].Ev = q;
			Sub[last - 1].Ep = t;
			std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
			Sub[last - 1].S_1.push(dummy_pair); 


			// initialize Sub[last - 2] -- the leaf case
			unsigned int hh = Sub[last - 2].Ei.size();
			std::vector<float> pp(hh, 0);
			std::vector<unsigned int> zz(hh, 0);
			Sub[last - 2].Ev = pp;
			Sub[last - 2].Ep = zz;
			Sub[last - 2].E = zz;
			std::iota(Sub[last - 2].E.begin(), Sub[last - 2].E.end(), 0);			
		}
		else if (!Sub[last - 2].Ei.empty()) {

			Sub.pop_back(); // delete subproblem ss

			// initialize Sub[last - 2] -- the leaf case
			unsigned int h = Sub[last - 2].Ei.size();
			std::vector<float> p(h, 0);
			std::vector<unsigned int> z(h, 0);
			Sub[last - 2].Ev = p;
			Sub[last - 2].Ep = z;
			Sub[last - 2].E = z;
			std::iota(Sub[last - 2].E.begin(), Sub[last - 2].E.end(), 0);				
		}
		else {
			Sub.pop_back(); // delete subproblem ss
			Sub.pop_back(); // delete subproblem s -- the leaf case
			--n;			
		}
	}
	else{

		Subproblem s = Subproblem(n);
		Sub.push_back(s);

		// scan the points to determine Di 
		unsigned int med = std::floor((start + end)/2);
		bool DE = 0; // DE == 0 means scan points to determin Di (find for end points); 
		ScanPoints_Col(V, H1, H2, (Sub.back()).Di, start, med, DE, n);

		// scan the points to determine Ei
		DE = 1;
		ScanPoints_Col(V, H1, H2, (Sub.back()).Ei, med, end, DE, n);

		if ((Sub.back()).Ei.empty() and (Sub.back()).Di.empty()) { // Di is empty and Ei is empty 
			Sub.pop_back();
			--n;
		}
		else if ((Sub.back()).Ei.empty() and !(Sub.back()).Di.empty()) { // Di is non-empty and Ei is empty
		}
		else if (!(Sub.back()).Ei.empty() and (Sub.back()).Di.empty()) { 

		}
		else { // non-leaf case

			// initialize (Sub.back()).Eb and (Sub.back()).Db
			unsigned int l = (Sub.back()).Di.size();
			unsigned int h = (Sub.back()).Ei.size();

			std::vector<long int> p(h, -1);
			std::vector<long int> z(l, -1);
			std::vector<unsigned int> t(h, 0);

			(Sub.back()).E = t;
			std::iota((Sub.back()).E.begin(), (Sub.back()).E.end(), 0);
			(Sub.back()).Eb = p;
			(Sub.back()).Db = z;
			Decide_Eb_Db_C((Sub.back()).Di, (Sub.back()).Ei, (Sub.back()).Db, (Sub.back()).Eb, (Sub.back()).E);

			// initialize other attributes of this subproblem
			std::vector<float> v(l, 0);
			std::vector<unsigned int> w(l, 0);
			(Sub.back()).Dv = v; 
			(Sub.back()).Dp = w;
			(Sub.back()).D = w;
			std::iota((Sub.back()).D.begin(), (Sub.back()).D.end(), 0);

			std::vector<float> q(h,0);
			(Sub.back()).Ev = q;
			(Sub.back()).Ep = t;
			std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
			(Sub.back()).S_1.push(dummy_pair); 
		}

		++n;
		DivideSubProbByCol(H1, H2, V, start, std::floor((start + end)/2), n, Sub);
		++n;
		DivideSubProbByCol(H1, H2, V, std::floor((start + end)/2), end, n, Sub);
	}
}

#endif