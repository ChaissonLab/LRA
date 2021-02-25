#ifndef DIVIDE_SUB_BY_COL2_H_
#define DIVIDE_SUB_BY_COL2_H_


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


//ScanPoints_col find Di and Ei array for non-leaf cases
// Note: H1[j].inv == 0 and backdiag
void 
ScanPoints_Col2(std::vector<info> & V, std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi, unsigned int & s, unsigned int & e,  bool & DE, unsigned int & n) {
	// elements in set are unique and follow an increasing order
	std::set<long int> ForwardIndex;
	for (unsigned int i = s; i < e; ++i) {
		unsigned int count = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {
			if (H1[H2[j]].ind == DE and H1[H2[j]].inv == 0) {
				long int l = (long int)(H1[H2[j]].se.second) + (long int)(H1[H2[j]].se.first);
				ForwardIndex.insert(l);		
				++count;		
			}
		}

		if (count != 0) {
			if (DE == 1) V[i].SS_B2.push_back(n);
			else V[i].SS_A2.push_back(n);
		}
	}	
	// elements in D/E array are in the ascending order
	for (std::set<long int>::iterator it = ForwardIndex.begin(); it != ForwardIndex.end(); ++it) { 
		Bi.push_back(*it);
	}
}


//ScanPoints_col find Di and Ei array for leaf cases
// Note: H1[j].inv == 0 and backdiag
void 
ScanPoints_Col2(std::vector<info> & V, std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi,  std::vector<long int> & Ci, 
						unsigned int & s, unsigned int & e, unsigned int & n) {
	
	std::set<long int> ForwardIndex1;
	std::set<long int> ForwardIndex2;

	for (unsigned int i = s; i < e; ++i) {
		unsigned int count1 = 0;
		unsigned int count2 = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {
			if (H1[H2[j]].ind == 1 and H1[H2[j]].inv == 0) {
				long int l = (long int)(H1[H2[j]].se.second) + (long int)(H1[H2[j]].se.first);
				ForwardIndex1.insert(l);		
				++count1;		
			}
			else if (H1[H2[j]].ind == 0 and H1[H2[j]].inv == 0) {
				long int r = (long int)(H1[H2[j]].se.second) + (long int)(H1[H2[j]].se.first);
				ForwardIndex2.insert(r);		
				++count2;		
			}
		}

		if (count1 != 0 and count2 != 0) {
			V[i].SS_B2.push_back(n);
			V[i].SS_A2.push_back(n);				
		}
	}	

	for (std::set<long int>::iterator it = ForwardIndex1.begin(); it != ForwardIndex1.end(); ++it) { // elements in D array are in the ascending order
		Bi.push_back(*it);
	}
	for (std::set<long int>::iterator it = ForwardIndex2.begin(); it != ForwardIndex2.end(); ++it) { // elements in D array are in the ascending order
		Ci.push_back(*it);
	}
}


void
Decide_Eb_Db_C2 (std::vector<long int> & Di, std::vector<long int> & Ei, std::vector<long int> & Db, std::vector<long int> & Eb, std::vector<unsigned int> & E) {

	for (unsigned int s = 0; s < Di.size(); ++s) {
		// find the index *t that Ei[*t] is the first element which is >= Di[s]
		std::vector<unsigned int>::iterator t = Lower_Bound<std::vector<unsigned int>::iterator,long int>(E.begin(), E.end(), Di[s], Ei); 
		if (t == E.end()) {
			break;
		}
		else{
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
DivideSubProbByCol2 (std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<info> & V, unsigned int start, unsigned int end, 
							unsigned int & n, StackOfSubProblems & Sub, int & eeC) { // [start, end) is a half open interval

	if (end == start + 1) { // subproblem A is empty, while B contains only one row. This is a leaf case.

		Subproblem ss = Subproblem(n);
		Sub.Push_Back(eeC, ss); // ss is a subproblem which Di and Ei coming from one row
		++eeC;
		ScanPoints_Col2(V, H1, H2, Sub[eeC - 1].Ei, Sub[eeC - 1].Di, start, end, n);	

		if (!Sub[eeC - 1].Ei.empty() and !Sub[eeC - 1].Di.empty()) { 
	
			// initialize Sub[eeC- 1]
			unsigned int l = Sub[eeC - 1].Di.size();
			unsigned int h = Sub[eeC - 1].Ei.size();

			Sub[eeC - 1].E.assign(h, 0); 
			std::iota(Sub[eeC - 1].E.begin(), Sub[eeC - 1].E.end(), 0);
			Sub[eeC - 1].Eb.assign(h, -1);
			Sub[eeC - 1].Db.assign(l, -1); 
			Decide_Eb_Db_C2(Sub[eeC - 1].Di, Sub[eeC - 1].Ei, Sub[eeC - 1].Db, Sub[eeC - 1].Eb, Sub[eeC - 1].E);

			// initialize other attributes of this subproblem
			Sub[eeC - 1].Dv.assign(l, 0);
			Sub[eeC - 1].Dp.assign(l, 0);
			Sub[eeC - 1].D.assign(l, 0);
			std::iota(Sub[eeC - 1].D.begin(), Sub[eeC - 1].D.end(), 0);

			Sub[eeC - 1].Ev.assign(h, 0);
			Sub[eeC - 1].Ep.assign(h, 0); 
			std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
			Sub[eeC - 1].S_1.push(dummy_pair); 
		
		}
		else {
			Sub.pop_back(); // delete subproblem ss
			// Sub.ClearSingle(eeC);
			--eeC;
			--n;			
		}
	}
	else{

		Subproblem s = Subproblem(n);
		Sub.Push_Back(eeC, s);
		++eeC;

		// scan the points to determine Di 
		unsigned int med = std::floor((start + end)/2);
		bool DE = 1; // DE == 1 means scan points to determin Ei (find for start points); 
		ScanPoints_Col2(V, H1, H2, Sub[eeC-1].Ei, start, med, DE, n);

		DE = 0; // scan the points to determine Di
		ScanPoints_Col2(V, H1, H2, Sub[eeC-1].Di, med, end, DE, n);

		if (Sub[eeC-1].Ei.empty() and Sub[eeC-1].Di.empty()) { // Di is empty and Ei is empty 
			Sub.pop_back(); // delete subproblem ss
			// Sub.ClearSingle(eeC);
			--eeC;
			--n;
		}
		else if (Sub[eeC-1].Ei.empty() and !Sub[eeC-1].Di.empty()) { // Di is non-empty and Ei is empty
			++n;
			DivideSubProbByCol2(H1, H2, V, std::floor((start + end)/2), end, n, Sub, eeC);
		}
		else if (!Sub[eeC-1].Ei.empty() and Sub[eeC-1].Di.empty()) { 
			++n;
			DivideSubProbByCol2(H1, H2, V, start, std::floor((start + end)/2), n, Sub, eeC);			
		}
		else { // non-leaf case

			// initialize Sub[eeC-1].Eb and Sub[eeC-1].Db
			unsigned int l = Sub[eeC-1].Di.size();
			unsigned int h = Sub[eeC-1].Ei.size();

			//std::vector<long int> p(h, -1);
			//std::vector<long int> z(l, -1);
			//std::vector<unsigned int> t(h, 0);

			Sub[eeC-1].E.assign(h, 0); 
			std::iota(Sub[eeC-1].E.begin(), Sub[eeC-1].E.end(), 0);
			Sub[eeC-1].Eb.assign(h, -1); 
			Sub[eeC-1].Db.assign(l, -1);
			Decide_Eb_Db_C2(Sub[eeC-1].Di, Sub[eeC-1].Ei, Sub[eeC-1].Db, Sub[eeC-1].Eb, Sub[eeC-1].E);

			// initialize other attributes of this subproblem
			//std::vector<float> v(l, 0);
			//std::vector<unsigned int> w(l, 0);
			Sub[eeC-1].Dv.assign(l, 0); 
			Sub[eeC-1].Dp.assign(l, 0); 
			Sub[eeC-1].D.assign(l, 0);
			std::iota(Sub[eeC-1].D.begin(), Sub[eeC-1].D.end(), 0);

			//std::vector<float> q(h,0);
			Sub[eeC-1].Ev.assign(h, 0);
			Sub[eeC-1].Ep.assign(h, 0);
			std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
			Sub[eeC-1].S_1.push(dummy_pair); 
			++n;
			DivideSubProbByCol2(H1, H2, V, std::floor((start + end)/2), end, n, Sub, eeC);
			++n;
			DivideSubProbByCol2(H1, H2, V, start, std::floor((start + end)/2), n, Sub, eeC);
		}
	}
}

#endif