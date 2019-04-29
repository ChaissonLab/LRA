#ifndef DIVIDE_SUB_BY_COL1_H_
#define DIVIDE_SUB_BY_COL1_H_


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
ScanPoints_Col1(std::vector<info> & V, std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi, unsigned int & s, unsigned int & e) {


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
ScanPoints_Col1(std::vector<info> & V, std::vector<Pair> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi, std::vector<unsigned int> & counter_D, unsigned int & s, unsigned int & e) {


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


//ScanPoints_col find Di and Ei array for non-leaf cases
//Note this function didn't count the number of points which have forward diagonal <= the current forward diagonal>
void 
ScanPoints_Col1(std::vector<info> & V, std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi, unsigned int & s, unsigned int & e,  bool & DE, unsigned int & n) {
	
	std::set<long int> ForwardIndex;
	for (unsigned int i = s; i < e; ++i) {
		unsigned int count = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {
			if (H1[H2[j]].ind == DE and H1[H2[j]].inv == 1) {
				long int l = static_cast<long int>(H1[H2[j]].se.second) - static_cast<long int>(H1[H2[j]].se.first);
				ForwardIndex.insert(l);		
				++count;		
			}
		}

		if (count != 0) {
			if (DE == 1) V[i].SS_B1.push_back(n);
			else V[i].SS_A1.push_back(n);
		}
	}	

	for (std::set<long int>::reverse_iterator it = ForwardIndex.rbegin(); it != ForwardIndex.rend(); ++it) { // elements in D array are in the decreasing order
		Bi.push_back(*it);
	}
}



//ScanPoints_col find Di and Ei array for leaf cases
//Note this function didn't count the number of points which have forward diagonal <= the current forward diagonal>
void 
ScanPoints_Col1(std::vector<info> & V, std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<long int> & Bi,  std::vector<long int> & Ci, 
						unsigned int & s, unsigned int & e, unsigned int & n) {
	
	std::set<long int> ForwardIndex1;
	std::set<long int> ForwardIndex2;

	for (unsigned int i = s; i < e; ++i) {
		unsigned int count1 = 0;
		unsigned int count2 = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {
			if (H1[H2[j]].ind == 1 and H1[H2[j]].inv == 1) {
				long int l = static_cast<long int>(H1[H2[j]].se.second) - static_cast<long int>(H1[H2[j]].se.first);
				ForwardIndex1.insert(l);		
				++count1;		
			}
			else if (H1[H2[j]].ind == 0 and H1[H2[j]].inv == 1) {
				long int r = static_cast<long int>(H1[H2[j]].se.second) - static_cast<long int>(H1[H2[j]].se.first);
				ForwardIndex2.insert(r);		
				++count2;		
			}
		}

		if (count1 != 0 and count2 != 0) {
			V[i].SS_B1.push_back(n);
			V[i].SS_A1.push_back(n);				
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
Decide_Eb_Db_C1 (std::vector<long int> & Di, std::vector<long int> & Ei, std::vector<long int> & Db, std::vector<long int> & Eb, std::vector<unsigned int> & E) {

	for (unsigned int s = 0; s < Di.size(); ++s) {
		std::vector<unsigned int>::reverse_iterator t = Lower_Bound<std::vector<unsigned int>::reverse_iterator,long int>(E.rbegin(), E.rend(), Di[s], Ei); // find the index *t that Ei[*t] 
																																							// is the first element which is >= Di[s]
																																							// Note: here we compare from the right
		if (t == E.rbegin()) {
			break;
		}
		else{
			//std::prev(t);
			--t; // move to right by one step
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
DivideSubProbByCol1 (std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<info> & V, unsigned int start, unsigned int end, 
							unsigned int & n, StackOfSubProblems & Sub, int & eeC) { // [start, end) is a half open interval

	if (end == start + 1) { // subproblem A is empty, while B contains only one row. This is a leaf case.

		Subproblem ss = Subproblem(n);
		Sub.Push_Back(eeC, ss); // ss is a subproblem which Di and Ei coming from one row
		++eeC;
		ScanPoints_Col1(V, H1, H2, Sub[eeC - 1].Ei, Sub[eeC - 1].Di, start, end, n);	

		if (!Sub[eeC - 1].Ei.empty() and !Sub[eeC - 1].Di.empty()) { 
	
			// initialize Sub[eeC- 1]
			unsigned int l = Sub[eeC - 1].Di.size();
			unsigned int h = Sub[eeC - 1].Ei.size();

			Sub[eeC - 1].E.assign(h, 0); 
			std::iota(Sub[eeC - 1].E.begin(), Sub[eeC - 1].E.end(), 0);
			Sub[eeC - 1].Eb.assign(h, -1);
			Sub[eeC - 1].Db.assign(l, -1); 
			Decide_Eb_Db_C1(Sub[eeC - 1].Di, Sub[eeC - 1].Ei, Sub[eeC - 1].Db, Sub[eeC - 1].Eb, Sub[eeC - 1].E);

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
			//Sub.pop_back(); // delete subproblem ss
			Sub.ClearSingle(eeC);
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
		bool DE = 0; // DE == 0 means scan points to determin Di (find for end points); 
		ScanPoints_Col1(V, H1, H2, Sub[eeC-1].Di, start, med, DE, n);

		// scan the points to determine Ei
		DE = 1;
		ScanPoints_Col1(V, H1, H2, Sub[eeC-1].Ei, med, end, DE, n);

		if (Sub[eeC-1].Ei.empty() and Sub[eeC-1].Di.empty()) { // Di is empty and Ei is empty 
			Sub.ClearSingle(eeC);
			--eeC;
			--n;
		}
		else if (Sub[eeC-1].Ei.empty() and !Sub[eeC-1].Di.empty()) { // Di is non-empty and Ei is empty
			++n;
			DivideSubProbByCol1(H1, H2, V, start, std::floor((start + end)/2), n, Sub, eeC);			
		}
		else if (!Sub[eeC-1].Ei.empty() and Sub[eeC-1].Di.empty()) { 
			++n;
			DivideSubProbByCol1(H1, H2, V, std::floor((start + end)/2), end, n, Sub, eeC);
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
			Decide_Eb_Db_C1(Sub[eeC-1].Di, Sub[eeC-1].Ei, Sub[eeC-1].Db, Sub[eeC-1].Eb, Sub[eeC-1].E);

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
			DivideSubProbByCol1(H1, H2, V, start, std::floor((start + end)/2), n, Sub, eeC);
			++n;
			DivideSubProbByCol1(H1, H2, V, std::floor((start + end)/2), end, n, Sub, eeC);
		}
	}
}



/*
void
DivideSubProbByCol1 (std::vector<Point> & H1, std::vector<unsigned int> & H2, std::vector<info> & V, unsigned int start, unsigned int end, 
							unsigned int & n, StackOfSubProblems & Sub, int & eeC) { // [start, end) is a half open interval

	if (end == start + 1) { // subproblem A is empty, while B contains only one row. This is a leaf case.

		Subproblem ss = Subproblem(n);
		Sub.Push_Back(eeC, ss); // ss is a subproblem which Di and Ei coming from one row
		++eeC;

		// scan the points to determine Ei and Di
		ScanPoints_Col1(V, H1, H2, Sub[eeC - 1].Ei, Sub[eeC - 1].Di, start, end, n);	


		if (!Sub[eeC - 1].Ei.empty() and !Sub[eeC - 1].Di.empty()) { 
	
			// initialize Sub[eeC - 1]
			unsigned int l = Sub[eeC - 1].Di.size();
			unsigned int h = Sub[eeC - 1].Ei.size();

			Sub[eeC - 1].E.assign(h, 0);
			std::iota(Sub[eeC - 1].E.begin(), Sub[eeC - 1].E.end(), 0);
			Sub[eeC - 1].Eb.assign(h, -1);
			Sub[eeC - 1].Db.assign(l, -1);
			Decide_Eb_Db_C1(Sub[eeC - 1].Di, Sub[eeC - 1].Ei, Sub[eeC - 1].Db, Sub[eeC - 1].Eb, Sub[eeC - 1].E);

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
			// delete subproblem ss			
			Sub.ClearSingle(eeC);
			--eeC;
			//--n;		
		}
	}
	else{

		Subproblem s = Subproblem(n);
		Sub.Push_Back(eeC, s);
		++eeC;

		// scan the points to determine Di 
		unsigned int med = std::floor((start + end)/2);
		bool DE = 0; // DE == 0 means scan points to determin Di (find for end points); 
		ScanPoints_Col1(V, H1, H2, Sub[eeC - 1].Di, start, med, DE, n);

		// scan the points to determine Ei
		DE = 1;
		ScanPoints_Col1(V, H1, H2, Sub[eeC - 1].Ei, med, end, DE, n);

		if (Sub[eeC - 1].Ei.empty() and Sub[eeC - 1].Di.empty()) { // Di is empty and Ei is empty 
			Sub.ClearSingle(eeC);
			--eeC;
			//--n;
		}
		else if (Sub[eeC - 1].Ei.empty() and !Sub[eeC - 1].Di.empty()) { // Di is non-empty and Ei is empty
			Sub.ClearSingle(eeC);
			--eeC;
			DivideSubProbByCol1(H1, H2, V, start, std::floor((start + end)/2), n, Sub, eeC);
		}
		else if (!Sub[eeC - 1].Ei.empty() and Sub[eeC - 1].Di.empty()) { 
			Sub.ClearSingle(eeC);
			--eeC;
			DivideSubProbByCol1(H1, H2, V, std::floor((start + end)/2), end, n, Sub, eeC);
		}
		else { // non-leaf case

			// initialize Sub[eeC - 1].Eb and Sub[eeC - 1].Db
			unsigned int l = Sub[eeC - 1].Di.size();
			unsigned int h = Sub[eeC - 1].Ei.size();


			Sub[eeC - 1].E.assign(h, 0);
			std::iota(Sub[eeC - 1].E.begin(), Sub[eeC - 1].E.end(), 0);
			Sub[eeC - 1].Eb.assign(h, -1);
			Sub[eeC - 1].Db.assign(l, -1);
			Decide_Eb_Db_C1(Sub[eeC - 1].Di, Sub[eeC - 1].Ei, Sub[eeC - 1].Db, Sub[eeC - 1].Eb, Sub[eeC - 1].E);

			// initialize other attributes of this subproblem
			Sub[eeC - 1].Dv.assign(l, 0);
			Sub[eeC - 1].Dp.assign(l, 0);
			Sub[eeC - 1].D.assign(l, 0); 
			std::iota(Sub[eeC - 1].D.begin(), Sub[eeC - 1].D.end(), 0);

			Sub[eeC - 1].Ev.assign(h, 0); 
			Sub[eeC - 1].Ep.assign(h, 0); 
			std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
			Sub[eeC - 1].S_1.push(dummy_pair); 

			++n;
			DivideSubProbByCol1(H1, H2, V, start, std::floor((start + end)/2), n, Sub, eeC);
			++n;
			DivideSubProbByCol1(H1, H2, V, std::floor((start + end)/2), end, n, Sub, eeC);

		}


	}
}
*/
#endif