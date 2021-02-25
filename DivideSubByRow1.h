#ifndef DIVIDE_SUB_BY_ROW1_H_
#define DIVIDE_SUB_BY_ROW1_H_


#include <iostream>
#include <utility> // std::pair
#include <numeric> //std::floor std::iota
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


// GetRowInfo summarize the row information of H1 in M
void
GetRowInfo (std::vector<Point> & H1, std::vector<info> & M) {

	unsigned int row = H1[0].se.first;
	unsigned int pstart = 0;
	unsigned int pend = 1;
	for (unsigned int i = 0; i < H1.size(); ++i) {
		//cerr << "i: " << i ;
		if (row == H1[i].se.first) {
			//cerr << "row: " << row << "H1[i]..se.first: " << H1[i].se.first << endl;
			pend = i + 1;
		}
		else {
			//cerr << "create a info\n";
			info p(pstart, pend, row);
			M.push_back(p);
			pstart = i;
			pend = i + 1;
			row = H1[i].se.first;
		}

		if (i == H1.size() - 1) {
			//cerr << "create a info\n";
			info p(pstart, pend, row);
			M.push_back(p);
		}
	}
}


// This function finds Di and Ei array
// Note this function didn't count the number of points which have forward diagonal <= the current forward diagonal
void 
ScanPoints_Row1 (std::vector<info> & V, std::vector<Point> & H1, std::vector<long int> & Bi,  unsigned int & s, unsigned int & e, bool & DE, unsigned int & n) {

	// elements in set are unique and follow an increasing order
	std::set<long int> ForwardIndex;
	for (unsigned int i = s; i < e; ++i) {

		unsigned int count = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {
			if (H1[j].ind == DE and H1[j].inv == 1) { // H1[j].ind == DE == 1 means finding start points
				long int l = (long int)(H1[j].se.second) - (long int)(H1[j].se.first);
				ForwardIndex.insert(l);
				++count;				
			}
		}

		if (count != 0) {
			if (DE == 1) V[i].SS_B1.push_back(n);
			else V[i].SS_A1.push_back(n);
		}

	}	

	for (std::set<long int>::iterator it = ForwardIndex.begin(); it != ForwardIndex.end(); ++it) {
		Bi.push_back(*it);
	}
}



// This function finds Di and Ei array for leaf-case
// Note this function didn't count the number of points which have forward diagonal <= the current forward diagonal
void 
ScanPoints_Row1 (std::vector<info> & V, std::vector<Point> & H1, std::vector<long int> & Bi, std::vector<long int> & Ci,
						 unsigned int & s, unsigned int & e, unsigned int & n) {

	std::set<long int> ForwardIndex1; // ForwardIndex1 is for Ei array
	std::set<long int> ForwardIndex2; // ForwardIndex2 is for Di array
	for (unsigned int i = s; i < e; ++i) {

		unsigned int count1 = 0;
		unsigned int count2 = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) { 
			if (H1[j].ind == 1 and H1[j].inv == 1) { // H1[j].ind  == 1 means finding start points
				long int l = (long int)(H1[j].se.second) - (long int)(H1[j].se.first);
				ForwardIndex1.insert(l);
				++count1;				  
			}
			else if (H1[j].ind == 0 and H1[j].inv == 1) { // H1[j].ind  == 0 means finding end points
				long int r = (long int)(H1[j].se.second) - (long int)(H1[j].se.first);
				ForwardIndex2.insert(r);
				++count2;				
			}			
		}

		if (count1 != 0 and count2 != 0) {
			V[i].SS_B1.push_back(n);
			V[i].SS_A1.push_back(n);		
		}

	}	

	for (std::set<long int>::iterator it = ForwardIndex1.begin(); it != ForwardIndex1.end(); ++it) {
		Bi.push_back(*it);
	}
	for (std::set<long int>::iterator it = ForwardIndex2.begin(); it != ForwardIndex2.end(); ++it) {
		Ci.push_back(*it);
	}
}




/*
// This function finds Di and Ei array
// Note: this function also count the number of points which have forward diagonal <= the current forward diagonal> 
void 
ScanPoints_Row1 (std::vector<info> & V, std::vector<Point> & H1, std::vector<long int> & Bi,  unsigned int & s, unsigned int & e) {

	std::map<long int, unsigned int> fmap; // <forward diagonal, number of points which have forward diagonal <= the current forward diagonal>
	for (unsigned int i = s; i < e; ++i) {
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {	
			long int l = (long int))(H1[j].second) - (long int))(H1[j].first);
			std::pair<std::map<long int, unsigned int>::iterator, bool> ret;
			ret = fmap.insert(std::pair<long int, unsigned int>(l, 1));
			if (ret.second == false) { // element with such forward diagonal l is already existed
				++ret.first->second;
			}

		}
	}
	for (std::map<long int, unsigned int>::iterator it = fmap.begin(); it != fmap.end(); ++it) {
		Bi.push_back(it->first);
	}	
}
*/



/*
// This function works for Di array
// Note: this function also count the number of points which have forward diagonal <= the current forward diagonal> 
// Note: requires input of counter_d array
void 
ScanPoints_Row1 (std::vector<info> & V, std::vector<Pair> & H1, std::vector<long int> & Bi, std::vector<unsigned int> & counter_D, unsigned int & s, unsigned int & e) {

	std::map<long int, unsigned int> fmap; // <forward diagonal, number of points which have forward diagonal <= the current forward diagonal>
	for (unsigned int i = s; i < e; ++i) {
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {	
			long int l = (long int)(H1[j].second) - (long int)(H1[j].first);
			std::pair<std::map<long int, unsigned int>::iterator, bool> ret;
			ret = fmap.insert(std::pair<long int, unsigned int>(l, 1));
			if (ret.second == false) { // element with such forward diagonal l is already existed
				++ret.first->second;
			}

		}
	}
	for (std::map<long int, unsigned int>::iterator it = fmap.begin(); it != fmap.end(); ++it) {
		Bi.push_back(it->first);
		counter_D.push_back(it->second);
	}	
	for (unsigned int p = 1; p < counter_D.size(); ++p) {
		counter_D[p] = counter_D[p - 1] + counter_D[p];
	}
}
*/


void
Decide_Eb_Db_R1 (std::vector<long int> & Di, std::vector<long int> & Ei, std::vector<long int> & Db, std::vector<long int> & Eb, std::vector<unsigned int> & E) {

	for (unsigned int s = 0; s < Di.size(); ++s) {
		// find the index *t that Ei[*t] is the first element which is >= Di[s]
		//
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


/*
void
DivideSubProbByRow1 (std::vector<Point> & H1, std::vector<info> & V, unsigned int start, unsigned int end, 
							unsigned int & n, StackOfSubProblems & Sub, int & eeR) { // [start, end) is a half open interval

	if (end == start + 1) { // subproblem A is empty, while B contains only one row. This is a leaf case.


		Subproblem s = Subproblem(n); // s is leaf subproblem
		Sub.Push_Back(eeR, s);
		++eeR;

		Subproblem ss = Subproblem(n + 1);
		Sub.Push_Back(eeR, ss); // ss is a subproblem which Di and Ei coming from one row
		++eeR;
		unsigned int last = eeR;

		// scan the points to determine Ei and Di
		ScanPoints_Row1(V, H1, Sub[last - 2].Ei, Sub[last - 1].Di, start, end, n);	

		if (!Sub[last - 2].Ei.empty() and !Sub[last - 1].Di.empty()) { 
			// initialize Sub[last - 1]
			Sub[last - 1].Ei = Sub[last - 2].Ei;

			unsigned int l = Sub[last - 1].Di.size();
			unsigned int h = Sub[last - 1].Ei.size();

			//std::vector<long int> p(h, -1);
			//std::vector<long int> z(l, -1);
			//std::vector<unsigned int> t(h, 0);
			Sub[last - 1].E.assign(h, 0);
			std::iota(Sub[last - 1].E.begin(), Sub[last - 1].E.end(), 0);
			Sub[last - 1].Eb.assign(h, -1);
			Sub[last - 1].Db.assign(l, -1);
			Decide_Eb_Db_R1(Sub[last - 1].Di, Sub[last - 1].Ei, Sub[last - 1].Db, Sub[last - 1].Eb, Sub[last - 1].E);

			// initialize other attributes of this subproblem
			//std::vector<float> v(l, 0);
			//std::vector<unsigned int> w(l, 0);
			Sub[last - 1].Dv.assign(l, 0); 
			Sub[last - 1].Dp.assign(l, 0);
			Sub[last - 1].D.assign(l, 0);
			std::iota(Sub[last - 1].D.begin(), Sub[last - 1].D.end(), 0);

			//std::vector<float> q(h, 0);
			Sub[last - 1].Ev.assign(h, 0);
			Sub[last - 1].Ep.assign(h, 0);
			std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
			Sub[last - 1].S_1.push(dummy_pair); 


			// initialize Sub[last - 2] -- the leaf case
			unsigned int hh = Sub[last - 2].Ei.size();
			//std::vector<float> pp(hh, 0);
			//std::vector<unsigned int> zz(hh, 0);
			Sub[last - 2].Ev.assign(hh, 0);
			Sub[last - 2].Ep.assign(hh, 0); 
			Sub[last - 2].E.assign(hh, 0);
			std::iota(Sub[last - 2].E.begin(), Sub[last - 2].E.end(), 0);		

		}
		else if (!Sub[last - 2].Ei.empty()) {
			//Sub.pop_back(); // delete subproblem ss
			Sub.ClearSingle(eeR);
			--eeR;		
			// initialize Sub[last - 2] -- the leaf case
			unsigned int h = Sub[last - 2].Ei.size();
			//std::vector<float> p(h, 0);
			//std::vector<unsigned int> z(h, 0);
			Sub[last - 2].Ev.assign(h, 0); 
			Sub[last - 2].Ep.assign(h, 0); 
			Sub[last - 2].E.assign(h, 0);
			std::iota(Sub[last - 2].E.begin(), Sub[last - 2].E.end(), 0);	
		}
		else {
			//Sub.pop_back(); // delete subproblem ss
			Sub.ClearSingle(eeR);
			--eeR;		
			Sub.ClearSingle(eeR);
			--eeR;		
			//Sub.pop_back(); // delete subproblem s -- the leaf case
			--n;			
		}
	}
	else{

		Subproblem s = Subproblem(n);
		Sub.Push_Back(eeR, s);
		++eeR;
		// scan the points to determine Di 
		unsigned int med = std::floor((start + end)/2);
		bool DE = 0; // DE == 0 means scan points to determin Di (find for end points); 
		//cerr << "scan points to determin Di in ["<< start << ", " << med << ")" << endl;
		ScanPoints_Row1(V, H1, Sub[eeR -1].Di, start, med, DE, n);
		// scan the points to determine Ei
		//cerr << "scan points to determine Ei in ["<< med << ", " << end << ")" << endl;
		DE = 1;
		ScanPoints_Row1(V, H1, Sub[eeR -1].Ei, med, end, DE, n);


		if (Sub[eeR -1].Ei.empty() and Sub[eeR -1].Di.empty()) { // Di is empty and Ei is empty  
			//Sub.pop_back();
			Sub.ClearSingle(eeR);
			--eeR;	
			--n;
		}
		else if (Sub[eeR -1].Ei.empty() and !Sub[eeR -1].Di.empty()) { // Di is non-empty and Ei is empty
			//cerr << "Di is non-empty and Ei is empty: " << n << "\n";
		}
		else if (!Sub[eeR -1].Ei.empty() and Sub[eeR -1].Di.empty()) { // Di is empty and Ei is non-empty 
			//cerr << "Di is empty and Ei is non-empty: " << n << "\n";

		}
		else { 

			// This is an non-leaf case
			// initialize Sub[eeR -1].Eb and Sub[eeR -1].Db
			unsigned int l = Sub[eeR -1].Di.size();
			unsigned int h = Sub[eeR -1].Ei.size();

			//std::vector<long int> p(h, -1);
			//std::vector<long int> z(l, -1);
			//std::vector<unsigned int> t(h, 0);
			Sub[eeR -1].E.assign(h, 0);
			std::iota(Sub[eeR -1].E.begin(), Sub[eeR -1].E.end(), 0);
			Sub[eeR -1].Eb.assign(h, -1); 
			Sub[eeR -1].Db.assign(l, -1);
			Decide_Eb_Db_R1(Sub[eeR -1].Di, Sub[eeR -1].Ei, Sub[eeR -1].Db, Sub[eeR -1].Eb, Sub[eeR -1].E);

			// initialize other attributes of this subproblem
			//std::vector<float> v(l, 0);
			//std::vector<unsigned int> w(l, 0);
			Sub[eeR -1].Dv.assign(l, 0); 
			Sub[eeR -1].Dp.assign(l, 0);
			Sub[eeR -1].D.assign(l, 0);
			std::iota(Sub[eeR -1].D.begin(), Sub[eeR -1].D.end(), 0);

			//std::vector<float> q(h, 0);
			Sub[eeR -1].Ev.assign(h, 0);
			Sub[eeR -1].Ep.assign(h, 0);
			std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
			Sub[eeR -1].S_1.push(dummy_pair); 
		}
		++n;
		//cerr <<"start: " << start+ 1 << ", med: " <<  std::floor((start + 1 + end + 1)/2) << ", n: " <<  n << "\n";
		DivideSubProbByRow1(H1, V, start, std::floor((start + end)/2), n, Sub, eeR);
		++n;
		//cerr <<"med: " << std::floor((start + 1 + end + 1)/2) << ", end: " << end + 1  << ", n: " <<  n << "\n";
		DivideSubProbByRow1(H1, V, std::floor((start + end)/2), end, n, Sub, eeR);
	}
}
*/




void
DivideSubProbByRow1 (std::vector<Point> & H1, std::vector<info> & V, unsigned int start, unsigned int end, 
							unsigned int & n, StackOfSubProblems & Sub, int & eeR) { // [start, end) is a half open interval

	if (end == start + 1) { // subproblem A is empty, while B contains only one row. This is a leaf case.


		Subproblem ss = Subproblem(n);
		Sub.Push_Back(eeR, ss); // ss is a subproblem which Di and Ei coming from one row
		++eeR;

		// scan the points to determine Ei and Di
		ScanPoints_Row1(V, H1, Sub[eeR - 1].Ei, Sub[eeR - 1].Di, start, end, n);	

		if (!Sub[eeR - 1].Ei.empty() and !Sub[eeR - 1].Di.empty()) { 
			
			// initialize Sub[eeR - 1]
			unsigned int l = Sub[eeR - 1].Di.size();
			unsigned int h = Sub[eeR - 1].Ei.size();

			Sub[eeR - 1].E.assign(h, 0);
			std::iota(Sub[eeR - 1].E.begin(), Sub[eeR - 1].E.end(), 0);
			Sub[eeR - 1].Eb.assign(h, -1);
			Sub[eeR - 1].Db.assign(l, -1);
			Decide_Eb_Db_R1(Sub[eeR - 1].Di, Sub[eeR - 1].Ei, Sub[eeR - 1].Db, Sub[eeR - 1].Eb, Sub[eeR - 1].E);

			// initialize other attributes of this subproblem
			Sub[eeR - 1].Dv.assign(l, 0); 
			Sub[eeR- 1].Dp.assign(l, 0);
			Sub[eeR - 1].D.assign(l, 0);
			std::iota(Sub[eeR - 1].D.begin(), Sub[eeR - 1].D.end(), 0);

			Sub[eeR - 1].Ev.assign(h, 0);
			Sub[eeR - 1].Ep.assign(h, 0);
			std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
			Sub[eeR - 1].S_1.push(dummy_pair); 
	
		}
		else {
			Sub.pop_back(); // delete subproblem ss
			// Sub.ClearSingle(eeR);
			--eeR;		
			--n;			
		}
	}
	else{

		Subproblem s = Subproblem(n);
		Sub.Push_Back(eeR, s);
		++eeR;
		// scan the points to determine Di 
		unsigned int med = std::floor((start + end)/2);
		bool DE = 0; // DE == 0 means scan points to determin Di (find for end points); 
		//cerr << "scan points to determin Di in ["<< start << ", " << med << ")" << endl;
		ScanPoints_Row1(V, H1, Sub[eeR -1].Di, start, med, DE, n);
		// scan the points to determine Ei
		//cerr << "scan points to determine Ei in ["<< med << ", " << end << ")" << endl;
		DE = 1;
		ScanPoints_Row1(V, H1, Sub[eeR -1].Ei, med, end, DE, n);


		if (Sub[eeR -1].Ei.empty() and Sub[eeR -1].Di.empty()) { // Di is empty and Ei is empty  
			Sub.pop_back();
			// Sub.ClearSingle(eeR);
			--eeR;	
			--n;
		}
		else if (Sub[eeR -1].Ei.empty() and !Sub[eeR -1].Di.empty()) { // Di is non-empty and Ei is empty
			//cerr << "Di is non-empty and Ei is empty: " << n << "\n";
			++n;
			//cerr <<"start: " << start+ 1 << ", med: " <<  std::floor((start + 1 + end + 1)/2) << ", n: " <<  n << "\n";
			DivideSubProbByRow1(H1, V, start, std::floor((start + end)/2), n, Sub, eeR);
		}
		else if (!Sub[eeR -1].Ei.empty() and Sub[eeR -1].Di.empty()) { // Di is empty and Ei is non-empty 
			//cerr << "Di is empty and Ei is non-empty: " << n << "\n";
			++n;
			//cerr <<"med: " << std::floor((start + 1 + end + 1)/2) << ", end: " << end + 1  << ", n: " <<  n << "\n";
			DivideSubProbByRow1(H1, V, std::floor((start + end)/2), end, n, Sub, eeR);
		}
		else { 

			// This is an non-leaf case
			// initialize Sub[eeR -1].Eb and Sub[eeR -1].Db
			unsigned int l = Sub[eeR -1].Di.size();
			unsigned int h = Sub[eeR -1].Ei.size();

			//std::vector<long int> p(h, -1);
			//std::vector<long int> z(l, -1);
			//std::vector<unsigned int> t(h, 0);
			Sub[eeR -1].E.assign(h, 0);
			std::iota(Sub[eeR -1].E.begin(), Sub[eeR -1].E.end(), 0);
			Sub[eeR -1].Eb.assign(h, -1); 
			Sub[eeR -1].Db.assign(l, -1);
			Decide_Eb_Db_R1(Sub[eeR -1].Di, Sub[eeR -1].Ei, Sub[eeR -1].Db, Sub[eeR -1].Eb, Sub[eeR -1].E);

			// initialize other attributes of this subproblem
			//std::vector<float> v(l, 0);
			//std::vector<unsigned int> w(l, 0);
			Sub[eeR -1].Dv.assign(l, 0); 
			Sub[eeR -1].Dp.assign(l, 0);
			Sub[eeR -1].D.assign(l, 0);
			std::iota(Sub[eeR -1].D.begin(), Sub[eeR -1].D.end(), 0);

			//std::vector<float> q(h, 0);
			Sub[eeR -1].Ev.assign(h, 0);
			Sub[eeR -1].Ep.assign(h, 0);
			std::pair<long int, long int> dummy_pair = std::make_pair(-1, h+1);
			Sub[eeR -1].S_1.push(dummy_pair); 
			++n;
			//cerr <<"start: " << start+ 1 << ", med: " <<  std::floor((start + 1 + end + 1)/2) << ", n: " <<  n << "\n";
			DivideSubProbByRow1(H1, V, start, std::floor((start + end)/2), n, Sub, eeR);
			++n;
			//cerr <<"med: " << std::floor((start + 1 + end + 1)/2) << ", end: " << end + 1  << ", n: " <<  n << "\n";
			DivideSubProbByRow1(H1, V, std::floor((start + end)/2), end, n, Sub, eeR);
		}

	}
}


#endif