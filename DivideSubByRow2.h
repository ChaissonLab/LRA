#ifndef DIVIDE_SUB_BY_ROW2_H_
#define DIVIDE_SUB_BY_ROW2_H_


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


// This function finds Di and Ei array for non-leaf case
// Note this function didn't count the number of points which have forward diagonal <= the current forward diagonal
// Note: H1[j].inv == 0
void 
ScanPoints_Row2 (std::vector<info> & V, std::vector<Point> & H1, std::vector<long int> & Bi,  unsigned int & s, unsigned int & e, bool & DE, unsigned int & n) {

	std::set<long int> ForwardIndex;
	for (unsigned int i = s; i < e; ++i) {

		unsigned int count = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) {
			if (H1[j].ind == DE and H1[j].inv == 0) { // H1[j].ind == DE == 1 means finding start points
				long int l = static_cast<long int>(H1[j].se.second) + static_cast<long int>(H1[j].se.first); // back diagonal
				ForwardIndex.insert(l);
				++count;				
			}
		}

		if (count != 0) {
			if (DE == 1) V[i].SS_B2.push_back(n);
			else V[i].SS_A2.push_back(n);
		}

	}	

	for (std::set<long int>::reverse_iterator it = ForwardIndex.rbegin(); it != ForwardIndex.rend(); ++it) { // elements in D/E array are in the decreasing order) {
		Bi.push_back(*it);
	}
}



// This function finds Di and Ei array for leaf-case
// Note this function didn't count the number of points which have forward diagonal <= the current forward diagonal
// Note: H1[j].inv == 0
void 
ScanPoints_Row2 (std::vector<info> & V, std::vector<Point> & H1, std::vector<long int> & Bi, std::vector<long int> & Ci,
						 unsigned int & s, unsigned int & e, unsigned int & n) {

	std::set<long int> ForwardIndex1; // ForwardIndex1 is for Ei array
	std::set<long int> ForwardIndex2; // ForwardIndex2 is for Di array
	for (unsigned int i = s; i < e; ++i) {

		unsigned int count1 = 0;
		unsigned int count2 = 0;
		for (unsigned int j = V[i].pstart; j < V[i].pend; ++j) { 
			if (H1[j].ind == 1 and H1[j].inv == 0) { // H1[j].ind  == 1 means finding start points
				long int l = static_cast<long int>(H1[j].se.second) + static_cast<long int>(H1[j].se.first); // back diagonal 
				ForwardIndex1.insert(l);
				++count1;				  
			}
			else if (H1[j].ind == 0 and H1[j].inv == 0) { // H1[j].ind  == 0 means finding end points
				long int r = static_cast<long int>(H1[j].se.second) + static_cast<long int>(H1[j].se.first);
				ForwardIndex2.insert(r);
				++count2;				
			}			
		}

		if (count1 != 0 and count2 != 0) {
			V[i].SS_B2.push_back(n);
			V[i].SS_A2.push_back(n);		
		}

	}	

	for (std::set<long int>::reverse_iterator it = ForwardIndex1.rbegin(); it != ForwardIndex1.rend(); ++it) { // elements in D array are in the decreasing order
		Bi.push_back(*it);
	}
	for (std::set<long int>::reverse_iterator it = ForwardIndex2.rbegin(); it != ForwardIndex2.rend(); ++it) { // elements in E array are in the decreasing order
		Ci.push_back(*it);
	}
}

// This function will decide Eb and Db array
void
Decide_Eb_Db_R2 (std::vector<long int> & Di, std::vector<long int> & Ei, std::vector<long int> & Db, std::vector<long int> & Eb, std::vector<unsigned int> & E) {

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


// This function Divide SubProblems by row for s2 and e2
void
DivideSubProbByRow2 (std::vector<Point> & H1, std::vector<info> & V, unsigned int start, unsigned int end, 
							unsigned int & n, StackOfSubProblems & Sub, int & eeR) { // [start, end) is a half open interval

	if (end == start + 1) { // subproblem A is empty, while B contains only one row. This is a leaf case.

		Subproblem ss = Subproblem(n);
		Sub.Push_Back(eeR, ss); // ss is a subproblem which Di and Ei coming from one row
		++eeR;

		// scan the points to determine Ei and Di
		ScanPoints_Row2(V, H1, Sub[eeR - 1].Ei, Sub[eeR - 1].Di, start, end, n);	

		if (!Sub[eeR - 1].Ei.empty() and !Sub[eeR - 1].Di.empty()) { 
			
			// initialize Sub[eeR - 1]
			unsigned int l = Sub[eeR - 1].Di.size();
			unsigned int h = Sub[eeR - 1].Ei.size();

			Sub[eeR - 1].E.assign(h, 0);
			std::iota(Sub[eeR - 1].E.begin(), Sub[eeR - 1].E.end(), 0);
			Sub[eeR - 1].Eb.assign(h, -1);
			Sub[eeR - 1].Db.assign(l, -1);
			Decide_Eb_Db_R2(Sub[eeR - 1].Di, Sub[eeR - 1].Ei, Sub[eeR - 1].Db, Sub[eeR - 1].Eb, Sub[eeR - 1].E);

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
			//Sub.pop_back(); // delete subproblem ss
			Sub.ClearSingle(eeR);
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
		ScanPoints_Row2(V, H1, Sub[eeR -1].Di, start, med, DE, n);
		// scan the points to determine Ei
		//cerr << "scan points to determine Ei in ["<< med << ", " << end << ")" << endl;
		DE = 1;
		ScanPoints_Row2(V, H1, Sub[eeR -1].Ei, med, end, DE, n);


		if (Sub[eeR -1].Ei.empty() and Sub[eeR -1].Di.empty()) { // Di is empty and Ei is empty  
			//Sub.pop_back();
			Sub.ClearSingle(eeR);
			--eeR;	
			--n;
		}
		else if (Sub[eeR -1].Ei.empty() and !Sub[eeR -1].Di.empty()) { // Di is non-empty and Ei is empty
			//cerr << "Di is non-empty and Ei is empty: " << n << "\n";
			++n;
			//cerr <<"start: " << start+ 1 << ", med: " <<  std::floor((start + 1 + end + 1)/2) << ", n: " <<  n << "\n";
			DivideSubProbByRow2(H1, V, start, std::floor((start + end)/2), n, Sub, eeR);
		}
		else if (!Sub[eeR -1].Ei.empty() and Sub[eeR -1].Di.empty()) { // Di is empty and Ei is non-empty 
			//cerr << "Di is empty and Ei is non-empty: " << n << "\n";
			++n;
			//cerr <<"med: " << std::floor((start + 1 + end + 1)/2) << ", end: " << end + 1  << ", n: " <<  n << "\n";
			DivideSubProbByRow2(H1, V, std::floor((start + end)/2), end, n, Sub, eeR);
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
			Decide_Eb_Db_R2(Sub[eeR -1].Di, Sub[eeR -1].Ei, Sub[eeR -1].Db, Sub[eeR -1].Eb, Sub[eeR -1].E);

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
			DivideSubProbByRow2(H1, V, start, std::floor((start + end)/2), n, Sub, eeR);
			++n;
			//cerr <<"med: " << std::floor((start + 1 + end + 1)/2) << ", end: " << end + 1  << ", n: " <<  n << "\n";
			DivideSubProbByRow2(H1, V, std::floor((start + end)/2), end, n, Sub, eeR);
		}

	}
}

#endif