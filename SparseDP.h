 // This program is implemented Sparse Dynamic Programming algorithm described in David Eppstein paper
// Author: Jingwen Ren

#ifndef SPARSE_DP_
#define SPARSE_DP_

#include <iostream>
#include <string>
#include <utility>
#include <algorithm> // std::lower_bound
#include <numeric> //std::floor
#include <cmath>
#include <set>
#include <iterator>
#include <assert.h>
#include <chrono> // generate random number
#include <random> // generate random number
#include <ctime>
#include <type_traits>
 

#include "SubProblem.h"
#include "Sorting.h"
#include "SubRountine.h"
#include "Fragment_Info.h"
#include "Info.h"
#include "overload.h"
#include "DivideSubByRow1.h"
#include "DivideSubByCol1.h"
#include "DivideSubByRow2.h"
#include "DivideSubByCol2.h"
#include "Point.h"
#include "TupleOps.h"
#include "Options.h"
#include "Clustering.h"


using std::cerr;
using std::cout;
using std::endl;
using std::iota;


typedef std::pair<unsigned int, unsigned int> Pair;
typedef std::pair<long int, long int> LPair;


// Find the first LPair s in [first, last) with s.second > val 
std::vector<LPair>::iterator
UPPERbound (std::vector<LPair>::iterator first, std::vector<LPair>::iterator last, unsigned int val) {
	
	std::vector<LPair>::iterator it;
	unsigned int count, step;
	count = std::distance(first, last);
	while (count > 0) {
		it = first; step = count/2; std::advance(it, step);
		if (val >= it->second) {
			first = ++it;
			count -= step + 1;
		}
		else count = step;
	}
	return first;
}



// TODO(Jingwen): Change this to first get Ev for all the points. Only use "lower_bound" to retrieve the index
void
FindValueInBlock (long int ForwardDiag, std::stack<LPair> S_1, std::vector<long int> & Ei, std::vector<LPair> & Block, unsigned int & i1, unsigned int & i2) {

	if (i1 >= Block.back().second and i1 < S_1.top().second) {
		i2 = S_1.top().first;
	}
	else {
		std::vector<LPair>::iterator it2 = UPPERbound(Block.begin(), Block.end(), i1); // Find the best candidate index for point Ei[i1]
		i2 = it2->first;
	}
}



void 
PassValueToD1 (unsigned int i, std::vector<Fragment_Info> & Value, const std::vector<Point> & H1, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, long int & ForwardDiag) {

	// If all the subproblems B which contain the point H1[i] have been processed, then we can pass the value of the point H1[i] to the D array for every subproblem A
	assert (Value[i].counter_B_R1 == 0 and Value[i].counter_B_C1 == 0);

	//cerr << "pass the value! Process the SubProblem A of the point " << H1[i] << endl;
	//cerr << "Value[" << i << "].SS_A_R1: " << Value[i].SS_A_R1 << "     Value[" << i << "].SS_A_C1: " << Value[i].SS_A_C1 << "\n";

	//cerr << "-------------Solve SS_A_R1-----------------------" << endl;
	unsigned int u;
	for (unsigned int tt = 0; tt < Value[i].SS_A_R1.size(); ++tt) {
		unsigned int u = Value[i].SS_A_R1[Value[i].SS_A_R1.size() - 1 - tt];
		//cerr << "u: " << u << "\n";

		if (SubR1[u].Ei.empty()) {
			//cerr << "SubR1[" << u << "].Ei is empty" << endl;
			continue;			
		}
		else {

			std::vector<unsigned int>::iterator it = Lower_Bound<std::vector<unsigned int>::iterator, long int>(SubR1[u].D.begin(), SubR1[u].D.end(), ForwardDiag, SubR1[u].Di);

			/*
			for (unsigned int q = *it; q < SubR1[u].Di.size(); ++q) {
				--SubR1[u].counter_D[q];
			}
			*/

			//cerr << "SubR1[" << u << "].Di: " << SubR1[u].Di << endl;
			//cerr << "The index of this point in D array: " << *it << "\n";		
			//cerr << "SubR1[u].Dv[" << *it << "]: " << SubR1[u].Dv[*it] << "   SubR1[u].Dp[" << *it << "]: " << SubR1[u].Dp[*it] << "\n";
			//cerr << "Update the D array Now\n";
			//cerr << "SubR1[u].Dv[" << *it << "]: " << SubR1[u].Dv[*it] << "  Value[" << i << "].val: " << Value[i].val << endl;
			if (SubR1[u].Dv[*it] < Value[i].val) { // TODO(Jingwen) < or <= ????
				SubR1[u].Dv[*it] = Value[i].val;
				SubR1[u].Dp[*it] = i;
				//cerr << "Updated! SubR1[" << u << "].Dv[*it]: " << SubR1[u].Dv[*it] << "   SubR1[" << u << "].Dp[*it]: " << SubR1[u].Dp[*it] << endl;
			}
			//else{cerr <<"Not Updated!" << endl;}
		}
		//cerr << endl;
	}


	//cerr << "----------Solve SS_A_C1-----------" << endl;
	for (unsigned int tt = 0; tt < Value[i].SS_A_C1.size(); ++tt) {
		unsigned int u = Value[i].SS_A_C1[Value[i].SS_A_C1.size() - 1 - tt];
		//cerr << "u: " << u << "\n";

		if (SubC1[u].Ei.empty()) {
			//cerr << "SubC1[" << u << "].Ei is empty" << endl;
			continue;
		}
		else {

			std::vector<unsigned int>::reverse_iterator it = Lower_Bound<std::vector<unsigned int>::reverse_iterator, long int>(SubC1[u].D.rbegin(), SubC1[u].D.rend(), ForwardDiag, SubC1[u].Di);
		
			/*
			for (unsigned int q = *it; q < SubC1[u].Di.size(); ++q) {
				--SubC1[u].counter_D[q];
			}
			*/


			//cerr << "SubC1[" << u << "].Di: " << SubC1[u].Di << endl;
			//cerr << "The index of this point in D array: " << *it << "\n";
			//cerr << "SubC1[u].Dv[" << *it << "]: " << SubC1[u].Dv[*it] << "   SubC1[u].Dp[" << *it << "]: " << SubC1[u].Dp[*it] << "\n";
			//cerr << "Update the D array Now\n";
			//cerr << "SubC1[u].Dv[*it]: " << SubC1[u].Dv[*it] << "  Value[i].val: " << Value[i].val << endl;
			if (SubC1[u].Dv[*it] < Value[i].val) {

				SubC1[u].Dv[*it] = Value[i].val;
				SubC1[u].Dp[*it] = i;

				//cerr << "Updated! SubC1[" << u << "].Dv[*it]: " << SubC1[u].Dv[*it] << "   SubC1[" << u << "].Dp[*it]: " << SubC1[u].Dp[*it] << endl;
			}
			//cerr <<"Not Updated!" << endl;
		}
		//cerr << endl;
	}
	//cerr << endl << endl;
}






void 
PassValueToD2 (unsigned int i, std::vector<Fragment_Info> & Value, const std::vector<Point> & H1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2, long int & BackDiag) {

	// If all the subproblems B which contain the point H1[i] have been processed, then we can pass the value of the point H1[i] to the D array for every subproblem A
	assert (Value[i].counter_B_R2 == 0 and Value[i].counter_B_C2 == 0);

	//cerr << "pass the value! Process the SubProblem A of the point " << H1[i] << endl;
	//cerr << "Value[" << i << "].SS_A_R2: " << Value[i].SS_A_R2 << "     Value[" << i << "].SS_A_C2: " << Value[i].SS_A_C2 << "\n";

	//cerr << "-------------Solve SS_A_R2-----------------------" << endl;
	unsigned int u;
	for (unsigned int tt = 0; tt < Value[i].SS_A_R2.size(); ++tt) {
		unsigned int u = Value[i].SS_A_R2[Value[i].SS_A_R2.size() - 1 - tt];
		//cerr << "u: " << u << "\n";

		if (SubR2[u].Ei.empty()) {
			//cerr << "SubR2[" << u << "].Ei is empty" << endl;
			continue;			
		}
		else {

			std::vector<unsigned int>::reverse_iterator it = Lower_Bound<std::vector<unsigned int>::reverse_iterator, long int>(SubR2[u].D.rbegin(), SubR2[u].D.rend(), BackDiag, SubR2[u].Di);
			/*
			for (unsigned int q = *it; q < SubR2[u].Di.size(); ++q) {
				--SubR2[u].counter_D[q];
			}
			*/

			//cerr << "SubR2[" << u << "].Di: " << SubR2[u].Di << endl;
			//cerr << "The index of this point in D array: " << *it << "\n";		
			//cerr << "SubR2[u].Dv[" << *it << "]: " << SubR2[u].Dv[*it] << "   SubR2[u].Dp[" << *it << "]: " << SubR2[u].Dp[*it] << "\n";
			//cerr << "Update the D array Now\n";
			//cerr << "SubR2[u].Dv[" << *it << "]: " << SubR2[u].Dv[*it] << "  Value[" << i << "].val: " << Value[i].val << endl;
			if (SubR2[u].Dv[*it] < Value[i].val) { // TODO(Jingwen) < or <= ????
				SubR2[u].Dv[*it] = Value[i].val;
				SubR2[u].Dp[*it] = i;
				//cerr << "Updated! SubR2[" << u << "].Dv[*it]: " << SubR2[u].Dv[*it] << "   SubR2[" << u << "].Dp[*it]: " << SubR2[u].Dp[*it] << endl;
			}
			//else{cerr <<"Not Updated!" << endl;}
		}
		//cerr << endl;
	}


	//cerr << "----------Solve SS_A_C2-----------" << endl;
	for (unsigned int tt = 0; tt < Value[i].SS_A_C2.size(); ++tt) {
		unsigned int u = Value[i].SS_A_C2[Value[i].SS_A_C2.size() - 1 - tt];
		//cerr << "u: " << u << "\n";

		if (SubC2[u].Ei.empty()) {
			//cerr << "SubC2[" << u << "].Ei is empty" << endl;
			continue;
		}
		else {

			std::vector<unsigned int>::iterator it = Lower_Bound<std::vector<unsigned int>::iterator, long int>(SubC2[u].D.begin(), SubC2[u].D.end(), BackDiag, SubC2[u].Di);
		
			/*
			for (unsigned int q = *it; q < SubC2[u].Di.size(); ++q) {
				--SubC2[u].counter_D[q];
			}
			*/


			//cerr << "SubC2[" << u << "].Di: " << SubC2[u].Di << endl;
			//cerr << "The index of this point in D array: " << *it << "\n";
			//cerr << "SubC2[u].Dv[" << *it << "]: " << SubC2[u].Dv[*it] << "   SubC2[u].Dp[" << *it << "]: " << SubC2[u].Dp[*it] << "\n";
			//cerr << "Update the D array Now\n";
			//cerr << "SubC2[u].Dv[*it]: " << SubC2[u].Dv[*it] << "  Value[i].val: " << Value[i].val << endl;
			if (SubC2[u].Dv[*it] < Value[i].val) {

				SubC2[u].Dv[*it] = Value[i].val;
				SubC2[u].Dp[*it] = i;

				//cerr << "Updated! SubC2[" << u << "].Dv[*it]: " << SubC2[u].Dv[*it] << "   SubC2[" << u << "].Dp[*it]: " << SubC2[u].Dp[*it] << endl;
			}
			//cerr <<"Not Updated!" << endl;
		}
		//cerr << endl;
	}
	//cerr << endl << endl;
}




/*
//This function is for fragments which are resulting from MergeSplit step. 
// Each fragment has different length, so we need to pass vector<Cluster> vt to the function
void ProcessPoint (const std::vector<Cluster> & FragInput, const std::vector<Point> & H1, const std::vector<unsigned int> & H3, std::vector<info> & V, 
						StackOfSubProblems & SubR, StackOfSubProblems & SubC, std::vector<Fragment_Info> & Value, const std::vector<float> & LookUpTable, Options &opts) {

	for (unsigned int i = 0; i < H3.size(); ++i) { // process points by back diagonal

		long int ForwardDiag = static_cast<long int>(H1[H3[i]].se.second) - static_cast<long int>(H1[H3[i]].se.first);

		//cerr << "\n\n\n\nprocessing point " << H1[H3[i]] << endl;
		//cerr << "the Value[" << H1[H3[i]].frag_num << "] of this point:  " <<  Value[H1[H3[i]].frag_num] << "\n";
		
		unsigned int ii = H1[H3[i]].frag_num;


		if (H1[H3[i]].ind == 1) { // H1[H3[i]] is a start point
			//cerr << "----------------------------------------this is a start point------------------------------------" << endl;
			//cerr << "dealing with the subproblem B first. Solve Value[" << H1[H3[i]].frag_num << "].SS_B_R and SS_B_C" << endl;
			//cerr <<"----------------Solve SS_B_R -------------------" << endl;

			// For each subproblem B_R that point H1[H3[i]] is in, get Ev[ForwardDiag] and update Value[H1[H3[i]].frag_num].val
			for (unsigned int k = 0; k < Value[ii].SS_B_R.size(); ++k) {
				unsigned int j = Value[ii].SS_B_R[Value[ii].SS_B_R.size() - 1 - k];
				//cerr << "Solving SubR[" << j << "]: " << SubR[j]<< "\n";

				// If subproblem SubR[j] is a leaf case, then 
				if (SubR[j].Di.empty()) {
					//cerr << "SubR[" << j << "] is a leaf case or non-leaf case but Di is empty" << "\n";
					--Value[ii].counter_B_R;
					continue;
				}
				else {
					// Then subproblem SubR[j] is non-leaf case. find the index of the point in E array
					std::vector<unsigned int>::iterator t = Lower_Bound<std::vector<unsigned int>::iterator,long int>(SubR[j].E.begin(), SubR[j].E.end(), ForwardDiag, SubR[j].Ei);
					//cerr << "SubR[" << j << "] is a non-leaf case!      Find the index of the point in E array" << "/n";
					//cerr << "SubR[" << j << "].Ei: " << SubR[j].Ei << "/n";
					//cerr << "The index of this point in E array: " << *t << "\n";

					if (SubR[j].Eb[*t] == -1) {
						//cerr << "SubR[" << j << "] is a non-leaf case but there is no forward diag smaller than the current" << endl;
						--Value[ii].counter_B_R;
						continue;
					}
					else {
						//assert(SubR[j].counter_D[SubR[j].Eb[*t]] == 0);
						//cerr << "SubR[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already." << endl;
						//cerr << "Start to compute the maximization structure.\n";

						SubR[j].now = SubR[j].Eb[*t];
						Maximization (SubR[j].now, SubR[j].last, SubR[j].Di, SubR[j].Ei, SubR[j].Dv, SubR[j].Db, SubR[j].Block, SubR[j].S_1, LookUpTable, opts); // TODO(Jingwen) anything change for SubC????
						SubR[j].last = SubR[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubR[j].S_1, SubR[j].Ei, SubR[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";

						 // FragInput is vector<Cluster> resulting from merge-split step
						SubR[j].Ev[i1] = SubR[j].Dv[i2] + w(SubR[j].Di[i2], SubR[j].Ei[i1], LookUpTable, opts) + 
											std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1); 
						SubR[j].Ep[i1] = i2;						

						//cerr << "SubR[" << j << "].Ev[" << i1 << "]: " << SubR[j].Ev[i1] << ", SubR[" << j << "].Ep[" << i1 << "]: " << SubR[j].Ep[i1] << "\n"; 


						// Update the value of this point
						if (Value[ii].val < SubR[j].Ev[i1]) {  
							Value[ii].val = SubR[j].Ev[i1];
							Value[ii].prev_sub = SubR[j].num;
							Value[ii].prev_ind = i1;
							Value[ii].prev = 1; // the best value comes from SubR
							//cerr << "update the value of this point\n";
							//cerr << "Value[ii].val: " << Value[ii].val << ", Value[ii].prev_sub: " << Value[ii].prev_sub << ", Value[ii].prev_ind: " << Value[ii].prev_ind << 
							//", Value[ii].prev: " << Value[ii].prev << endl;
						}			
						--Value[ii].counter_B_R;	
						//cerr << "Do not update the value of this point\n";
					}
				}
			}		


			//cerr << endl << endl;
			//cerr << "--------------Solve SS_B_C ----------------" << endl;
			// For each subproblem B_C that point H1[i] is in, get Ev[H1[i].first - H1[i].second] and update Value[H1[i]].val
			for (unsigned int k = 0; k < Value[ii].SS_B_C.size(); ++k) {
				unsigned int j = Value[ii].SS_B_C[Value[ii].SS_B_C.size() - 1 - k];
				//cerr << "SubC[" << j << "]: " << SubC[j]<< "\n";

				// If subproblem SubC[j] is a leaf case, then 
				if (SubC[j].Di.empty()) {

					//cerr << "SubC[" << j << "] is a leaf case or it's a non-leaf case with an empty Di" << "\n";
					--Value[ii].counter_B_C;
					continue;
				}
				else {

					// find the index of this point in E array
					std::vector<unsigned int>::reverse_iterator t = Lower_Bound<std::vector<unsigned int>::reverse_iterator,long int>(SubC[j].E.rbegin(), SubC[j].E.rend(), ForwardDiag, SubC[j].Ei);
					//cerr << "SubC[" << j << "] is a non-leaf case!      Find the index of the point in E array" << "/n";
					//cerr << "SubC[" << j << "].Ei: " << SubC[j].Ei << "/n";
					//cerr << "The index of this point in E array: " << *t << "\n";

					if (SubC[j].Eb[*t] == -1) {
						//cerr << "SubC[" << j << "] is a non-leaf case but there is no forward diag larger than the current\n ";
						--Value[ii].counter_B_C;
						continue;
					}
					else {

						//assert(SubC[j].counter_D[SubC[j].Eb[*t]] == 0);
						//cerr << "SubC[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already.\n";
						//cerr << "Start to compute the maximization structure.\n";
						SubC[j].now = SubC[j].Eb[*t]; 
						Maximization (SubC[j].now, SubC[j].last, SubC[j].Di, SubC[j].Ei, SubC[j].Dv, SubC[j].Db, SubC[j].Block, SubC[j].S_1, LookUpTable, opts); 
						SubC[j].last = SubC[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores 
						unsigned int i2;
						FindValueInBlock(ForwardDiag, SubC[j].S_1, SubC[j].Ei, SubC[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					
						// FragInput is vector<Cluster> resulting from merge-split step
						SubC[j].Ev[i1] = SubC[j].Dv[i2] + w(SubC[j].Di[i2], SubC[j].Ei[i1], LookUpTable, opts) + 
											std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1); 
						SubC[j].Ep[i1] = i2;							

						//cerr << "SubC[" << j << "].Ev[" << i1 << "]: " << SubC[j].Ev[i1] << ", SubC[" << j << "].Ep[" << i1 << "]: " << SubC[j].Ep[i1] << "\n"; 

						// Update the value of this point
						if (Value[ii].val < SubC[j].Ev[i1]) {  
							Value[ii].val = SubC[j].Ev[i1];
							Value[ii].prev_sub = SubC[j].num;
							Value[ii].prev_ind = i1;
							Value[ii].prev = 0; // the best value comes from SubC
							//cerr << "update the value of this point\n";
							//cerr << "Value[ii].val: " << Value[ii].val << ", Value[ii].prev_sub: " << Value[ii].prev_sub << ", Value[ii].prev_ind: " << Value[ii].prev_ind << 
							//", Value[ii].prev: " << Value[ii].prev << endl;
						}			
						--Value[ii].counter_B_C;	
						//cerr << "Do not update the value of this point\n";
					}
				}
			}

		}

		else { // H1[H3[i]] is an end point
			//cerr << "---------------------------------This is an end point--------------------------------"<< endl;
			//cerr << "Pass the value Value[ii].val to SS_A_R and SS_A_C" << endl;

			// If all the subproblems B that the point H1[i] is in are already processed, then Value[H1[i]].val is ready
			// For each subproblem A that the point H1[i] is in, Pass the Value[H1[i]].val to the D array 
			PassValueToD(ii, Value, H1, SubR, SubC, ForwardDiag);
		}
	}
}
*/


// This function is for fragments which are not resulting from MergeSplit step. 
// Note: Each fragment has the same length
void 
ProcessPoint (const std::vector<Point> & H1, std::vector<info> & V, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1,
				 StackOfSubProblems & SubR2, StackOfSubProblems & SubC2, std::vector<Fragment_Info> & Value, Options & opts, const std::vector<float> & LookUpTable, int rate) {

//ProcessPoint (const std::vector<Point> & H1, const std::vector<unsigned int> & H3, std::vector<info> & V, StackOfSubProblems & SubR, StackOfSubProblems & SubC,
//				  std::vector<Fragment_Info> & Value, Options & opts, const std::vector<float> & LookUpTable, int rate) {


	for (unsigned int i = 0; i < H1.size(); ++i) { // process points by row

		// TODO(Jingwen): change this
		long int ForwardDiag = static_cast<long int>(H1[i].se.second) - static_cast<long int>(H1[i].se.first);
		long int BackDiag = static_cast<long int>(H1[i].se.second) + static_cast<long int>(H1[i].se.first);

		//cerr << "\n\n\n\nprocessing point " << H1[H3[i]] << endl;
		//cerr << "the Value[" << H1[H3[i]].frag_num << "] of this point:  " <<  Value[H1[H3[i]].frag_num] << "\n";
		
		unsigned int ii = H1[i].frag_num;


		if (H1[i].ind == 1 and H1[i].inv == 1) { // H1[i] is a start point s1
			//cerr << "----------------------------------------this is a start point (s1) ------------------------------------" << endl;
			//cerr << "dealing with the subproblem B first. Solve Value[" << H1[i].frag_num << "].SS_B_R1 and SS_B_C1" << endl;
			//cerr <<"----------------Solve SS_B_R1 -------------------" << endl;

			// For each subproblem B_R1 that point H1[i] is in, get Ev[ForwardDiag] and update Value[H1[i]].frag_num].val
			for (unsigned int k = 0; k < Value[ii].SS_B_R1.size(); ++k) {
				unsigned int j = Value[ii].SS_B_R1[Value[ii].SS_B_R1.size() - 1 - k];
				//cerr << "Solving SubR1[" << j << "]: " << SubR1[j]<< "\n";

				// If subproblem SubR1[j] is a leaf case, then 
				if (SubR1[j].Di.empty()) {
					//cerr << "SubR1[" << j << "] is a leaf case or non-leaf case but Di is empty" << "\n";
					--Value[ii].counter_B_R1;
					continue;
				}
				else {
					// Then subproblem SubR1[j] is non-leaf case. find the index of the point in E array
					std::vector<unsigned int>::iterator t = Lower_Bound<std::vector<unsigned int>::iterator,long int>(SubR1[j].E.begin(), SubR1[j].E.end(), ForwardDiag, SubR1[j].Ei);
					//cerr << "SubR1[" << j << "] is a non-leaf case!      Find the index of the point in E array" << "/n";
					//cerr << "SubR1[" << j << "].Ei: " << SubR1[j].Ei << "/n";
					//cerr << "The index of this point in E array: " << *t << "\n";

					if (SubR1[j].Eb[*t] == -1) {
						//cerr << "SubR1[" << j << "] is a non-leaf case but there is no forward diag smaller than the current" << endl;
						--Value[ii].counter_B_R1;
						continue;
					}
					else {
						//assert(SubR1[j].counter_D[SubR1[j].Eb[*t]] == 0);
						//cerr << "SubR1[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already." << endl;
						//cerr << "Start to compute the maximization structure.\n";

						SubR1[j].now = SubR1[j].Eb[*t];
						Maximization (SubR1[j].now, SubR1[j].last, SubR1[j].Di, SubR1[j].Ei, SubR1[j].Dv, SubR1[j].Db, SubR1[j].Block, SubR1[j].S_1, LookUpTable, opts); // TODO(Jingwen) anything change for SubC????
						SubR1[j].last = SubR1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubR1[j].S_1, SubR1[j].Ei, SubR1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";


//						SubR1[j].Ev[i1] = SubR1[j].Dv[i2] + w(SubR1[j].Di[i2], SubR1[j].Ei[i1], LookUpTable, opts) + opts.globalK + 1; 
//						SubR1[j].Ep[i1] = i2;							

						SubR1[j].Ev[i1] = SubR1[j].Dv[i2] + w(SubR1[j].Di[i2], SubR1[j].Ei[i1], LookUpTable, opts) + (opts.globalK + 1) * rate; 
						SubR1[j].Ep[i1] = i2;							

						//cerr << "SubR1[" << j << "].Ev[" << i1 << "]: " << SubR1[j].Ev[i1] << ", SubR1[" << j << "].Ep[" << i1 << "]: " << SubR1[j].Ep[i1] << "\n"; 


						// Update the value of this point
						if (Value[ii].val < SubR1[j].Ev[i1]) {  
							Value[ii].val = SubR1[j].Ev[i1];
							Value[ii].prev_sub = SubR1[j].num;
							Value[ii].prev_ind = i1;
							Value[ii].prev = 1; // the best value comes from row subproblem
							Value[ii].inv = 1; // the best value comes from SubR1
							//cerr << "update the value of this point\n";
							//cerr << "Value[ii].val: " << Value[ii].val << ", Value[ii].prev_sub: " << Value[ii].prev_sub << ", Value[ii].prev_ind: " << Value[ii].prev_ind << 
							//", Value[ii].prev: " << Value[ii].prev << endl;
						}			
						--Value[ii].counter_B_R1;	
						//cerr << "Do not update the value of this point\n";
					}
				}
			}		


			//cerr << endl << endl;
			//cerr << "--------------Solve SS_B_C1 ----------------" << endl;
			// For each subproblem B_C1 that point H1[i] is in, get Ev[H1[i].first - H1[i].second] and update Value[H1[i]].val
			for (unsigned int k = 0; k < Value[ii].SS_B_C1.size(); ++k) {
				unsigned int j = Value[ii].SS_B_C1[Value[ii].SS_B_C1.size() - 1 - k];
				//cerr << "SubC1[" << j << "]: " << SubC1[j]<< "\n";

				// If subproblem SubC1[j] is a leaf case, then 
				if (SubC1[j].Di.empty()) {

					//cerr << "SubC1[" << j << "] is a leaf case or it's a non-leaf case with an empty Di" << "\n";
					--Value[ii].counter_B_C1;
					continue;
				}
				else {

					// find the index of this point in E array
					std::vector<unsigned int>::reverse_iterator t = Lower_Bound<std::vector<unsigned int>::reverse_iterator,long int>(SubC1[j].E.rbegin(), SubC1[j].E.rend(), ForwardDiag, SubC1[j].Ei);
					//cerr << "SubC1[" << j << "] is a non-leaf case!      Find the index of the point in E array" << "/n";
					//cerr << "SubC1[" << j << "].Ei: " << SubC1[j].Ei << "/n";
					//cerr << "The index of this point in E array: " << *t << "\n";

					if (SubC1[j].Eb[*t] == -1) {
						//cerr << "SubC1[" << j << "] is a non-leaf case but there is no forward diag larger than the current\n ";
						--Value[ii].counter_B_C1;
						continue;
					}
					else {

						//assert(SubC1[j].counter_D[SubC1[j].Eb[*t]] == 0);
						//cerr << "SubC1[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already.\n";
						//cerr << "Start to compute the maximization structure.\n";
						SubC1[j].now = SubC1[j].Eb[*t]; 
						Maximization (SubC1[j].now, SubC1[j].last, SubC1[j].Di, SubC1[j].Ei, SubC1[j].Dv, SubC1[j].Db, SubC1[j].Block, SubC1[j].S_1, LookUpTable, opts); 
						SubC1[j].last = SubC1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores 
						unsigned int i2;
						FindValueInBlock(ForwardDiag, SubC1[j].S_1, SubC1[j].Ei, SubC1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					
	//					SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts) + opts.globalK + 1; 
	//					SubC1[j].Ep[i1] = i2;							

						SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts) + (opts.globalK + 1) * rate; 
						SubC1[j].Ep[i1] = i2;							

						//cerr << "SubC1[" << j << "].Ev[" << i1 << "]: " << SubC1[j].Ev[i1] << ", SubC1[" << j << "].Ep[" << i1 << "]: " << SubC1[j].Ep[i1] << "\n"; 

						// Update the value of this point
						if (Value[ii].val < SubC1[j].Ev[i1]) {  
							Value[ii].val = SubC1[j].Ev[i1];
							Value[ii].prev_sub = SubC1[j].num;
							Value[ii].prev_ind = i1;
							Value[ii].prev = 0; // the best value comes from col subproblem
							Value[ii].inv = 1; // the best value comes from SubC1
							//cerr << "update the value of this point\n";
							//cerr << "Value[ii].val: " << Value[ii].val << ", Value[ii].prev_sub: " << Value[ii].prev_sub << ", Value[ii].prev_ind: " << Value[ii].prev_ind << 
							//", Value[ii].prev: " << Value[ii].prev << endl;
						}			
						--Value[ii].counter_B_C1;	
						//cerr << "Do not update the value of this point\n";
					}
				}
			}
		}
		else if (H1[i].ind == 1 and H1[i].inv == 0) { // H1[i] is a start point s2

			//cerr << "----------------------------------------this is a start point (s2) ------------------------------------" << endl;
			//cerr << "dealing with the subproblem B first. Solve Value[" << H1[i].frag_num << "].SS_B_R2 and SS_B_C2" << endl;
			//cerr <<"----------------Solve SS_B_R2 -------------------" << endl;

			// For each subproblem B_R2 that point H1[i] is in, get Ev[BackDiag] and update Value[H1[i]].frag_num].val
			for (unsigned int k = 0; k < Value[ii].SS_B_R2.size(); ++k) {
				unsigned int j = Value[ii].SS_B_R2[Value[ii].SS_B_R2.size() - 1 - k];
				//cerr << "Solving SubR2[" << j << "]: " << SubR2[j]<< "\n";

				// If subproblem SubR2[j] is a leaf case, then 
				if (SubR2[j].Di.empty()) {
					//cerr << "SubR2[" << j << "] is a leaf case or non-leaf case but Di is empty" << "\n";
					--Value[ii].counter_B_R2;
					continue;
				}
				else {
					// Then subproblem SubR2[j] is non-leaf case. find the index of the point in E array
					std::vector<unsigned int>::reverse_iterator t = Lower_Bound<std::vector<unsigned int>::reverse_iterator,long int>(SubR2[j].E.rbegin(), SubR2[j].E.rend(), BackDiag, SubR2[j].Ei);
					//cerr << "SubR2[" << j << "] is a non-leaf case!      Find the index of the point in E array" << "/n";
					//cerr << "SubR2[" << j << "].Ei: " << SubR2[j].Ei << "/n";
					//cerr << "The index of this point in E array: " << *t << "\n";

					if (SubR2[j].Eb[*t] == -1) {
						//cerr << "SubR2[" << j << "] is a non-leaf case but there is no back diag larger than the current" << endl;
						--Value[ii].counter_B_R2;
						continue;
					}
					else {
						//assert(SubR2[j].counter_D[SubR2[j].Eb[*t]] == 0);
						//cerr << "SubR2[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already." << endl;
						//cerr << "Start to compute the maximization structure.\n";

						SubR2[j].now = SubR2[j].Eb[*t];
						Maximization (SubR2[j].now, SubR2[j].last, SubR2[j].Di, SubR2[j].Ei, SubR2[j].Dv, SubR2[j].Db, SubR2[j].Block, SubR2[j].S_1, LookUpTable, opts); // TODO(Jingwen) anything change for SubC????
						SubR2[j].last = SubR2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that BackDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for BackDiag
						FindValueInBlock(BackDiag, SubR2[j].S_1, SubR2[j].Ei, SubR2[j].Block, i1, i2);
						//cerr << "the index in Ei that BackDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for BackDiag ---- i2: " << i2 << "\n";

						SubR2[j].Ev[i1] = SubR2[j].Dv[i2] + w(SubR2[j].Di[i2], SubR2[j].Ei[i1], LookUpTable, opts) + opts.globalK + 1; 
						SubR2[j].Ep[i1] = i2;							

						//cerr << "SubR2[" << j << "].Ev[" << i1 << "]: " << SubR2[j].Ev[i1] << ", SubR2[" << j << "].Ep[" << i1 << "]: " << SubR2[j].Ep[i1] << "\n"; 


						// Update the value of this point
						if (Value[ii].val < SubR2[j].Ev[i1]) {  
							Value[ii].val = SubR2[j].Ev[i1];
							Value[ii].prev_sub = SubR2[j].num;
							Value[ii].prev_ind = i1;
							Value[ii].prev = 1; // the best value comes from row subproblem
							Value[ii].inv = 0; // the best value comes from SubR2
							//cerr << "update the value of this point\n";
							//cerr << "Value[ii].val: " << Value[ii].val << ", Value[ii].prev_sub: " << Value[ii].prev_sub << ", Value[ii].prev_ind: " << Value[ii].prev_ind << 
							//", Value[ii].prev: " << Value[ii].prev << endl;
						}			
						--Value[ii].counter_B_R2;	
						//cerr << "Do not update the value of this point\n";
					}
				}
			}		


			//cerr << endl << endl;
			//cerr << "--------------Solve SS_B_C2 ----------------" << endl;
			// For each subproblem B_C2 that point H1[i] is in, get Ev[H1[i].first - H1[i].second] and update Value[H1[i]].val
			for (unsigned int k = 0; k < Value[ii].SS_B_C2.size(); ++k) {
				//cerr << "k: " << k << endl;
				unsigned int j = Value[ii].SS_B_C2[Value[ii].SS_B_C2.size() - 1 - k];
				//cerr << "j: " << j << endl;
				//cerr << "SubC2[" << j << "]: " << SubC2[j]<< "\n";

				// If subproblem SubC2[j] is a leaf case, then 
				if (SubC2[j].Di.empty()) {

					//cerr << "SubC2[" << j << "] is a leaf case or it's a non-leaf case with an empty Di" << "\n";
					--Value[ii].counter_B_C2;
					continue;
				}
				else {

					// find the index of this point in E array
					std::vector<unsigned int>::iterator t = Lower_Bound<std::vector<unsigned int>::iterator,long int>(SubC2[j].E.begin(), SubC2[j].E.end(), BackDiag, SubC2[j].Ei);
					//cerr << "SubC2[" << j << "] is a non-leaf case!      Find the index of the point in E array" << "/n";
					//cerr << "SubC2[" << j << "].Ei: " << SubC2[j].Ei << "/n";
					//cerr << "The index of this point in E array: " << *t << "\n";

					if (SubC2[j].Eb[*t] == -1) {
						//cerr << "SubC2[" << j << "] is a non-leaf case but there is no forward diag larger than the current\n ";
						--Value[ii].counter_B_C2;
						continue;
					}
					else {

						//assert(SubC2[j].counter_D[SubC2[j].Eb[*t]] == 0);
						//cerr << "SubC2[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already.\n";
						//cerr << "Start to compute the maximization structure.\n";
						SubC2[j].now = SubC2[j].Eb[*t]; 
						Maximization (SubC2[j].now, SubC2[j].last, SubC2[j].Di, SubC2[j].Ei, SubC2[j].Dv, SubC2[j].Db, SubC2[j].Block, SubC2[j].S_1, LookUpTable, opts); 
						SubC2[j].last = SubC2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores 
						unsigned int i2;
						FindValueInBlock(BackDiag, SubC2[j].S_1, SubC2[j].Ei, SubC2[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					

						SubC2[j].Ev[i1] = SubC2[j].Dv[i2] + w(SubC2[j].Di[i2], SubC2[j].Ei[i1], LookUpTable, opts) + opts.globalK + 1; 
						SubC2[j].Ep[i1] = i2;							

						//cerr << "SubC2[" << j << "].Ev[" << i1 << "]: " << SubC2[j].Ev[i1] << ", SubC2[" << j << "].Ep[" << i1 << "]: " << SubC2[j].Ep[i1] << "\n"; 

						// Update the value of this point
						if (Value[ii].val < SubC2[j].Ev[i1]) {  
							Value[ii].val = SubC2[j].Ev[i1];
							Value[ii].prev_sub = SubC2[j].num;
							Value[ii].prev_ind = i1;
							Value[ii].prev = 0; // the best value comes from col subproblem
							Value[ii].inv = 0; // the best value comes from SubC2
							//cerr << "update the value of this point\n";
							//cerr << "Value[ii].val: " << Value[ii].val << ", Value[ii].prev_sub: " << Value[ii].prev_sub << ", Value[ii].prev_ind: " << Value[ii].prev_ind << 
							//", Value[ii].prev: " << Value[ii].prev << endl;
						}			
						--Value[ii].counter_B_C2;	
						//cerr << "Do not update the value of this point\n";
					}
				}
			}
		}
		else if (H1[i].ind == 0 and H1[i].inv == 1) { // H1[i] is an end point (e1)
			//cerr << "---------------------------------This is an end point (e1) --------------------------------"<< endl;
			//cerr << "Pass the value Value[ii].val to SS_A_R1 and SS_A_C1" << endl;

			// If all the subproblems B that the point H1[i] is in are already processed, then Value[H1[i]].val is ready
			// For each subproblem A that the point H1[i] is in, Pass the Value[H1[i]].val to the D array 
			PassValueToD1(ii, Value, H1, SubR1, SubC1, ForwardDiag);
		}
		else { // H1[i] is an end point (e2)
			//cerr << "---------------------------------This is an end point (e2) --------------------------------"<< endl;
			//cerr << "Pass the value Value[ii].val to SS_A_R2 and SS_A_C2" << endl;

			// If all the subproblems B that the point H1[i] is in are already processed, then Value[H1[i]].val is ready
			// For each subproblem A that the point H1[i] is in, Pass the Value[H1[i]].val to the D array 
			PassValueToD2(ii, Value, H1, SubR2, SubC2, BackDiag);
		}
	}
}


void 
TraceBack (StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
					const std::vector<Fragment_Info> & Value, unsigned int i, std::vector<unsigned int> & FinalChain) {

	long int prev_sub = Value[i].prev_sub;
	long int prev_ind = Value[i].prev_ind;
	FinalChain.push_back(i);

	if (prev_sub != -1 and prev_ind != -1) {
		
		if (Value[i].prev == 1 and Value[i].inv == 1) { // The previous subproblem is SubR1
			unsigned int ind = SubR1[prev_sub].Ep[prev_ind];
			TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR1[prev_sub].Dp[ind], FinalChain);
		}
		else if (Value[i].prev == 1 and Value[i].inv == 0) { // The previous subproblem is SubR2
			unsigned int ind = SubR2[prev_sub].Ep[prev_ind];
			TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR2[prev_sub].Dp[ind], FinalChain);
		}
		else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
			unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
			TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC1[prev_sub].Dp[ind], FinalChain);
		}		
		else { // The previous subproblem is SubC2
			unsigned int ind = SubC2[prev_sub].Ep[prev_ind];
			TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC2[prev_sub].Dp[ind], FinalChain);
		}

	}
}

/*
// The input for this function is vector<Cluster> vt which is from Merge_Split step
int SparseDP (const std::vector<Cluster> & FragInput, std::vector<unsigned int> & chain, Options & opts, const std::vector<float> & LookUpTable) {

	std::vector<Point>  H1;
	// get points from FragInput and store them in H1
	for (unsigned int i = 0; i < FragInput.size(); ++i) {
		// insert start point into H1
		Point s;
		H1.push_back(s);
		H1.back().ind = 1; 
		H1.back().frag_num = i;
		H1.back().se.first = FragInput[i].qStart;
		H1.back().se.second = FragInput[i].tStart;

		// insert end point into H1
		Point e;
		H1.push_back(e);
		H1.back().ind = 0;
		H1.back().frag_num = i;
		H1.back().se.first = FragInput[i].qEnd;
		H1.back().se.second = FragInput[i].tEnd;

	}			

	clock_t begin = std::clock();

	//Sort the point by row
	sort(H1.begin(), H1.end(), SortByRowOp<Point>());

	//cerr << "H1: " << H1 << endl;
	std::vector<unsigned int> H2(H1.size());
	std::vector<unsigned int> H3(H1.size());
	iota(H2.begin(), H2.end(), 0);
	iota(H3.begin(), H3.end(), 0);

	//Sort the point by column and by back diagonal
	sort(H2.begin(), H2.end(), SortByColOp<Point, unsigned int>(H1));
	sort(H3.begin(), H3.end(), SortByBackDiagOp<Point, unsigned int>(H1));

	
	// print out H2, and H3

	std::vector<info> Row;
	std::vector<info> Col;
	GetRowInfo(H1, Row);
	GetColInfo(H1, H2, Col);
	//cerr << "Row: " << Row << "\n";
	//cerr << "Col: " << Col << "\n";

	unsigned int n = 0;
	unsigned int m = 0;
	//std::vector<Subproblem> SubR;
	//std::vector<Subproblem> SubC;
	StackOfSubProblems SubR;
	StackOfSubProblems SubC;
	int eeR=0, eeC=0;

	//cerr << "DivideSubByRow\n";
	DivideSubProbByRow1(H1, Row, 0, Row.size(), n, SubR, eeR);
	//cerr << "SubR: " << SubR << endl;


	//cerr << "DivideSubByCol\n";
	DivideSubProbByCol1(H1, H2, Col, 0, Col.size(), m, SubC, eeC);
	//cerr << "SubC: " << SubC << endl;

	std::vector<Fragment_Info> Value(FragInput.size());
	for (unsigned int t = 0; t < Row.size(); ++t) {
		for (unsigned int tt = Row[t].pstart; tt < Row[t].pend; ++tt) {
			unsigned int ii = H1[tt].frag_num;
			if (H1[tt].ind == 1) { //H1[tt] is a start point
				Value[ii].SS_B_R = Row[t].SS_B;
				Value[ii].counter_B_R = Row[t].SS_B.size();
				Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1);
			}
			else { // H1[tt] is an end point
				Value[ii].SS_A_R = Row[t].SS_A;
				Value[ii].counter_A_R = Row[t].SS_A.size();
				Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1);
			}

		}
	}



	for (unsigned int t = 0; t < Col.size(); ++t) {
		for (unsigned int tt = Col[t].pstart; tt < Col[t].pend; ++tt) { 

			unsigned int ii = H1[H2[tt]].frag_num;
			if (H1[H2[tt]].ind == 1) { //H1[H2[tt]] a start point
				Value[ii].SS_B_C = Col[t].SS_B;
				Value[ii].counter_B_C = Col[t].SS_B.size();
				Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1);
			}
			else { // H1[H2[tt]] is an end point
				Value[ii].SS_A_C = Col[t].SS_A;
				Value[ii].counter_A_C = Col[t].SS_A.size();
				Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1);
			}

		}
	}


	//cerr << "Value: " << Value << endl;
	

	//cerr << "ProcessPoint\n";

	ProcessPoint(FragInput, H1, H3, Row, SubR, SubC, Value, LookUpTable, opts);		

	//cerr << "end\n";

	// find the max_value for the FinalChain 
	unsigned int  l = 0;
	float max_value = 0;
	unsigned int max_pos = 0;
	while (l < Value.size()) {
		if (Value[l].val > max_value) {
			max_value = Value[l].val;
			max_pos = l;
		}
		++l;
	}
	
	//cerr << "TraceBack\n";
	// Trace back to get the FinalChain
	// store the index of points
	chain.clear();
	TraceBack(SubR, SubC, Value, max_pos, chain);
	std::reverse(chain.begin(), chain.end());


	// Clear SubR and SubC
	SubR.Clear(eeR);
	SubC.Clear(eeC);

	// get the time for the program
	clock_t end = std::clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//cerr << "Time: " << elapsed_secs << endl;

		// print out the final chain

	return 0;
}
*/


// The input for this function is GenomePairs which is not from Merge_Split step
// Each fragment has the same length
int SparseDP (const GenomePairs & FragInput, std::vector<unsigned int> & chain, Options & opts, const std::vector<float> & LookUpTable, int rate = 1) {
	std::vector<Point> H1;

	// FragInput is vector<GenomePair>
	// get points from FragInput and store them in H1		
	for (unsigned int i = 0; i < FragInput.size(); ++i) {
	
		// insert start point s1 into H1
		Point s1;
		H1.push_back(s1);
		H1.back().ind = 1; // start
		H1.back().inv = 1; // forward direction
		H1.back().frag_num = i;
		H1.back().se.first = FragInput[i].first.pos; 
		H1.back().se.second = FragInput[i].second.pos;	

		// insert start point s2 into H1
		Point s2;
		H1.push_back(s2);
		H1.back().ind = 1; // start
		H1.back().inv = 0; // backward direction
		H1.back().frag_num = i;
		H1.back().se.first = FragInput[i].first.pos; 
		H1.back().se.second = FragInput[i].second.pos + opts.globalK;	

		// insert end point e1 into H1
		Point e1;
		H1.push_back(e1);
		H1.back().ind = 0; // end
		H1.back().inv = 1; // forward direction		
		H1.back().frag_num = i;
		H1.back().se.first = FragInput[i].first.pos + opts.globalK;
		H1.back().se.second = FragInput[i].second.pos + opts.globalK;	

		// insert end point e2 into H1
		Point e2;
		H1.push_back(e2);
		H1.back().ind = 0; // end
		H1.back().inv = 0; // backward direction		
		H1.back().frag_num = i;
		H1.back().se.first = FragInput[i].first.pos + opts.globalK;
		H1.back().se.second = FragInput[i].second.pos;	
	}

	clock_t begin = std::clock();

	//Sort the point by row
	sort(H1.begin(), H1.end(), SortByRowOp<Point>()); // with same q and t coordinates, end point < start point

	//cerr << "H1: " << H1 << endl;
	std::vector<unsigned int> H2(H1.size());
	//std::vector<unsigned int> H3(H1.size()); // TODO(Jingwen): Probably don't need this
	iota(H2.begin(), H2.end(), 0);
	//iota(H3.begin(), H3.end(), 0);

	//Sort the point by column and by back diagonal
	sort(H2.begin(), H2.end(), SortByColOp<Point, unsigned int>(H1));
	//sort(H3.begin(), H3.end(), SortByBackDiagOp<Point, unsigned int>(H1));

	
	// print out H2, and H3
	/*
	cerr << "H2: [";
	for (unsigned int t = 0; t < H2.size(); ++t) {
		cerr << H1[H2[t]] << endl;
	}	
	cerr << "]\n";


	cerr << "H3: [";
	for (unsigned int t = 0; t < H2.size(); ++t) {
		cerr << H1[H3[t]] << endl;
	}	
	cerr << "]\n";
	*/

	std::vector<info> Row;
	std::vector<info> Col;
	GetRowInfo(H1, Row);
	GetColInfo(H1, H2, Col);
	//cerr << "Row: " << Row << "\n";
	//cerr << "Col: " << Col << "\n";

	unsigned int n1 = 0;
	unsigned int m1 = 0;
	unsigned int n2 = 0;
	unsigned int m2 = 0;
	//std::vector<Subproblem> SubR;
	//std::vector<Subproblem> SubC;
	StackOfSubProblems SubR1;
	StackOfSubProblems SubC1;
	int eeR1 = 0, eeC1 = 0;

	StackOfSubProblems SubR2;
	StackOfSubProblems SubC2;
	int eeR2 = 0, eeC2 = 0;

	//cerr << "DivideSubByRow\n";
	DivideSubProbByRow1(H1, Row, 0, Row.size(), n1, SubR1, eeR1);
	//cerr << "SubR: " << SubR << endl;

	//cerr << "DivideSubByCol\n";
	DivideSubProbByCol1(H1, H2, Col, 0, Col.size(), m1, SubC1, eeC1);
	//cerr << "SubC: " << SubC << endl;

	DivideSubProbByRow2(H1, Row, 0, Row.size(), n2, SubR2, eeR2);	
	DivideSubProbByCol2(H1, H2, Col, 0, Col.size(), m2, SubC2, eeC2);


	// Get SS_A_R1, SS_B_R1, SS_A_R2 and SS_B_R2 for each fragment
	std::vector<Fragment_Info> Value(FragInput.size());
	for (unsigned int t = 0; t < Row.size(); ++t) {
		for (unsigned int tt = Row[t].pstart; tt < Row[t].pend; ++tt) {

			unsigned int ii = H1[tt].frag_num;

			if (H1[tt].ind == 1 and H1[tt].inv == 1) { //H1[tt] is a start point (s1)
				Value[ii].SS_B_R1 = Row[t].SS_B1;
				Value[ii].counter_B_R1 = Row[t].SS_B1.size();
				Value[ii].val = (opts.globalK + 1);
			}
			else if (H1[tt].ind == 0 and H1[tt].inv == 1) { // H1[tt] is an end point (e1)
				Value[ii].SS_A_R1 = Row[t].SS_A1;
				Value[ii].counter_A_R1 = Row[t].SS_A1.size();
				//Value[ii].val = (opts.globalK + 1); // TODO(Jingwen): make sure that this is redundant
			}
			else if (H1[tt].ind == 1 and H1[tt].inv == 0) { //H1[tt] is a start point (s2)
				Value[ii].SS_B_R2 = Row[t].SS_B2;
				Value[ii].counter_B_R2 = Row[t].SS_B2.size();
				//Value[ii].val = (opts.globalK + 1);
			}
			else { // H1[tt] is an end point (e2)
				Value[ii].SS_A_R2 = Row[t].SS_A2;
				Value[ii].counter_A_R2 = Row[t].SS_A2.size();
				//Value[ii].val = (opts.globalK + 1);				
			}
		}
	}


	// Get SS_A_C1, SS_B_C1, SS_A_C2 and SS_B_C2 for each fragment
	for (unsigned int t = 0; t < Col.size(); ++t) {
		for (unsigned int tt = Col[t].pstart; tt < Col[t].pend; ++tt) { 

			unsigned int ii = H1[H2[tt]].frag_num;

			if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 1) { //H1[H2[tt]] a start point (s1)
				Value[ii].SS_B_C1 = Col[t].SS_B1;
				Value[ii].counter_B_C1 = Col[t].SS_B1.size();
				Value[ii].val = (opts.globalK + 1);
			}
			else if (H1[H2[tt]].ind == 0 and H1[H2[tt]].inv == 1) { // H1[H2[tt]] is an end point (e1)
				Value[ii].SS_A_C1 = Col[t].SS_A1;
				Value[ii].counter_A_C1 = Col[t].SS_A1.size();
				//Value[ii].val = (opts.globalK + 1);
			}
			else if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 0) { //H1[H2[tt]] a start point (s2)
				Value[ii].SS_B_C2 = Col[t].SS_B2;
				Value[ii].counter_B_C2 = Col[t].SS_B2.size();
				//Value[ii].val = (opts.globalK + 1);				
			}
			else { // H1[H2[tt]] is an end point (e2)
				Value[ii].SS_A_C2 = Col[t].SS_A2;
				Value[ii].counter_A_C2 = Col[t].SS_A2.size();
				//Value[ii].val = (opts.globalK + 1);				
			}

		}
	}


	//cerr << "Value: " << Value << endl;
	

	//cerr << "ProcessPoint\n";

	ProcessPoint(H1, Row, SubR1, SubC1, SubR2, SubC2, Value, opts, LookUpTable, rate);

	
	//cerr << "end\n";

	// find the max_value for the FinalChain 
	unsigned int  l = 0;
	float max_value = 0;
	unsigned int max_pos = 0;
	while (l < Value.size()) {
		if (Value[l].val > max_value) {
			max_value = Value[l].val;
			max_pos = l;
		}
		++l;
	}
	
	//cerr << "TraceBack\n";
	// Trace back to get the FinalChain
	// store the index of points
	chain.clear();
	TraceBack(SubR1, SubC1, SubR2,SubC2, Value, max_pos, chain);
	std::reverse(chain.begin(), chain.end());

	// Clear SubR and SubC
	SubR1.Clear(eeR1);
	SubC1.Clear(eeC1);
	SubR2.Clear(eeR2);
	SubC2.Clear(eeC2);

	// get the time for the program
	clock_t end = std::clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//cerr << "Time: " << elapsed_secs << endl;

/*	// print out the final chain
	cerr << "FinalChain: (";
	for (unsigned int i = 0; i < FinalChain.size(); ++i) {
		cerr << A[FinalChain[i]] << ", ";
	}
	cerr << ")" << endl;
*/

	return 0;
}


#endif