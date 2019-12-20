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
#include "Read.h"
#include "Chain.h"


using std::cerr;
using std::cout;
using std::endl;
using std::iota;


void 
ComputeMatchStart(vector<unsigned int> & MatchNum, int & totalMatch, const vector<Cluster> & FragInput, SplitChain & inputChain){
	totalMatch = FragInput[inputChain[0]].matches.size();
	for (int c = 1; c < inputChain.size(); c++) {
		totalMatch += FragInput[inputChain[c]].matches.size();
		MatchNum[c] = FragInput[inputChain[c - 1]].matches.size();
	}

	for (int c = 1; c < MatchNum.size(); c++) {
		MatchNum[c] = MatchNum[c-1] + MatchNum[c]; 
	}
}


void
insertPointsPair(const vector<Cluster> & FragInput, vector<Point> & H1, unsigned int frag_num, GenomePos & qs, GenomePos & ts, 
									int & length, int & clusterNum, int & matchstartNum, int Pair, int strand) {
	if (Pair == 0) {
		// 
		// Insert pair s1 and e1;
		//
		Point s1;
		H1.push_back(s1);
		H1.back().ind = 1; // start
		H1.back().inv = 1; // forward direction
		H1.back().frag_num = frag_num;
		H1.back().se.first = qs;
		H1.back().se.second = ts;										
		H1.back().orient = strand; // forward strand
		H1.back().clusterNum = clusterNum; 
		H1.back().matchstartNum = matchstartNum;

		// insert end point e1 into H1
		Point e1;
		H1.push_back(e1);
		H1.back().ind = 0; // end
		H1.back().inv = 1; // forward direction		
		H1.back().frag_num = frag_num;
		H1.back().se.first = qs + length;
		H1.back().se.second = ts + length;
		H1.back().orient = strand; 	
		H1.back().clusterNum = clusterNum; 
		H1.back().matchstartNum = matchstartNum;		
	}
	else {
		// 
		// Insert pair s2 and e2;
		//
		Point s2;
		H1.push_back(s2);
		H1.back().ind = 1; // start
		H1.back().inv = 0; // backward direction
		H1.back().frag_num = frag_num;
		H1.back().se.first = qs; 
		H1.back().se.second = ts + length - 1; //// TODO(Jingwen): fix this "-1" in the other SparseDP
		H1.back().orient = strand; // reversed strand			
		H1.back().clusterNum = clusterNum; 
		H1.back().matchstartNum = matchstartNum;	

		// insert end point e2 into H1
		Point e2;
		H1.push_back(e2);
		H1.back().ind = 0; // end
		H1.back().inv = 0; // backward direction		
		H1.back().frag_num = frag_num;
		H1.back().se.first = qs + length;
		assert(ts != 0); // IF this happens, just delete all the "-1" in the section of inserting points;
		H1.back().se.second = ts - 1;
		H1.back().orient = strand; // reversed strand			
		H1.back().clusterNum = clusterNum; 
		H1.back().matchstartNum = matchstartNum;		
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
template<typename Tup>
void ProcessPoint (const std::vector<Point> & H1, std::vector<info> & V, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1,
				 StackOfSubProblems & SubR2, StackOfSubProblems & SubC2, std::vector<Fragment_Info> & Value, Options & opts, 
				 const std::vector<float> & LookUpTable, const Tup & FragInput, int rate) { // std::vector<ClusterCoordinates> & FragInput

	for (unsigned int i = 0; i < H1.size(); i++) { // process points by row

		long int ForwardDiag = static_cast<long int>(H1[i].se.second) - static_cast<long int>(H1[i].se.first);
		long int BackDiag = static_cast<long int>(H1[i].se.second) + static_cast<long int>(H1[i].se.first);

		//cerr << "\n\n\n\nprocessing point " << H1[H3[i]] << endl;
		//cerr << "the Value[" << H1[H3[i]].frag_num << "] of this point:  " <<  Value[H1[H3[i]].frag_num] << "\n";
		
		unsigned int ii = H1[i].frag_num;


		if (H1[i].ind == 1 and H1[i].inv == 1) { // H1[i] is a start point s1
			//cerr << "----------------------------------------this is a start point (s1) ------------------------------------" << endl;
			//cerr << "dealing with the subproblem B first. Solve Value[" << H1[i].frag_num << "].SS_B_R1 and SS_B_C1" << endl;
			//cerr <<"----------------Solve SS_B_R1 -------------------" << endl;
			//
			// For each subproblem B_R1 that point H1[i] is in, get Ev[ForwardDiag] and update Value[H1[i]].frag_num].val
			//
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
						//cerr << "SubR1[" << j << "] is a non-leaf case but there is no forward diag smaller than the current in D array" << endl;
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
						

						SubR1[j].Ev[i1] = SubR1[j].Dv[i2] + w(SubR1[j].Di[i2], SubR1[j].Ei[i1], LookUpTable, opts) + FragInput.matchesLengths[ii] * rate; 
						SubR1[j].Ep[i1] = i2;							

						//cerr << "SubR1[" << j << "].Ev[" << i1 << "]: " << SubR1[j].Ev[i1] << ", SubR1[" << j << "].Ep[" << i1 << "]: " << SubR1[j].Ep[i1] << "\n"; 


						// Update the value of this point
						//TODO(Jingwen): if this point is a s1 of a reverse orientated anchor, only update the value when  SubR1[j].Dp[Ep[i1]] points to a forward orientated anchor

						//TODO(Jingwen): only for debug
						assert(Value[ii].orient == H1[i].orient);
						int p = SubR1[j].Dp[SubR1[j].Ep[i1]];
						//
						// If the current anchor is reverse oriented and the previous anchor is also reverse oriented, then do not update Value[ii]
						// Otherwise update Value[ii]
						//
						//if (! (Value[ii].orient == 0 and Value[p].orient == 0)) { 
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
						//}	
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
			
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubC1[j].S_1, SubC1[j].Ei, SubC1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					
	//					SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts) + opts.globalK; 
	//					SubC1[j].Ep[i1] = i2;							

						SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts) + FragInput.matchesLengths[ii] * rate; 
						SubC1[j].Ep[i1] = i2;							

						//cerr << "SubC1[" << j << "].Ev[" << i1 << "]: " << SubC1[j].Ev[i1] << ", SubC1[" << j << "].Ep[" << i1 << "]: " << SubC1[j].Ep[i1] << "\n"; 

						// Update the value of this point
						//TODO(Jingwen): if this point is a s1 of a reverse orientated anchor, only update the value when  SubR1[j].Dp[Ep[i1]] points to a forward orientated anchor
						//TODO(Jingwen): Only for debug
						assert(Value[ii].orient == H1[i].orient);
						int p = SubC1[j].Dp[SubC1[j].Ep[i1]];
						//
						// If the current anchor is reverse oriented and the previous anchor is also reverse oriented, then do not update Value[ii]
						// Otherwise update Value[ii]
						//if (! (Value[ii].orient == 0 and Value[p].orient == 0)) {
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
						//}
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

						SubR2[j].Ev[i1] = SubR2[j].Dv[i2] + w(SubR2[j].Di[i2], SubR2[j].Ei[i1], LookUpTable, opts) + FragInput.matchesLengths[ii] * rate; 
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
					

						SubC2[j].Ev[i1] = SubC2[j].Dv[i2] + w(SubC2[j].Di[i2], SubC2[j].Ei[i1], LookUpTable, opts) + FragInput.matchesLengths[ii] * rate; 
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
*/

		
//
// Each fragment has different length, so we need to pass FragInput to this function
// This function is for splitCluster
//
template<typename Tup>
void ProcessPoint (const std::vector<Point> & H1, std::vector<info> & V, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1,
				 StackOfSubProblems & SubR2, StackOfSubProblems & SubC2, std::vector<Fragment_Info> & Value, Options & opts, 
				 const std::vector<float> & LookUpTable, const std::vector<Tup> & FragInput, float & rate, bool step) { // std::vector<ClusterCoordinates> & FragInput

	for (unsigned int i = 0; i < H1.size(); ++i) { // process points by row

		long int ForwardDiag = static_cast<long int>(H1[i].se.second) - static_cast<long int>(H1[i].se.first);
		long int BackDiag = static_cast<long int>(H1[i].se.second) + static_cast<long int>(H1[i].se.first);

		//cerr << "\n\n\n\nprocessing point " << H1[H3[i]] << endl;
		//cerr << "the Value[" << H1[H3[i]].frag_num << "] of this point:  " <<  Value[H1[H3[i]].frag_num] << "\n";
		
		unsigned int ii = H1[i].frag_num;


		if (H1[i].ind == 1 and H1[i].inv == 1) { // H1[i] is a start point s1
			//cerr << "----------------------------------------this is a start point (s1) ------------------------------------" << endl;
			//cerr << "dealing with the subproblem B first. Solve Value[" << H1[i].frag_num << "].SS_B_R1 and SS_B_C1" << endl;
			//cerr <<"----------------Solve SS_B_R1 -------------------" << endl;
			//
			// For each subproblem B_R1 that point H1[i] is in, get Ev[ForwardDiag] and update Value[H1[i]].frag_num].val
			//
			for (unsigned int k = 0; k < Value[ii].SS_B_R1.size(); ++k) {
				unsigned int j = Value[ii].SS_B_R1[Value[ii].SS_B_R1.size() - 1 - k];
				//cerr << "Solving SubR1[" << j << "]: " << SubR1[j]<< "\n";

				// If subproblem SubR1[j] is a leaf case, then 
				if (SubR1[j].Di.empty()) {
					//cerr << "SubR1[" << j << "] is a leaf case or non-leaf case but Di is empty" << "\n";
					assert(Value[ii].counter_B_R1 != 0);
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
						//cerr << "SubR1[" << j << "] is a non-leaf case but there is no forward diag smaller than the current in D array" << endl;
						assert(Value[ii].counter_B_R1 != 0);
						--Value[ii].counter_B_R1;
						continue;
					}
					else {
						//assert(SubR1[j].counter_D[SubR1[j].Eb[*t]] == 0);
						//cerr << "SubR1[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already." << endl;
						//cerr << "Start to compute the maximization structure.\n";

						SubR1[j].now = SubR1[j].Eb[*t];
						Maximization (SubR1[j].now, SubR1[j].last, SubR1[j].Di, SubR1[j].Ei, SubR1[j].Dv, SubR1[j].Db, SubR1[j].Block, SubR1[j].S_1, LookUpTable, opts, step); // TODO(Jingwen) anything change for SubC????
						SubR1[j].last = SubR1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubR1[j].S_1, SubR1[j].Ei, SubR1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
						

						SubR1[j].Ev[i1] = SubR1[j].Dv[i2] + w(SubR1[j].Di[i2], SubR1[j].Ei[i1], LookUpTable, opts, step) + 
																std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
						SubR1[j].Ep[i1] = i2;							

						//cerr << "SubR1[" << j << "].Ev[" << i1 << "]: " << SubR1[j].Ev[i1] << ", SubR1[" << j << "].Ep[" << i1 << "]: " << SubR1[j].Ep[i1] << "\n"; 


						// Update the value of this point
						//TODO(Jingwen): if this point is a s1 of a reverse orientated anchor, only update the value when  SubR1[j].Dp[Ep[i1]] points to a forward orientated anchor

						//TODO(Jingwen): only for debug
						assert(Value[ii].orient == H1[i].orient);
						int p = SubR1[j].Dp[SubR1[j].Ep[i1]];
						// if the current anchor is reverse oriented and the previous anchor is also reverse oriented, then do not update Value[ii]
						// Otherwise update Value[ii]
						//if (! (Value[ii].orient == 0 and Value[p].orient == 0)) { 
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
						//}	
						assert(Value[ii].counter_B_R1 != 0);
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
					assert(Value[ii].counter_B_C1 != 0);
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
						assert(Value[ii].counter_B_C1 != 0);
						--Value[ii].counter_B_C1;
						continue;
					}
					else {

						//assert(SubC1[j].counter_D[SubC1[j].Eb[*t]] == 0);
						//cerr << "SubC1[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already.\n";
						//cerr << "Start to compute the maximization structure.\n";
						SubC1[j].now = SubC1[j].Eb[*t]; 
						Maximization (SubC1[j].now, SubC1[j].last, SubC1[j].Di, SubC1[j].Ei, SubC1[j].Dv, SubC1[j].Db, SubC1[j].Block, SubC1[j].S_1, LookUpTable, opts, step); 
						SubC1[j].last = SubC1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubC1[j].S_1, SubC1[j].Ei, SubC1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					
	//					SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts) + opts.globalK; 
	//					SubC1[j].Ep[i1] = i2;							

						SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts, step) + 
																std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
						SubC1[j].Ep[i1] = i2;							

						//cerr << "SubC1[" << j << "].Ev[" << i1 << "]: " << SubC1[j].Ev[i1] << ", SubC1[" << j << "].Ep[" << i1 << "]: " << SubC1[j].Ep[i1] << "\n"; 

						// Update the value of this point
						//TODO(Jingwen): if this point is a s1 of a reverse orientated anchor, only update the value when  SubR1[j].Dp[Ep[i1]] points to a forward orientated anchor
						//TODO(Jingwen): Only for debug
						assert(Value[ii].orient == H1[i].orient);
						int p = SubC1[j].Dp[SubC1[j].Ep[i1]];
						//
						// If the current anchor is reverse oriented and the previous anchor is also reverse oriented, then do not update Value[ii]
						// Otherwise update Value[ii]
						//if (! (Value[ii].orient == 0 and Value[p].orient == 0)) {
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
						//}
						assert(Value[ii].counter_B_C1 != 0);
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
					assert(Value[ii].counter_B_R2 != 0);
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
						assert(Value[ii].counter_B_R2 != 0);
						--Value[ii].counter_B_R2;
						continue;
					}
					else {
						//assert(SubR2[j].counter_D[SubR2[j].Eb[*t]] == 0);
						//cerr << "SubR2[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already." << endl;
						//cerr << "Start to compute the maximization structure.\n";

						SubR2[j].now = SubR2[j].Eb[*t];
						Maximization (SubR2[j].now, SubR2[j].last, SubR2[j].Di, SubR2[j].Ei, SubR2[j].Dv, SubR2[j].Db, SubR2[j].Block, SubR2[j].S_1, LookUpTable, opts, step); // TODO(Jingwen) anything change for SubC????
						SubR2[j].last = SubR2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that BackDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for BackDiag
						FindValueInBlock(BackDiag, SubR2[j].S_1, SubR2[j].Ei, SubR2[j].Block, i1, i2);
						//cerr << "the index in Ei that BackDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for BackDiag ---- i2: " << i2 << "\n";

						SubR2[j].Ev[i1] = SubR2[j].Dv[i2] + w(SubR2[j].Di[i2], SubR2[j].Ei[i1], LookUpTable, opts, step) + 
																	std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
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
						assert(Value[ii].counter_B_R2 != 0);
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
					assert(Value[ii].counter_B_C2 != 0);
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
						assert(Value[ii].counter_B_C2 != 0);
						--Value[ii].counter_B_C2;
						continue;
					}
					else {

						//assert(SubC2[j].counter_D[SubC2[j].Eb[*t]] == 0);
						//cerr << "SubC2[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already.\n";
						//cerr << "Start to compute the maximization structure.\n";
						SubC2[j].now = SubC2[j].Eb[*t]; 
						Maximization (SubC2[j].now, SubC2[j].last, SubC2[j].Di, SubC2[j].Ei, SubC2[j].Dv, SubC2[j].Db, SubC2[j].Block, SubC2[j].S_1, LookUpTable, opts, step); 
						SubC2[j].last = SubC2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores 
						unsigned int i2;
						FindValueInBlock(BackDiag, SubC2[j].S_1, SubC2[j].Ei, SubC2[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					

						SubC2[j].Ev[i1] = SubC2[j].Dv[i2] + w(SubC2[j].Di[i2], SubC2[j].Ei[i1], LookUpTable, opts, step) + 
																	std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
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
						assert(Value[ii].counter_B_C2 != 0);	
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


//
// Each fragment has different length, so we need to pass vector<Cluster> vt to the function
// This function is for ExtendClusters;
//
template<typename Tup>
void ProcessPoint (const std::vector<Point> & H1, std::vector<info> & V, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1,
				 StackOfSubProblems & SubR2, StackOfSubProblems & SubC2, std::vector<Fragment_Info> & Value, Options & opts, 
				 const std::vector<float> & LookUpTable, const std::vector<Tup> & FragInput, SplitChain & inputChain, 
				 	const vector<unsigned int> & MatchStart, float & rate, bool step) { // std::vector<ClusterCoordinates> & FragInput

	for (unsigned int i = 0; i < H1.size(); ++i) { // process points by row

		long int ForwardDiag = static_cast<long int>(H1[i].se.second) - static_cast<long int>(H1[i].se.first);
		long int BackDiag = static_cast<long int>(H1[i].se.second) + static_cast<long int>(H1[i].se.first);

		//cerr << "\n\n\n\nprocessing point " << H1[H3[i]] << endl;
		//cerr << "the Value[" << H1[H3[i]].frag_num << "] of this point:  " <<  Value[H1[H3[i]].frag_num] << "\n";
		
		unsigned int ii = H1[i].frag_num;
		int mi = H1[i].matchstartNum;


		if (H1[i].ind == 1 and H1[i].inv == 1) { // H1[i] is a start point s1
			//cerr << "----------------------------------------this is a start point (s1) ------------------------------------" << endl;
			//cerr << "dealing with the subproblem B first. Solve Value[" << H1[i].frag_num << "].SS_B_R1 and SS_B_C1" << endl;
			//cerr <<"----------------Solve SS_B_R1 -------------------" << endl;
			//
			// For each subproblem B_R1 that point H1[i] is in, get Ev[ForwardDiag] and update Value[H1[i]].frag_num].val
			//
			for (unsigned int k = 0; k < Value[ii].SS_B_R1.size(); ++k) {
				unsigned int j = Value[ii].SS_B_R1[Value[ii].SS_B_R1.size() - 1 - k];
				//cerr << "Solving SubR1[" << j << "]: " << SubR1[j]<< "\n";

				// If subproblem SubR1[j] is a leaf case, then 
				if (SubR1[j].Di.empty()) {
					//cerr << "SubR1[" << j << "] is a leaf case or non-leaf case but Di is empty" << "\n";
					assert(Value[ii].counter_B_R1 != 0);
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
						//cerr << "SubR1[" << j << "] is a non-leaf case but there is no forward diag smaller than the current in D array" << endl;
						assert(Value[ii].counter_B_R1 != 0);
						--Value[ii].counter_B_R1;
						continue;
					}
					else {
						//assert(SubR1[j].counter_D[SubR1[j].Eb[*t]] == 0);
						//cerr << "SubR1[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already." << endl;
						//cerr << "Start to compute the maximization structure.\n";

						SubR1[j].now = SubR1[j].Eb[*t];
						Maximization (SubR1[j].now, SubR1[j].last, SubR1[j].Di, SubR1[j].Ei, SubR1[j].Dv, SubR1[j].Db, SubR1[j].Block, SubR1[j].S_1, LookUpTable, opts, step); // TODO(Jingwen) anything change for SubC????
						SubR1[j].last = SubR1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubR1[j].S_1, SubR1[j].Ei, SubR1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
						

						SubR1[j].Ev[i1] = SubR1[j].Dv[i2] + w(SubR1[j].Di[i2], SubR1[j].Ei[i1], LookUpTable, opts, step) + 
																rate * FragInput[inputChain[mi]].matchesLengths[ii - MatchStart[mi]];
																//std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
						SubR1[j].Ep[i1] = i2;							

						//cerr << "SubR1[" << j << "].Ev[" << i1 << "]: " << SubR1[j].Ev[i1] << ", SubR1[" << j << "].Ep[" << i1 << "]: " << SubR1[j].Ep[i1] << "\n"; 


						// Update the value of this point
						//TODO(Jingwen): if this point is a s1 of a reverse orientated anchor, only update the value when  SubR1[j].Dp[Ep[i1]] points to a forward orientated anchor

						//TODO(Jingwen): only for debug
						assert(Value[ii].orient == H1[i].orient);
						int p = SubR1[j].Dp[SubR1[j].Ep[i1]];
						// if the current anchor is reverse oriented and the previous anchor is also reverse oriented, then do not update Value[ii]
						// Otherwise update Value[ii]
						//if (! (Value[ii].orient == 0 and Value[p].orient == 0)) { 
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
						//}	
						assert(Value[ii].counter_B_R1 != 0);
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
					assert(Value[ii].counter_B_C1 != 0);
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
						assert(Value[ii].counter_B_C1 != 0);
						--Value[ii].counter_B_C1;
						continue;
					}
					else {

						//assert(SubC1[j].counter_D[SubC1[j].Eb[*t]] == 0);
						//cerr << "SubC1[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already.\n";
						//cerr << "Start to compute the maximization structure.\n";
						SubC1[j].now = SubC1[j].Eb[*t]; 
						Maximization (SubC1[j].now, SubC1[j].last, SubC1[j].Di, SubC1[j].Ei, SubC1[j].Dv, SubC1[j].Db, SubC1[j].Block, SubC1[j].S_1, LookUpTable, opts, step); 
						SubC1[j].last = SubC1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubC1[j].S_1, SubC1[j].Ei, SubC1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					
	//					SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts) + opts.globalK; 
	//					SubC1[j].Ep[i1] = i2;							

						SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts, step) + 
													rate * FragInput[inputChain[mi]].matchesLengths[ii - MatchStart[mi]];
																//std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
						SubC1[j].Ep[i1] = i2;							

						//cerr << "SubC1[" << j << "].Ev[" << i1 << "]: " << SubC1[j].Ev[i1] << ", SubC1[" << j << "].Ep[" << i1 << "]: " << SubC1[j].Ep[i1] << "\n"; 

						// Update the value of this point
						//TODO(Jingwen): if this point is a s1 of a reverse orientated anchor, only update the value when  SubR1[j].Dp[Ep[i1]] points to a forward orientated anchor
						//TODO(Jingwen): Only for debug
						assert(Value[ii].orient == H1[i].orient);
						int p = SubC1[j].Dp[SubC1[j].Ep[i1]];
						//
						// If the current anchor is reverse oriented and the previous anchor is also reverse oriented, then do not update Value[ii]
						// Otherwise update Value[ii]
						//if (! (Value[ii].orient == 0 and Value[p].orient == 0)) {
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
						//}
						assert(Value[ii].counter_B_C1 != 0);
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
					assert(Value[ii].counter_B_R2 > 0);
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
						assert(Value[ii].counter_B_R2 > 0);
						--Value[ii].counter_B_R2;
						continue;
					}
					else {
						//assert(SubR2[j].counter_D[SubR2[j].Eb[*t]] == 0);
						//cerr << "SubR2[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already." << endl;
						//cerr << "Start to compute the maximization structure.\n";

						SubR2[j].now = SubR2[j].Eb[*t];
						Maximization (SubR2[j].now, SubR2[j].last, SubR2[j].Di, SubR2[j].Ei, SubR2[j].Dv, SubR2[j].Db, SubR2[j].Block, SubR2[j].S_1, LookUpTable, opts, step); // TODO(Jingwen) anything change for SubC????
						SubR2[j].last = SubR2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that BackDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for BackDiag
						FindValueInBlock(BackDiag, SubR2[j].S_1, SubR2[j].Ei, SubR2[j].Block, i1, i2);
						//cerr << "the index in Ei that BackDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for BackDiag ---- i2: " << i2 << "\n";

						SubR2[j].Ev[i1] = SubR2[j].Dv[i2] + w(SubR2[j].Di[i2], SubR2[j].Ei[i1], LookUpTable, opts, step) + 
														rate * FragInput[inputChain[mi]].matchesLengths[ii - MatchStart[mi]];
																	//std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
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
						assert(Value[ii].counter_B_R2 > 0);
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
					assert(Value[ii].counter_B_C2 > 0);
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
						assert(Value[ii].counter_B_C2 > 0);
						--Value[ii].counter_B_C2;
						continue;
					}
					else {

						//assert(SubC2[j].counter_D[SubC2[j].Eb[*t]] == 0);
						//cerr << "SubC2[" << j << "] is a non-leaf case and The part of D array where forward diags are smaller than current is filled out already.\n";
						//cerr << "Start to compute the maximization structure.\n";
						SubC2[j].now = SubC2[j].Eb[*t]; 
						Maximization (SubC2[j].now, SubC2[j].last, SubC2[j].Di, SubC2[j].Ei, SubC2[j].Dv, SubC2[j].Db, SubC2[j].Block, SubC2[j].S_1, LookUpTable, opts, step); 
						SubC2[j].last = SubC2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores 
						unsigned int i2;
						FindValueInBlock(BackDiag, SubC2[j].S_1, SubC2[j].Ei, SubC2[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					

						SubC2[j].Ev[i1] = SubC2[j].Dv[i2] + w(SubC2[j].Di[i2], SubC2[j].Ei[i1], LookUpTable, opts, step) + 
														rate * FragInput[inputChain[mi]].matchesLengths[ii - MatchStart[mi]];
																	//std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
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
						assert(Value[ii].counter_B_C2 > 0);
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


//
// This function is for tracing back such a chain that every anchor on this chain is unused;
//
void 
TraceBack (StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
					const vector<Fragment_Info> & Value, unsigned int & i, vector<unsigned int> & onechain, vector<bool> & link, vector<bool> & used) {

	long int prev_sub = Value[i].prev_sub;
	long int prev_ind = Value[i].prev_ind;
	if (used[i] == 0) {
		onechain.push_back(i);
		used[i] = 1;
	
		if (prev_sub != -1 and prev_ind != -1) {

			if (Value[i].prev == 1 and Value[i].inv == 1) { // The previous subproblem is SubR1
				unsigned int ind = SubR1[prev_sub].Ep[prev_ind];
				if (used[SubR1[prev_sub].Dp[ind]] == 0) {
					link.push_back(0);
					TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR1[prev_sub].Dp[ind], onechain, link, used);					
				}
				else {
					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
						used[onechain[lu]] = 0;
					}
					onechain.clear();
					link.clear();
				}
			}
			else if (Value[i].prev == 1 and Value[i].inv == 0) { // The previous subproblem is SubR2
				unsigned int ind = SubR2[prev_sub].Ep[prev_ind];
				if (used[SubR2[prev_sub].Dp[ind]] == 0) {
					link.push_back(1);
					TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR2[prev_sub].Dp[ind], onechain, link, used);
				}
				else {
					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
						used[onechain[lu]] = 0;
					}
					onechain.clear();	
					link.clear();				
				}
			}
			else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
				unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
				if (used[SubC1[prev_sub].Dp[ind]] == 0) {
					link.push_back(0);
					TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC1[prev_sub].Dp[ind], onechain, link, used);
				}
				else {
					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
						used[onechain[lu]] = 0;
					}
					onechain.clear();	
					link.clear();					
				}
			}		
			else { // The previous subproblem is SubC2
				unsigned int ind = SubC2[prev_sub].Ep[prev_ind];
				if (used[SubC2[prev_sub].Dp[ind]] == 0) {
					link.push_back(1);
					TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC2[prev_sub].Dp[ind], onechain, link, used);
				}
				else {
					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
						used[onechain[lu]] = 0;
					}
					onechain.clear();	
					link.clear();				
				}
			}

		}
	}
}

/*
//
// This function is for tracing back a chain;
//
void 
TraceBack (StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
					const vector<Fragment_Info> & Value, unsigned int & i, vector<unsigned int> & Chain) {

	long int prev_sub = Value[i].prev_sub;
	long int prev_ind = Value[i].prev_ind;
	Chain.push_back(i);
	//cerr << "i: " << i << endl;

	if (prev_sub != -1 and prev_ind != -1) {
		
		if (Value[i].prev == 1 and Value[i].inv == 1) { // The previous subproblem is SubR1
			assert(prev_sub <  SubR1.StackSub.size());
			assert(prev_ind < SubR1[prev_sub].Ep.size());
			unsigned int ind = SubR1[prev_sub].Ep[prev_ind];
			assert(ind < SubR1[prev_sub].Dp.size());
			TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR1[prev_sub].Dp[ind], Chain);
		}
		else if (Value[i].prev == 1 and Value[i].inv == 0) { // The previous subproblem is SubR2
			assert(prev_sub <  SubR2.StackSub.size());
			assert(prev_ind < SubR2[prev_sub].Ep.size());
			unsigned int ind = SubR2[prev_sub].Ep[prev_ind];
			assert(ind < SubR2[prev_sub].Dp.size());
			TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR2[prev_sub].Dp[ind], Chain);
		}
		else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
			assert(prev_sub <  SubC1.StackSub.size());
			assert(prev_ind < SubC1[prev_sub].Ep.size());
			unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
			assert(ind < SubC1[prev_sub].Dp.size());
			TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC1[prev_sub].Dp[ind], Chain);
		}		
		else { // The previous subproblem is SubC2
			assert(prev_sub < SubC2.StackSub.size());
			assert(prev_ind < SubC2[prev_sub].Ep.size());
			unsigned int ind = SubC2[prev_sub].Ep[prev_ind];
			assert(ind < SubC2[prev_sub].Dp.size());
			TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC2[prev_sub].Dp[ind], Chain);
		}

	}
}
*/


//
// This function is for tracing back a chain;
//
void 
TraceBack (StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
					const vector<Fragment_Info> & Value, unsigned int i, vector<unsigned int> & Chain) {

	long int prev_sub = Value[i].prev_sub;
	long int prev_ind = Value[i].prev_ind;
	Chain.push_back(i);
	//cerr << "i: " << i << endl;

	while (prev_sub != -1 and prev_ind != -1) {
		
		if (Value[i].prev == 1 and Value[i].inv == 1) { // The previous subproblem is SubR1
			assert(prev_sub <  SubR1.StackSub.size());
			assert(prev_ind < SubR1[prev_sub].Ep.size());
			unsigned int ind = SubR1[prev_sub].Ep[prev_ind];
			assert(ind < SubR1[prev_sub].Dp.size());
			i = SubR1[prev_sub].Dp[ind];
			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR1[prev_sub].Dp[ind], Chain);
		}
		else if (Value[i].prev == 1 and Value[i].inv == 0) { // The previous subproblem is SubR2
			assert(prev_sub <  SubR2.StackSub.size());
			assert(prev_ind < SubR2[prev_sub].Ep.size());
			unsigned int ind = SubR2[prev_sub].Ep[prev_ind];
			assert(ind < SubR2[prev_sub].Dp.size());
			i = SubR2[prev_sub].Dp[ind];
			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR2[prev_sub].Dp[ind], Chain);
		}
		else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
			assert(prev_sub <  SubC1.StackSub.size());
			assert(prev_ind < SubC1[prev_sub].Ep.size());
			unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
			assert(ind < SubC1[prev_sub].Dp.size());
			i = SubC1[prev_sub].Dp[ind];
			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC1[prev_sub].Dp[ind], Chain);
		}		
		else { // The previous subproblem is SubC2
			assert(prev_sub < SubC2.StackSub.size());
			assert(prev_ind < SubC2[prev_sub].Ep.size());
			unsigned int ind = SubC2[prev_sub].Ep[prev_ind];
			assert(ind < SubC2[prev_sub].Dp.size());
			i = SubC2[prev_sub].Dp[ind];
			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC2[prev_sub].Dp[ind], Chain);
		}
		prev_sub = Value[i].prev_sub;
		prev_ind = Value[i].prev_ind;		
		Chain.push_back(i);
	}
}

//
// This function is for splitClusters;
//
void 
DecidePrimaryChains(const vector<Cluster> & FragInput, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
					const std::vector<Fragment_Info> & Value, vector<Primary_chain> &Primary_chains, Read & read, Options & opts) {

	std::vector<bool> used(Value.size(), 0);
	Fragment_valueOrder fragments_valueOrder(&Value);

	for (unsigned int fv = 0; fv < fragments_valueOrder.size(); fv++) {

		unsigned int d = fragments_valueOrder.index[fv];

		vector<unsigned int> onechain;
		vector<bool> link;
		TraceBack(SubR1, SubC1, SubR2, SubC2, Value, d, onechain, link, used);

		if (onechain.size() != 0) {
			// Note: onechain store index from the last one to the first one
			GenomePos qEnd = FragInput[onechain[0]].qEnd;
			GenomePos tEnd = FragInput[onechain[0]].tEnd;
			GenomePos qStart = FragInput[onechain.back()].qStart;
			GenomePos tStart = FragInput[onechain.back()].tStart;	

			// Somethimes not all onechain[i] align to the same chromosom
			for (int c = 0; c < onechain.size(); c++) {
				qEnd = max(FragInput[onechain[c]].qEnd, qEnd);
				tEnd = max(FragInput[onechain[c]].tEnd, tEnd);
				qStart = min(FragInput[onechain[c]].qStart, qStart);
				tStart = min(FragInput[onechain[c]].tStart, tStart);
			}

			//
			// If this chain overlap with read greater than 10%, insert it to chains
			//
			if (((float)(qEnd - qStart)/read.length) > 0.1) {
				//
				// Compare onechain to all the primary chains we've found. 
				// If onechain overlaps with one primary chain over 70% ---> onechain is a secondary chain 
				// If onechain overlaps with all the primary chains less than 50% ---> onechain is another primary chain
				//
				if (Primary_chains.size() == 0) {
					Primary_chain Pc(CHain(qStart, qEnd, tStart, tEnd, onechain, link));
					Primary_chains.push_back(Pc);
				} 
				else {
					bool newpr = 1, inserted = 0;
					int p = 0;
					while (p < Primary_chains.size()) {
						if (Primary_chains[p].chains[0].Overlaps(qStart, qEnd, 0.7)) {
							if (Primary_chains[p].chains.size() < opts.SecondaryAln + 1) {
								Primary_chains[p].chains.push_back(CHain(qStart, qEnd, tStart, tEnd, onechain, link));
								inserted = 1;								
							}
							break;
						}
						else if (Primary_chains[p].chains[0].Overlaps(qStart, qEnd, 0.5)) {
							newpr = 0;
						}
						++p;
					}			
					if (p == Primary_chains.size() - 1 and inserted == 0 and newpr == 1) {		
						if (Primary_chains.size() < opts.PrimaryAln) {
							Primary_chain Pc(CHain(qStart, qEnd, tStart, tEnd, onechain, link));
							Primary_chains.push_back(Pc);
						}	
						else break;		
					}	
				}
			}
			else break;				
		}
		onechain.clear();
		link.clear();
	}
}

/*
// (KEEPKEEP)
// The input for this function is GenomePairs which is NOT from Merge_Split step
// Each fragment has the same length ----> change to different length
int SparseDP (const Cluster &FragInput, std::vector<unsigned int> &chain, Options &opts, const std::vector<float> &LookUpTable, bool ReverseOnly = 1, int rate = 1) {

	if (FragInput.matches.size() == 0) return 0;

	std::vector<Point> H1;
	// get points from FragInput and store them in H1			
	for (unsigned int i = 0; i < FragInput.matches.size(); i++) {

		if (FragInput.strand == 0) {
			// insert start point s1 into H1
			Point s1;
			H1.push_back(s1);
			H1.back().ind = 1; // start
			H1.back().inv = 1; // forward direction
			H1.back().frag_num = i;
			H1.back().se.first = FragInput.matches[i].first.pos;
			H1.back().se.second = FragInput.matches[i].second.pos;	
			H1.back().orient = 1; // the point comes from a forward oriented anchor

			// insert end point e1 into H1
			Point e1;
			H1.push_back(e1);
			H1.back().ind = 0; // end
			H1.back().inv = 1; // forward direction		
			H1.back().frag_num = i;
			H1.back().se.first = FragInput.matches[i].first.pos + FragInput.matchesLengths[i];
			H1.back().se.second = FragInput.matches[i].second.pos + FragInput.matchesLengths[i];	
			H1.back().orient = 1; // the point comes from a forward oriented anchor				
		}
		else if (ReverseOnly == 0) {
			// insert start point s1 into H1
			Point s1;
			H1.push_back(s1);
			H1.back().ind = 1; // start
			H1.back().inv = 1; // forward direction
			H1.back().frag_num = i;
			H1.back().se.first = FragInput.matches[i].first.pos;
			H1.back().se.second = FragInput.matches[i].second.pos;		
			H1.back().orient = 0; // the point comes from a reverse oriented anchor


			// insert end point e1 into H1
			Point e1;
			H1.push_back(e1);
			H1.back().ind = 0; // end
			H1.back().inv = 1; // forward direction		
			H1.back().frag_num = i;
			H1.back().se.first = FragInput.matches[i].first.pos + FragInput.matchesLengths[i];
			H1.back().se.second = FragInput.matches[i].second.pos + FragInput.matchesLengths[i];
			H1.back().orient = 0; // the point comes from a reverse oriented anchor


			// insert start point s2 into H1
			Point s2;
			H1.push_back(s2);
			H1.back().ind = 1; // start
			H1.back().inv = 0; // backward direction
			H1.back().frag_num = i;
			H1.back().se.first = FragInput.matches[i].first.pos;
			H1.back().se.second = FragInput.matches[i].second.pos + FragInput.matchesLengths[i]; 
			H1.back().orient = 0; // the point comes from a reverse oriented anchor


			// insert end point e2 into H1
			Point e2;
			H1.push_back(e2);
			H1.back().ind = 0; // end
			H1.back().inv = 0; // backward direction		
			H1.back().frag_num = i;
			H1.back().se.first = FragInput.matches[i].first.pos + FragInput.matchesLengths[i];
			H1.back().se.second = FragInput.matches[i].second.pos;	
			H1.back().orient = 0; // the point comes from a reverse oriented anchor
		} 
		else {
			// insert start point s2 into H1
			Point s2;
			H1.push_back(s2);
			H1.back().ind = 1; // start
			H1.back().inv = 0; // backward direction
			H1.back().frag_num = i;
			H1.back().se.first = FragInput.matches[i].first.pos;
			H1.back().se.second = FragInput.matches[i].second.pos + FragInput.matchesLengths[i]; 
			H1.back().orient = 0; // the point comes from a reverse oriented anchor


			// insert end point e2 into H1
			Point e2;
			H1.push_back(e2);
			H1.back().ind = 0; // end
			H1.back().inv = 0; // backward direction		
			H1.back().frag_num = i;
			H1.back().se.first = FragInput.matches[i].first.pos + FragInput.matchesLengths[i];
			H1.back().se.second = FragInput.matches[i].second.pos;	
			H1.back().orient = 0; // the point comes from a reverse oriented anchor			
		}
	}


	clock_t begin = std::clock();

	//Sort the point by row
	sort(H1.begin(), H1.end(), SortByRowOp<Point>()); // with same q and t coordinates, end point < start point

	//cerr << "H1: " << H1 << endl;
	std::vector<unsigned int> H2(H1.size());
	//std::vector<unsigned int> H3(H1.size()); // TODO(Jingwen): Probably don't need this
	iota(H2.begin(), H2.end(), 0);
	//iota(H3.begin(), H3.end(), 0);

	//Sort the point by column 
	sort(H2.begin(), H2.end(), SortByColOp<Point, unsigned int>(H1));
	//sort(H3.begin(), H3.end(), SortByBackDiagOp<Point, unsigned int>(H1));

	

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
	std::vector<Fragment_Info> Value(FragInput.matches.size());
	for (unsigned int t = 0; t < Row.size(); ++t) {
		for (unsigned int tt = Row[t].pstart; tt < Row[t].pend; ++tt) {
			//cerr << "Row H tt: " << tt << endl;
			unsigned int ii = H1[tt].frag_num;
			//cerr << "Row Value ii: " << ii << endl<< endl;

			if (H1[tt].ind == 1 and H1[tt].inv == 1) { //H1[tt] is a start point (s1)
				Value[ii].SS_B_R1 = Row[t].SS_B1;
				Value[ii].counter_B_R1 = Row[t].SS_B1.size();
				Value[ii].val = FragInput.matchesLengths[ii] * rate;
				Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 0 and H1[tt].inv == 1) { // H1[tt] is an end point (e1)
				Value[ii].SS_A_R1 = Row[t].SS_A1;
				Value[ii].counter_A_R1 = Row[t].SS_A1.size();
				//Value[ii].val = opts.globalK * rate; 
				//Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 1 and H1[tt].inv == 0) { //H1[tt] is a start point (s2)
				Value[ii].SS_B_R2 = Row[t].SS_B2;
				Value[ii].counter_B_R2 = Row[t].SS_B2.size();
				Value[ii].val = FragInput.matchesLengths[ii] * rate; //This is not redundant, because not every match has four points, some may have only two (s2, e2)
				Value[ii].orient = H1[tt].orient;
			}
			else { // H1[tt] is an end point (e2)
				Value[ii].SS_A_R2 = Row[t].SS_A2;
				Value[ii].counter_A_R2 = Row[t].SS_A2.size();
				//Value[ii].val = opts.globalK * rate;
				//Value[ii].orient = H1[tt].orient;			
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
				Value[ii].val = FragInput.matchesLengths[ii] * rate;
				Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 0 and H1[H2[tt]].inv == 1) { // H1[H2[tt]] is an end point (e1)
				Value[ii].SS_A_C1 = Col[t].SS_A1;
				Value[ii].counter_A_C1 = Col[t].SS_A1.size();
				//Value[ii].val = opts.globalK * rate;
				//Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 0) { //H1[H2[tt]] a start point (s2)
				Value[ii].SS_B_C2 = Col[t].SS_B2;
				Value[ii].counter_B_C2 = Col[t].SS_B2.size();
				Value[ii].val = FragInput.matchesLengths[ii] * rate;
				Value[ii].orient = H1[H2[tt]].orient;				
			}
			else { // H1[H2[tt]] is an end point (e2)
				Value[ii].SS_A_C2 = Col[t].SS_A2;
				Value[ii].counter_A_C2 = Col[t].SS_A2.size();
				//Value[ii].val = opts.globalK * rate;
				//Value[ii].orient = H1[H2[tt]].orient;			
			}

		}
	}


	//cerr << "Value: " << Value << endl;
	

	//cerr << "ProcessPoint\n";

	//ProcessPoint(H1, Row, SubR1, SubC1, SubR2, SubC2, Value, opts, LookUpTable, rate);
	ProcessPoint<Cluster>(H1, Row, SubR1, SubC1, SubR2, SubC2, Value, opts, LookUpTable, FragInput, rate);		

	
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
	TraceBack(SubR1, SubC1, SubR2, SubC2, Value, max_pos, chain);
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
	return 0;
}
*/


//
// The input for this function is splitchains[st] and ExtendClusters
// Tackle with fragments of different lengths;
// This SDP needs to insert 4 points for only anchors in the overlapping region between Clusters;
// This SDP needs to increase the cost for linking 2 anchors of different directions;
//
int SparseDP (SplitChain & inputChain, vector<Cluster> & FragInput, FinalChain & finalchain, Options & opts, const vector<float> & LookUpTable, Read & read) {

	//
	// Input: vector<Cluster> & FragInput and vector<unsigned int> inputChain;
	// Generally input 2 points for each anchor; s1 and e1 for forward anchor, while s2 and e2 for reversed anchor;
	// For anchors in the overlapping region between Clusters, we need to insert 4 points for them;
	//

	// Compute the matches start in each Cluster;
	//
	vector<unsigned int> MatchStart(inputChain.size(), 0);
	int totalMatch = 0;
	ComputeMatchStart(MatchStart, totalMatch, FragInput, inputChain);

	// Decide the rate;
	//
	float rate = 1;
	//if ((float)totalMatch / (float) read.length < 0.1) rate = 3;

	// get points from FragInput and store them in H1;
	//
	std::vector<Point> H1;		
	int next, prev;
	for (int c = 0; c < inputChain.size(); c++) {

		int cm = inputChain[c];
		prev = c - 1;
		next = c + 1;
		int curstrand = FragInput[cm].strand;
		int prevstrand, nextstrand;

		// qsb, qeb, tsb, teb form the unoverlapping region in the current Cluster;
		GenomePos qsb = FragInput[cm].qStart;
		GenomePos qeb = FragInput[cm].qEnd;
		GenomePos tsb = FragInput[cm].tStart;
		GenomePos teb = FragInput[cm].tEnd;

		if (prev != -1) {
			prevstrand = FragInput[inputChain[prev]].strand;
			GenomePos q1 = FragInput[inputChain[prev]].qStart;
			if (q1 > qsb and q1 < qeb) qeb = q1;

			if (inputChain.link[prev] == 0) { // prev is on the upper right corner of the current Cluster;
				GenomePos t1 = FragInput[inputChain[prev]].tStart;
				if (t1 > tsb and t1 < teb) teb = t1;				
			}
			else {
				GenomePos t2 = FragInput[inputChain[prev]].tEnd;
				if (t2 > tsb and t2 < teb) tsb = t2;
			}
		}

		if (next < inputChain.size()) {
 			nextstrand = FragInput[inputChain[next]].strand;
			GenomePos q2 = FragInput[inputChain[next]].qEnd; 
			if (q2 > qsb and q2 < qeb) qsb = q2;

			if (inputChain.link[cm] == 0) { // next is on the lower left corner of the current Cluster;
				GenomePos t2 = FragInput[inputChain[next]].tEnd;
				if (t2 > tsb and t2 < teb) tsb = t2;
			}
			else {
				GenomePos t1 = FragInput[inputChain[next]].tStart;
				if (t1 > tsb and t1 < teb) teb = t1;

			}
		}
		assert(qsb <= qeb);
		assert(tsb <= teb);
/*
		//
		// Check the previous Cluster to get the overlapping boundary;
		//
		if (prev != -1) {
			prevstrand = FragInput[inputChain[prev]].strand;
			GenomePos q1 = FragInput[inputChain[prev]].qStart;
			if (q1 > qsb and q1 < qeb) qeb = q1;
			GenomePos t1 = FragInput[inputChain[prev]].tStart;
			GenomePos t2 = FragInput[inputChain[prev]].tEnd;
			if (t1 > tsb and t1 < teb) { // prev is on the upper right corner of the current Cluster;
				teb = t1;
			}
			else if (t2 > tsb and t2 <= teb) { // prev is on the lower right corner of the current Cluster; //t1 <= tsb and t2 > tsb and t2 <= teb
				tsb = t2;
			}
		}

		//
		// Check the next Cluster to get the overlapping region;
		//
		if (next < inputChain.size()) {
 			nextstrand = FragInput[inputChain[next]].strand;
			GenomePos q2 = FragInput[inputChain[next]].qEnd; 
			if (q2 > qsb and q2 < qeb) qsb = q2;
			GenomePos t1 = FragInput[inputChain[next]].tStart;
			GenomePos t2 = FragInput[inputChain[next]].tEnd;
			if (t2 > teb and t1 > tsb and t1 < teb) { // next is on the upper left corner of the current Cluster;
				teb = t1;
			}
			else if (t1 < tsb and t2 > tsb and t2 <= teb) { // next is on the lower left corner of the current Cluster;
				tsb = t2;
			}			
		}
*/
		//
		// Insert points into H1;
		//
		for (unsigned int i = 0; i < FragInput[cm].matches.size(); i++) {

			// If the current anchor is inside the non-overlapping region, then insert two points for it;
			//
			// If it is forward stranded, then insert s1 and s2 for it;
			//
			if (curstrand == 0) { 

				// insert start point s1 into H1
				insertPointsPair(FragInput, H1, i + MatchStart[c], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, 
											FragInput[cm].matchesLengths[i], cm, c, 0, 1);
			/*
				Point s1;
				H1.push_back(s1);
				H1.back().ind = 1; // start
				H1.back().inv = 1; // forward direction
				H1.back().frag_num = i + MatchStart[c];
				H1.back().se.first = FragInput[cm].matches[i].first.pos;
				H1.back().se.second = FragInput[cm].matches[i].second.pos;										
				H1.back().orient = 1; // forward strand
				H1.back().clusterNum = cm; 
				H1.back().matchstartNum = c;

				// insert end point e1 into H1
				Point e1;
				H1.push_back(e1);
				H1.back().ind = 0; // end
				H1.back().inv = 1; // forward direction		
				H1.back().frag_num = i + MatchStart[c];
				H1.back().se.first = FragInput[cm].matches[i].first.pos + FragInput[cm].matchesLengths[i];
				H1.back().se.second = FragInput[cm].matches[i].second.pos + FragInput[cm].matchesLengths[i];
				H1.back().orient = 1; 	
				H1.back().clusterNum = cm; 
				H1.back().matchstartNum = c;
				*/
						
			}
			else { 
				// insert start point s2 into H1
				insertPointsPair(FragInput, H1, i + MatchStart[c], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, 
											FragInput[cm].matchesLengths[i], cm, c, 1, 0);
				/*
				Point s2;
				H1.push_back(s2);
				H1.back().ind = 1; // start
				H1.back().inv = 0; // backward direction
				H1.back().frag_num = i + MatchStart[c];
				H1.back().se.first = FragInput[cm].matches[i].first.pos; 
				H1.back().se.second = FragInput[cm].matches[i].second.pos + FragInput[cm].matchesLengths[i] - 1; //// TODO(Jingwen): fix this "-1" in the other SparseDP
				H1.back().orient = 0; // reversed strand			
				H1.back().clusterNum = cm; 
				H1.back().matchstartNum = c;	

				// insert end point e2 into H1
				Point e2;
				H1.push_back(e2);
				H1.back().ind = 0; // end
				H1.back().inv = 0; // backward direction		
				H1.back().frag_num = i + MatchStart[c];
				H1.back().se.first = FragInput[cm].matches[i].first.pos + FragInput[cm].matchesLengths[i];
				assert(FragInput[cm].matches[i].second.pos != 0); // IF this happens, just delete all the "-1" in the section of inserting points;
				H1.back().se.second = FragInput[cm].matches[i].second.pos - 1;
				H1.back().orient = 0; // reversed strand			
				H1.back().clusterNum = cm; 
				H1.back().matchstartNum = c;
				*/		
			}

		/*
			//
			// If the anchor is in the overlapping region, then insert the rest 2 points;
			if (!(FragInput[cm].matches[i].first.pos >= qsb and FragInput[cm].matches[i].first.pos + FragInput[cm].matchesLengths[i] <= qeb
					and FragInput[cm].matches[i].second.pos >= tsb and FragInput[cm].matches[i].second.pos + FragInput[cm].matchesLengths[i] <= teb)) {

				//
				// If it is forward stranded, then insert s2 and e2 for it;
				if (FragInput[cm].strand == 0) {
					// insert start point s2 into H1
					insertPointsPair(FragInput, H1, i + MatchStart[c], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, 
												FragInput[cm].matchesLengths[i], cm, c, 1, 1);
				}
				else {
					// 
					// If it is reversed stranded, then insert s1 and e1 for it;
					// insert start point s1 into H1
					insertPointsPair(FragInput, H1, i + MatchStart[c], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, 
												FragInput[cm].matchesLengths[i], cm, c, 0, 0);
				}
			}
		*/
			
			//
			// If 1. the anchor is in the left overlapping region, 2. curstrand != prevstrand, insert the rest 2 points;
			// If 1. the anchor is in the right overlapping region, 2. curstrand != nextstrand, insert the rest 2 points;
			//
			// The left overlapping region;
			//
			int add = 0;
			if (prev != -1 and curstrand != prevstrand) {
				if (curstrand == 0 and 
					(FragInput[cm].matches[i].first.pos + FragInput[cm].matchesLengths[i] > qeb or
					 FragInput[cm].matches[i].second.pos + FragInput[cm].matchesLengths[i] > teb)) {
					//
					// If it is forward stranded, then insert s2 and e2 for it;					
					// insert start point s2 into H1
					insertPointsPair(FragInput, H1, i + MatchStart[c], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, 
												FragInput[cm].matchesLengths[i], cm, c, 1, 1);
					add++;
				}
				else if (curstrand == 1 and 
					(FragInput[cm].matches[i].first.pos + FragInput[cm].matchesLengths[i] > qeb or 
					 FragInput[cm].matches[i].second.pos < tsb)) {
					// 
					// If it is reversed stranded, then insert s1 and e1 for it;
					// insert start point s1 into H1
					insertPointsPair(FragInput, H1, i + MatchStart[c], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, 
												FragInput[cm].matchesLengths[i], cm, c, 0, 0);
					add++;
				}
			}

			//
			// The right overlapping region;
			//
			if (next < inputChain.size() and add < 1 and curstrand != nextstrand) {

				if (curstrand == 0 and (FragInput[cm].matches[i].first.pos < qsb or FragInput[cm].matches[i].second.pos < tsb)) {
					//
					// If it is forward stranded, then insert s2 and e2 for it;					
					// insert start point s2 into H1
					insertPointsPair(FragInput, H1, i + MatchStart[c], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, 
												FragInput[cm].matchesLengths[i], cm, c, 1, 1);	
					add++;			
				}
				else if (curstrand == 1 and (FragInput[cm].matches[i].first.pos < qsb or 
							FragInput[cm].matches[i].second.pos + FragInput[cm].matchesLengths[i] > teb)) {
					// 
					// If it is reversed stranded, then insert s1 and e1 for it;
					// insert start point s1 into H1
					insertPointsPair(FragInput, H1, i + MatchStart[c], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, 
												FragInput[cm].matchesLengths[i], cm, c, 0, 0);
					add++;
				}				
			}
			assert(add <= 1);
		}
	}

	//clock_t begin = std::clock();

	//Sort the point by row
	//
	sort(H1.begin(), H1.end(), SortByRowOp<Point>()); // with same q and t coordinates, end point < start point
	std::vector<unsigned int> H2(H1.size());
	iota(H2.begin(), H2.end(), 0);
	sort(H2.begin(), H2.end(), SortByColOp<Point, unsigned int>(H1));

	
	vector<info> Row;
	vector<info> Col;
	GetRowInfo(H1, Row);
	GetColInfo(H1, H2, Col);
	//cerr << "Row: " << Row << "\n";
	//cerr << "Col: " << Col << "\n";

	unsigned int n1 = 0;
	unsigned int m1 = 0;
	unsigned int n2 = 0;
	unsigned int m2 = 0;

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

	//
	// Get SS_A_R1, SS_B_R1, SS_A_R2 and SS_B_R2 for each fragment 
	//
	vector<Fragment_Info> Value(totalMatch);
	for (unsigned int t = 0; t < Row.size(); ++t) {

		for (unsigned int tt = Row[t].pstart; tt < Row[t].pend; ++tt) {

			unsigned int ii = H1[tt].frag_num;
			int fi = H1[tt].clusterNum;
			int mi = H1[tt].matchstartNum;

			if (H1[tt].ind == 1 and H1[tt].inv == 1) { //H1[tt] is a start point (s1)
				Value[ii].SS_B_R1 = Row[t].SS_B1;
				Value[ii].counter_B_R1 = Row[t].SS_B1.size();
				Value[ii].val = FragInput[inputChain[mi]].matchesLengths[ii - MatchStart[mi]];
				Value[ii].clusterNum = fi;
				Value[ii].matchstartNum = mi;
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 0 and H1[tt].inv == 1) { // H1[tt] is an end point (e1)
				Value[ii].SS_A_R1 = Row[t].SS_A1;
				Value[ii].counter_A_R1 = Row[t].SS_A1.size();
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 1 and H1[tt].inv == 0) { //H1[tt] is a start point (s2)
				Value[ii].SS_B_R2 = Row[t].SS_B2;
				Value[ii].counter_B_R2 = Row[t].SS_B2.size();
				Value[ii].val = FragInput[inputChain[mi]].matchesLengths[ii - MatchStart[mi]];
				Value[ii].clusterNum = fi;
				Value[ii].matchstartNum = mi;
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[tt].orient;
			}
			else { // H1[tt] is an end point (e2)
				Value[ii].SS_A_R2 = Row[t].SS_A2;
				Value[ii].counter_A_R2 = Row[t].SS_A2.size();
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[tt].orient;				
			}
		}
	}

	//
	// Get SS_A_C1, SS_B_C1, SS_A_C2 and SS_B_C2 for each fragment
	//
	for (unsigned int t = 0; t < Col.size(); ++t) {
		for (unsigned int tt = Col[t].pstart; tt < Col[t].pend; ++tt) { 

			unsigned int ii = H1[H2[tt]].frag_num;
			int fi = H1[H2[tt]].clusterNum;
			int mi = H1[H2[tt]].matchstartNum;

			if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 1) { //H1[H2[tt]] a start point (s1)
				Value[ii].SS_B_C1 = Col[t].SS_B1;
				Value[ii].counter_B_C1 = Col[t].SS_B1.size();
				Value[ii].val = FragInput[inputChain[mi]].matchesLengths[ii - MatchStart[mi]];
				Value[ii].clusterNum = fi;
				Value[ii].matchstartNum = mi;
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 0 and H1[H2[tt]].inv == 1) { // H1[H2[tt]] is an end point (e1)
				Value[ii].SS_A_C1 = Col[t].SS_A1;
				Value[ii].counter_A_C1 = Col[t].SS_A1.size();
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 0) { //H1[H2[tt]] a start point (s2)
				Value[ii].SS_B_C2 = Col[t].SS_B2;
				Value[ii].counter_B_C2 = Col[t].SS_B2.size();
				Value[ii].val = FragInput[inputChain[mi]].matchesLengths[ii - MatchStart[mi]];
				Value[ii].clusterNum = fi;
				Value[ii].matchstartNum = mi;
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[H2[tt]].orient;				
			}
			else { // H1[H2[tt]] is an end point (e2)
				Value[ii].SS_A_C2 = Col[t].SS_A2;
				Value[ii].counter_A_C2 = Col[t].SS_A2.size();
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[H2[tt]].orient;			
			}

		}
	}

	//cerr << "Value: " << Value << endl;
	//cerr << "ProcessPoint\n";
	finalchain.InitializeOtherParts (MatchStart, totalMatch, Value);
	bool step = 1;
	ProcessPoint<Cluster>(H1, Row, SubR1, SubC1, SubR2, SubC2, Value, opts, LookUpTable, FragInput, inputChain, MatchStart, rate, step); 

	//
	// find the max_value for the FinalChain 
	//
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
	finalchain.chain.clear();
	TraceBack(SubR1, SubC1, SubR2, SubC2, Value, max_pos, finalchain.chain); // NOTICE: This chain is from the last anchors to the first anchor;

	// Clear SubR and SubC
	SubR1.Clear(eeR1);
	SubC1.Clear(eeC1);
	SubR2.Clear(eeR2);
	SubC2.Clear(eeC2);

	// get the time for the program
	//clock_t end = std::clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//cerr << "Time: " << elapsed_secs << endl;

	return 0;
}


//
// The input for this function is splitclusters
// fragments of different lengths
// This SDP needs to insert 4 points for any anchors
//
int SparseDP (const vector<Cluster> & FragInput, vector<Primary_chain> & Primary_chains, Options & opts, const vector<float> & LookUpTable, Read & read, float & rate) {

	std::vector<Point>  H1;
	// FragInput is vector<Cluster>
	// get points from FragInput and store them in H1		

	for (unsigned int i = 0; i < FragInput.size(); i++) {

			// insert start point s1 into H1
			Point s1;
			H1.push_back(s1);
			H1.back().ind = 1; // start
			H1.back().inv = 1; // forward direction
			H1.back().frag_num = i;
			H1.back().se.first = FragInput[i].qStart; 
			H1.back().se.second = FragInput[i].tStart;	
			if (FragInput[i].strand == 0) {
				H1.back().orient = 1; 
			}
			else {
				H1.back().orient = 0; 				
			}


			// insert end point e1 into H1
			Point e1;
			H1.push_back(e1);
			H1.back().ind = 0; // end
			H1.back().inv = 1; // forward direction		
			H1.back().frag_num = i;
			H1.back().se.first = FragInput[i].qEnd;
			H1.back().se.second = FragInput[i].tEnd;
			if (FragInput[i].strand == 0) {
				H1.back().orient = 1; 
			}
			else {
				H1.back().orient = 0; 				
			}

			// insert start point s2 into H1
			Point s2;
			H1.push_back(s2);
			H1.back().ind = 1; // start
			H1.back().inv = 0; // backward direction
			H1.back().frag_num = i;
			H1.back().se.first = FragInput[i].qStart; 
			H1.back().se.second = FragInput[i].tEnd;	
			if (FragInput[i].strand == 0) {
				H1.back().orient = 1; 
			}
			else {
				H1.back().orient = 0; 				
			}


			// insert end point e2 into H1
			Point e2;
			H1.push_back(e2);
			H1.back().ind = 0; // end
			H1.back().inv = 0; // backward direction		
			H1.back().frag_num = i;
			H1.back().se.first = FragInput[i].qEnd;
			H1.back().se.second = FragInput[i].tStart;	
			if (FragInput[i].strand == 0) {
				H1.back().orient = 1; 
			}
			else {
				H1.back().orient = 0; 				
			}

	}

	//clock_t begin = std::clock();

	//Sort the point by row
	sort(H1.begin(), H1.end(), SortByRowOp<Point>()); // with same q and t coordinates, end point < start point
	std::vector<unsigned int> H2(H1.size());
	iota(H2.begin(), H2.end(), 0);
	sort(H2.begin(), H2.end(), SortByColOp<Point, unsigned int>(H1));

	
	vector<info> Row;
	vector<info> Col;
	GetRowInfo(H1, Row);
	GetColInfo(H1, H2, Col);
	//cerr << "Row: " << Row << "\n";
	//cerr << "Col: " << Col << "\n";

	unsigned int n1 = 0;
	unsigned int m1 = 0;
	unsigned int n2 = 0;
	unsigned int m2 = 0;

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

	//
	// Get SS_A_R1, SS_B_R1, SS_A_R2 and SS_B_R2 for each fragment
	//
	std::vector<Fragment_Info> Value(FragInput.size());
	for (unsigned int t = 0; t < Row.size(); ++t) {

		for (unsigned int tt = Row[t].pstart; tt < Row[t].pend; ++tt) {

			unsigned int ii = H1[tt].frag_num;

			if (H1[tt].ind == 1 and H1[tt].inv == 1) { //H1[tt] is a start point (s1)
				Value[ii].SS_B_R1 = Row[t].SS_B1;
				Value[ii].counter_B_R1 = Row[t].SS_B1.size();
				Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 0 and H1[tt].inv == 1) { // H1[tt] is an end point (e1)
				Value[ii].SS_A_R1 = Row[t].SS_A1;
				Value[ii].counter_A_R1 = Row[t].SS_A1.size();
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 1 and H1[tt].inv == 0) { //H1[tt] is a start point (s2)
				Value[ii].SS_B_R2 = Row[t].SS_B2;
				Value[ii].counter_B_R2 = Row[t].SS_B2.size();
				Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[tt].orient;
			}
			else { // H1[tt] is an end point (e2)
				Value[ii].SS_A_R2 = Row[t].SS_A2;
				Value[ii].counter_A_R2 = Row[t].SS_A2.size();
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[tt].orient;				
			}
		}
	}

	//
	// Get SS_A_C1, SS_B_C1, SS_A_C2 and SS_B_C2 for each fragment
	//
	for (unsigned int t = 0; t < Col.size(); ++t) {
		for (unsigned int tt = Col[t].pstart; tt < Col[t].pend; ++tt) { 

			unsigned int ii = H1[H2[tt]].frag_num;

			if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 1) { //H1[H2[tt]] a start point (s1)
				Value[ii].SS_B_C1 = Col[t].SS_B1;
				Value[ii].counter_B_C1 = Col[t].SS_B1.size();
				Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 0 and H1[H2[tt]].inv == 1) { // H1[H2[tt]] is an end point (e1)
				Value[ii].SS_A_C1 = Col[t].SS_A1;
				Value[ii].counter_A_C1 = Col[t].SS_A1.size();
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 0) { //H1[H2[tt]] a start point (s2)
				Value[ii].SS_B_C2 = Col[t].SS_B2;
				Value[ii].counter_B_C2 = Col[t].SS_B2.size();
				Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[H2[tt]].orient;				
			}
			else { // H1[H2[tt]] is an end point (e2)
				Value[ii].SS_A_C2 = Col[t].SS_A2;
				Value[ii].counter_A_C2 = Col[t].SS_A2.size();
				//Value[ii].val = std::min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[H2[tt]].orient;			
			}

		}
	}


	//cerr << "Value: " << Value << endl;
	//cerr << "ProcessPoint\n";
	bool step = 0;
	ProcessPoint<Cluster>(H1, Row, SubR1, SubC1, SubR2, SubC2, Value, opts, LookUpTable, FragInput, rate, step);	

	//
	// Trace back chains; There are two parameters: PrimaryAln, SecondaryAln
	//
	DecidePrimaryChains(FragInput, SubR1, SubC1, SubR2, SubC2, Value, Primary_chains, read, opts);

	// Clear SubR and SubC
	SubR1.Clear(eeR1);
	SubC1.Clear(eeC1);
	SubR2.Clear(eeR2);
	SubC2.Clear(eeC2);

	// get the time for the program
	//clock_t end = std::clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//cerr << "Time: " << elapsed_secs << endl;

	return 0;
}








#endif