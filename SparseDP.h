// This program is implemented Sparse Dynamic Programming algorithm described in David Eppstein paper
// Author: Jingwen Ren

#ifndef SPARSE_DP_
#define SPARSE_DP_

#include <iostream>
#include <string>
#include <utility>
#include <algorithm> // lower_bound
#include <numeric> //floor
#include <cmath>
#include <set>
#include <iterator>
#include <assert.h>
#include <chrono> // generate random number
#include <random> // generate random number
#include <type_traits>
 

#include "SubProblem.h"
#include "Sorting.h"
#include "SubRountine.h"
#include "Fragment_Info.h"
#include "Info.h"
// #include "overload.h"
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

void clear(vector<Point> &H1, vector<unsigned int> &H2, StackOfSubProblems &SubR1, StackOfSubProblems &SubC1, StackOfSubProblems &SubR2, 
		StackOfSubProblems &SubC2, vector<Fragment_Info> &Value, vector<info> &Row, vector<info> &Col) {
	SubR1.clear();
	SubC1.clear();
	SubR2.clear();
	SubC2.clear();
	Value.clear();
	Row.clear();
	Col.clear();
	H1.clear();
	H2.clear();
}

template <typename Tup>
void ComputeMatchStart(vector<int> & MatchNum, int & totalMatch, const vector<Tup *> & FragInput, SplitChain & inputChain){
	totalMatch = FragInput[inputChain[0]]->size();
	for (int c = 1; c < inputChain.size(); c++) {
		totalMatch += FragInput[inputChain[c]]->size();
		MatchNum[c] = FragInput[inputChain[c - 1]]->size();
	}
	for (int c = 1; c < MatchNum.size(); c++) {
		MatchNum[c] = MatchNum[c-1] + MatchNum[c]; 
	}
}

template <typename Tup>
void ComputeMatchStart(vector<int> & MatchNum, int & totalMatch, const vector<Tup> & FragInput, SplitChain & inputChain){
	totalMatch = FragInput[inputChain[0]].size();
	for (int c = 1; c < inputChain.size(); c++) {
		totalMatch += FragInput[inputChain[c]].size();
		MatchNum[c] = FragInput[inputChain[c - 1]].size();
	}
	for (int c = 1; c < MatchNum.size(); c++) {
		MatchNum[c] = MatchNum[c-1] + MatchNum[c]; 
	}
}
void
insertPointsPair(vector<Point> & H1, unsigned int frag_num, GenomePos qs, GenomePos ts, int length, int &clusterNum, int Pair, int strand) {
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
		// H1.back().matchstartNum = matchstartNum;

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
		// H1.back().matchstartNum = matchstartNum;		
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
		H1.back().se.second = ts + length; //// TODO(Jingwen): fix this "-1" in the other SparseDP
		H1.back().orient = strand; // reversed strand			
		H1.back().clusterNum = clusterNum; 
		// H1.back().matchstartNum = matchstartNum;	

		// insert end point e2 into H1
		Point e2;
		H1.push_back(e2);
		H1.back().ind = 0; // end
		H1.back().inv = 0; // backward direction		
		H1.back().frag_num = frag_num;
		H1.back().se.first = qs + length;
		// assert(ts != 0); // IF this happens, just delete all the "-1" in the section of inserting points;
		H1.back().se.second = ts;
		H1.back().orient = strand; // reversed strand			
		H1.back().clusterNum = clusterNum; 
		// H1.back().matchstartNum = matchstartNum;		
	}

}


void 
PassValueToD1 (unsigned int i, vector<Fragment_Info> & Value, const vector<Point> & H1, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, long int & ForwardDiag) {

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

			vector<unsigned int>::iterator it = Lower_Bound<vector<unsigned int>::iterator, long int>(SubR1[u].D.begin(), SubR1[u].D.end(), ForwardDiag, SubR1[u].Di);

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

			vector<unsigned int>::reverse_iterator it = Lower_Bound<vector<unsigned int>::reverse_iterator, long int>(SubC1[u].D.rbegin(), SubC1[u].D.rend(), ForwardDiag, SubC1[u].Di);
		
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
PassValueToD2 (unsigned int i, vector<Fragment_Info> & Value, const vector<Point> & H1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2, long int & BackDiag) {

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

			vector<unsigned int>::reverse_iterator it = Lower_Bound<vector<unsigned int>::reverse_iterator, long int>(SubR2[u].D.rbegin(), SubR2[u].D.rend(), BackDiag, SubR2[u].Di);
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

			vector<unsigned int>::iterator it = Lower_Bound<vector<unsigned int>::iterator, long int>(SubC2[u].D.begin(), SubC2[u].D.end(), BackDiag, SubC2[u].Di);
		
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
		
//
// Each fragment has different length, so we need to pass FragInput to this function
// This function is for splitCluster
//
template<typename Tup>
void ProcessPoint (const vector<Point> & H1, vector<info> & V, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1,
				 StackOfSubProblems & SubR2, StackOfSubProblems & SubC2, vector<Fragment_Info> & Value, Options & opts, 
				 const vector<float> & LookUpTable, const vector<Tup> & FragInput, float & rate) { // vector<ClusterCoordinates> & FragInput

	bool step_sdp = 0;
	for (unsigned int i = 0; i < H1.size(); ++i) { // process points by row

		long int ForwardDiag = (long int)(H1[i].se.second) - (long int)(H1[i].se.first);
		long int BackDiag = (long int)(H1[i].se.second) + (long int)(H1[i].se.first);

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
					vector<unsigned int>::iterator t = Lower_Bound<vector<unsigned int>::iterator,long int>(SubR1[j].E.begin(), SubR1[j].E.end(), ForwardDiag, SubR1[j].Ei);
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
						Maximization (SubR1[j].now, SubR1[j].last, SubR1[j].Di, SubR1[j].Ei, SubR1[j].Dv, SubR1[j].Db, SubR1[j].Block, SubR1[j].S_1, LookUpTable, opts, step_sdp); // TODO(Jingwen) anything change for SubC????
						SubR1[j].last = SubR1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubR1[j].S_1, SubR1[j].Ei, SubR1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
						

						SubR1[j].Ev[i1] = SubR1[j].Dv[i2] + w(SubR1[j].Di[i2], SubR1[j].Ei[i1], LookUpTable, opts, step_sdp) + FragInput[ii].Val*rate;
																//min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
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
					vector<unsigned int>::reverse_iterator t = Lower_Bound<vector<unsigned int>::reverse_iterator,long int>(SubC1[j].E.rbegin(), SubC1[j].E.rend(), ForwardDiag, SubC1[j].Ei);
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
						Maximization (SubC1[j].now, SubC1[j].last, SubC1[j].Di, SubC1[j].Ei, SubC1[j].Dv, SubC1[j].Db, SubC1[j].Block, SubC1[j].S_1, LookUpTable, opts, step_sdp); 
						SubC1[j].last = SubC1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubC1[j].S_1, SubC1[j].Ei, SubC1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					
	//					SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts) + opts.globalK; 
	//					SubC1[j].Ep[i1] = i2;							

						SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts, step_sdp) + FragInput[ii].Val*rate;
																//min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
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
					vector<unsigned int>::reverse_iterator t = Lower_Bound<vector<unsigned int>::reverse_iterator,long int>(SubR2[j].E.rbegin(), SubR2[j].E.rend(), BackDiag, SubR2[j].Ei);
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
						Maximization (SubR2[j].now, SubR2[j].last, SubR2[j].Di, SubR2[j].Ei, SubR2[j].Dv, SubR2[j].Db, SubR2[j].Block, SubR2[j].S_1, LookUpTable, opts, step_sdp); // TODO(Jingwen) anything change for SubC????
						SubR2[j].last = SubR2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that BackDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for BackDiag
						FindValueInBlock(BackDiag, SubR2[j].S_1, SubR2[j].Ei, SubR2[j].Block, i1, i2);
						//cerr << "the index in Ei that BackDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for BackDiag ---- i2: " << i2 << "\n";

						SubR2[j].Ev[i1] = SubR2[j].Dv[i2] + w(SubR2[j].Di[i2], SubR2[j].Ei[i1], LookUpTable, opts, step_sdp) + FragInput[ii].Val*rate;
																	//min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
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
					vector<unsigned int>::iterator t = Lower_Bound<vector<unsigned int>::iterator,long int>(SubC2[j].E.begin(), SubC2[j].E.end(), BackDiag, SubC2[j].Ei);
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
						Maximization (SubC2[j].now, SubC2[j].last, SubC2[j].Di, SubC2[j].Ei, SubC2[j].Dv, SubC2[j].Db, SubC2[j].Block, SubC2[j].S_1, LookUpTable, opts, step_sdp); 
						SubC2[j].last = SubC2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores 
						unsigned int i2;
						FindValueInBlock(BackDiag, SubC2[j].S_1, SubC2[j].Ei, SubC2[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					

						SubC2[j].Ev[i1] = SubC2[j].Dv[i2] + w(SubC2[j].Di[i2], SubC2[j].Ei[i1], LookUpTable, opts, step_sdp) + FragInput[ii].Val*rate;
																	//min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate; 
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
void ProcessPoint (const vector<Point> & H1, vector<info> & V, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1,
				 StackOfSubProblems & SubR2, StackOfSubProblems & SubC2, vector<Fragment_Info> & Value, Options & opts, 
				 const vector<float> & LookUpTable, const vector<Tup*> &FragInput, float & rate) { // vector<ClusterCoordinates> & FragInput
	bool step_sdp = 1;
	for (unsigned int i = 0; i < H1.size(); ++i) { // process points by row

		long int ForwardDiag = (long int)(H1[i].se.second) - (long int)(H1[i].se.first);
		long int BackDiag = (long int)(H1[i].se.second) + (long int)(H1[i].se.first);

		//cerr << "\n\n\n\nprocessing point " << H1[H3[i]] << endl;
		//cerr << "the Value[" << H1[H3[i]].frag_num << "] of this point:  " <<  Value[H1[H3[i]].frag_num] << "\n";
		
		unsigned int ii = H1[i].frag_num;
		// int mi = H1[i].matchstartNum;
		int fi = H1[i].clusterNum;
		int ms = FragInput[fi]->matchStart;

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
					vector<unsigned int>::iterator t = Lower_Bound<vector<unsigned int>::iterator,long int>(SubR1[j].E.begin(), SubR1[j].E.end(), ForwardDiag, SubR1[j].Ei);
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
						Maximization (SubR1[j].now, SubR1[j].last, SubR1[j].Di, SubR1[j].Ei, SubR1[j].Dv, SubR1[j].Db, SubR1[j].Block, SubR1[j].S_1, LookUpTable, opts, step_sdp); // TODO(Jingwen) anything change for SubC????
						SubR1[j].last = SubR1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubR1[j].S_1, SubR1[j].Ei, SubR1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
						

						SubR1[j].Ev[i1] = SubR1[j].Dv[i2] + w(SubR1[j].Di[i2], SubR1[j].Ei[i1], LookUpTable, opts, step_sdp) + 
																rate * FragInput[fi]->length(ii - ms);
																//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]];
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
					vector<unsigned int>::reverse_iterator t = Lower_Bound<vector<unsigned int>::reverse_iterator,long int>(SubC1[j].E.rbegin(), SubC1[j].E.rend(), ForwardDiag, SubC1[j].Ei);
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
						Maximization (SubC1[j].now, SubC1[j].last, SubC1[j].Di, SubC1[j].Ei, SubC1[j].Dv, SubC1[j].Db, SubC1[j].Block, SubC1[j].S_1, LookUpTable, opts, step_sdp); 
						SubC1[j].last = SubC1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubC1[j].S_1, SubC1[j].Ei, SubC1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					
	//					SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts) + opts.globalK; 
	//					SubC1[j].Ep[i1] = i2;							

						SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts, step_sdp) + 
																	rate * FragInput[fi]->length(ii - ms);
																//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]];
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
					vector<unsigned int>::reverse_iterator t = Lower_Bound<vector<unsigned int>::reverse_iterator,long int>(SubR2[j].E.rbegin(), SubR2[j].E.rend(), BackDiag, SubR2[j].Ei);
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
						Maximization (SubR2[j].now, SubR2[j].last, SubR2[j].Di, SubR2[j].Ei, SubR2[j].Dv, SubR2[j].Db, SubR2[j].Block, SubR2[j].S_1, LookUpTable, opts, step_sdp); // TODO(Jingwen) anything change for SubC????
						SubR2[j].last = SubR2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that BackDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for BackDiag
						FindValueInBlock(BackDiag, SubR2[j].S_1, SubR2[j].Ei, SubR2[j].Block, i1, i2);
						//cerr << "the index in Ei that BackDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for BackDiag ---- i2: " << i2 << "\n";

						SubR2[j].Ev[i1] = SubR2[j].Dv[i2] + w(SubR2[j].Di[i2], SubR2[j].Ei[i1], LookUpTable, opts, step_sdp) + 
														rate * FragInput[fi]->length(ii - ms);
																	//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]];
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
					vector<unsigned int>::iterator t = Lower_Bound<vector<unsigned int>::iterator,long int>(SubC2[j].E.begin(), SubC2[j].E.end(), BackDiag, SubC2[j].Ei);
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
						Maximization (SubC2[j].now, SubC2[j].last, SubC2[j].Di, SubC2[j].Ei, SubC2[j].Dv, SubC2[j].Db, SubC2[j].Block, SubC2[j].S_1, LookUpTable, opts, step_sdp); 
						SubC2[j].last = SubC2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores 
						unsigned int i2;
						FindValueInBlock(BackDiag, SubC2[j].S_1, SubC2[j].Ei, SubC2[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
					

						SubC2[j].Ev[i1] = SubC2[j].Dv[i2] + w(SubC2[j].Di[i2], SubC2[j].Ei[i1], LookUpTable, opts, step_sdp) + 
														rate * FragInput[fi]->length(ii - ms);
																	// FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]];
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
// Each fragment has different length,
// This function is for pure matches;
//
template<typename Tup>
void ProcessPoint (const vector<Point> & H1, vector<info> & V, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1,
				 StackOfSubProblems & SubR2, StackOfSubProblems & SubC2, vector<Fragment_Info> & Value, Options & opts, 
				 const vector<float> & LookUpTable, const vector<Tup> &FragInput, const vector<int> & MatchStart, float & rate) { 
	bool step_sdp = 0;
	for (unsigned int i = 0; i < H1.size(); ++i) { // process points by row

		long int ForwardDiag = (long int)(H1[i].se.second) - (long int)(H1[i].se.first);
		long int BackDiag = (long int)(H1[i].se.second) + (long int)(H1[i].se.first);

		//cerr << "\n\n\n\nprocessing point " << H1[H3[i]] << endl;
		//cerr << "the Value[" << H1[H3[i]].frag_num << "] of this point:  " <<  Value[H1[H3[i]].frag_num] << "\n";
		
		unsigned int ii = H1[i].frag_num;
		// int mi = H1[i].matchstartNum;
		int fi = H1[i].clusterNum;
		int ms = MatchStart[fi];

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
					vector<unsigned int>::iterator t = Lower_Bound<vector<unsigned int>::iterator,long int>(SubR1[j].E.begin(), SubR1[j].E.end(), ForwardDiag, SubR1[j].Ei);
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
						Maximization (SubR1[j].now, SubR1[j].last, SubR1[j].Di, SubR1[j].Ei, SubR1[j].Dv, SubR1[j].Db, SubR1[j].Block, SubR1[j].S_1, LookUpTable, opts, step_sdp); // TODO(Jingwen) anything change for SubC????
						SubR1[j].last = SubR1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubR1[j].S_1, SubR1[j].Ei, SubR1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
						

						SubR1[j].Ev[i1] = SubR1[j].Dv[i2] + w(SubR1[j].Di[i2], SubR1[j].Ei[i1], LookUpTable, opts, step_sdp) + rate * FragInput[fi].matchesLengths[ii - ms];
						SubR1[j].Ep[i1] = i2;							

						//cerr << "SubR1[" << j << "].Ev[" << i1 << "]: " << SubR1[j].Ev[i1] << ", SubR1[" << j << "].Ep[" << i1 << "]: " << SubR1[j].Ep[i1] << "\n"; 


						// Update the value of this point
						//TODO(Jingwen): if this point is a s1 of a reverse orientated anchor, only update the value when  SubR1[j].Dp[Ep[i1]] points to a forward orientated anchor
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
					vector<unsigned int>::reverse_iterator t = Lower_Bound<vector<unsigned int>::reverse_iterator,long int>(SubC1[j].E.rbegin(), SubC1[j].E.rend(), ForwardDiag, SubC1[j].Ei);
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
						Maximization (SubC1[j].now, SubC1[j].last, SubC1[j].Di, SubC1[j].Ei, SubC1[j].Dv, SubC1[j].Db, SubC1[j].Block, SubC1[j].S_1, LookUpTable, opts, step_sdp); 
						SubC1[j].last = SubC1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubC1[j].S_1, SubC1[j].Ei, SubC1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
						SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts, step_sdp) + rate * FragInput[fi].matchesLengths[ii - ms];
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
					vector<unsigned int>::reverse_iterator t = Lower_Bound<vector<unsigned int>::reverse_iterator,long int>(SubR2[j].E.rbegin(), SubR2[j].E.rend(), BackDiag, SubR2[j].Ei);
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
						Maximization (SubR2[j].now, SubR2[j].last, SubR2[j].Di, SubR2[j].Ei, SubR2[j].Dv, SubR2[j].Db, SubR2[j].Block, SubR2[j].S_1, LookUpTable, opts, step_sdp); // TODO(Jingwen) anything change for SubC????
						SubR2[j].last = SubR2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that BackDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for BackDiag
						FindValueInBlock(BackDiag, SubR2[j].S_1, SubR2[j].Ei, SubR2[j].Block, i1, i2);
						//cerr << "the index in Ei that BackDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for BackDiag ---- i2: " << i2 << "\n";

						SubR2[j].Ev[i1] = SubR2[j].Dv[i2] + w(SubR2[j].Di[i2], SubR2[j].Ei[i1], LookUpTable, opts, step_sdp) + rate * FragInput[fi].matchesLengths[ii - ms];
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
					vector<unsigned int>::iterator t = Lower_Bound<vector<unsigned int>::iterator,long int>(SubC2[j].E.begin(), SubC2[j].E.end(), BackDiag, SubC2[j].Ei);
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
						Maximization (SubC2[j].now, SubC2[j].last, SubC2[j].Di, SubC2[j].Ei, SubC2[j].Dv, SubC2[j].Db, SubC2[j].Block, SubC2[j].S_1, LookUpTable, opts, step_sdp); 
						SubC2[j].last = SubC2[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
			
						unsigned int i1 = *t; // i1 stores 
						unsigned int i2;
						FindValueInBlock(BackDiag, SubC2[j].S_1, SubC2[j].Ei, SubC2[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
						SubC2[j].Ev[i1] = SubC2[j].Dv[i2] + w(SubC2[j].Di[i2], SubC2[j].Ei[i1], LookUpTable, opts, step_sdp) + rate * FragInput[fi].matchesLengths[ii - ms];
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
// No recursive; 
void 
TraceBack (StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
		   const vector<Fragment_Info> & Value, unsigned int & i, vector<unsigned int> & onechain, vector<bool> & link, 
		   vector<bool> & used) {

	long int prev_sub = Value[i].prev_sub;
	long int prev_ind = Value[i].prev_ind;
	if (used[i] == 0) {
		onechain.push_back(i);
		used[i] = 1;

		while (prev_sub != -1 and prev_ind != -1) {
			if (Value[i].prev == 1 and Value[i].inv == 1) { // The previous subproblem is SubR1
				unsigned int ind = SubR1[prev_sub].Ep[prev_ind];
				if (used[SubR1[prev_sub].Dp[ind]] == 0) {
					link.push_back(0);
					i = SubR1[prev_sub].Dp[ind];
				}
				else {
					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
						used[onechain[lu]] = 0;
					}
					onechain.clear();
					link.clear();
					break;
				}
			}	
			else if (Value[i].prev == 1 and Value[i].inv == 0) { // The previous subproblem is SubR2
				unsigned int ind = SubR2[prev_sub].Ep[prev_ind];
				if (used[SubR2[prev_sub].Dp[ind]] == 0) {
					link.push_back(1);
					i = SubR2[prev_sub].Dp[ind];
				}
				else {
					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
						used[onechain[lu]] = 0;
					}
					onechain.clear();	
					link.clear();	
					break;			
				}
			}
			else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
				unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
				if (used[SubC1[prev_sub].Dp[ind]] == 0) {
					link.push_back(0);
					i = SubC1[prev_sub].Dp[ind];
				}
				else {
					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
						used[onechain[lu]] = 0;
					}
					onechain.clear();	
					link.clear();	
					break;				
				}
			}		
			else { // The previous subproblem is SubC2
				unsigned int ind = SubC2[prev_sub].Ep[prev_ind];
				if (used[SubC2[prev_sub].Dp[ind]] == 0) {
					link.push_back(1);
					i = SubC2[prev_sub].Dp[ind];
				}
				else {
					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
						used[onechain[lu]] = 0;
					}
					onechain.clear();	
					link.clear();	
					break;			
				}
			}		
			prev_sub = Value[i].prev_sub;
			prev_ind = Value[i].prev_ind;
			if (used[i] == 0) {
				onechain.push_back(i);
				used[i] = 1;
			}
			else {
				for (unsigned int lu = 0; lu < onechain.size(); lu++) {
					used[onechain[lu]] = 0;
				}
				onechain.clear();	
				link.clear();	
				break;		
			}	
		}	
	}	
}

// //
// // This function is for tracing back such a chain that every anchor on this chain is unused;
// // recursive form 
// void 
// TraceBack (StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
// 		   const vector<Fragment_Info> & Value, unsigned int & i, vector<unsigned int> & onechain, vector<bool> & link, 
// 		   vector<bool> & used) {

// 	long int prev_sub = Value[i].prev_sub;
// 	long int prev_ind = Value[i].prev_ind;
// 	if (used[i] == 0) {
// 		onechain.push_back(i);
// 		used[i] = 1;
	
// 		if (prev_sub != -1 and prev_ind != -1) {

// 			if (Value[i].prev == 1 and Value[i].inv == 1) { // The previous subproblem is SubR1
// 				unsigned int ind = SubR1[prev_sub].Ep[prev_ind];
// 				if (used[SubR1[prev_sub].Dp[ind]] == 0) {
// 					link.push_back(0);
// 					TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR1[prev_sub].Dp[ind], onechain, link, used);					
// 				}
// 				else {
// 					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
// 						used[onechain[lu]] = 0;
// 					}
// 					onechain.clear();
// 					link.clear();
// 				}
// 			}
// 			else if (Value[i].prev == 1 and Value[i].inv == 0) { // The previous subproblem is SubR2
// 				unsigned int ind = SubR2[prev_sub].Ep[prev_ind];
// 				if (used[SubR2[prev_sub].Dp[ind]] == 0) {
// 					link.push_back(1);
// 					TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR2[prev_sub].Dp[ind], onechain, link, used);
// 				}
// 				else {
// 					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
// 						used[onechain[lu]] = 0;
// 					}
// 					onechain.clear();	
// 					link.clear();				
// 				}
// 			}
// 			else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
// 				unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
// 				if (used[SubC1[prev_sub].Dp[ind]] == 0) {
// 					link.push_back(0);
// 					TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC1[prev_sub].Dp[ind], onechain, link, used);
// 				}
// 				else {
// 					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
// 						used[onechain[lu]] = 0;
// 					}
// 					onechain.clear();	
// 					link.clear();					
// 				}
// 			}		
// 			else { // The previous subproblem is SubC2
// 				unsigned int ind = SubC2[prev_sub].Ep[prev_ind];
// 				if (used[SubC2[prev_sub].Dp[ind]] == 0) {
// 					link.push_back(1);
// 					TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC2[prev_sub].Dp[ind], onechain, link, used);
// 				}
// 				else {
// 					for (unsigned int lu = 0; lu < onechain.size(); lu++) {
// 						used[onechain[lu]] = 0;
// 					}
// 					onechain.clear();	
// 					link.clear();				
// 				}
// 			}

// 		}
// 	}
// }

//
// This function is for tracing back a chain;
//
void 
TraceBack (StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
					const vector<Fragment_Info> & Value, unsigned int i, vector<unsigned int> & Chain, vector<bool> & link) {

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
			link.push_back(0);
			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR1[prev_sub].Dp[ind], Chain);
		}
		else if (Value[i].prev == 1 and Value[i].inv == 0) { // The previous subproblem is SubR2
			assert(prev_sub <  SubR2.StackSub.size());
			assert(prev_ind < SubR2[prev_sub].Ep.size());
			unsigned int ind = SubR2[prev_sub].Ep[prev_ind];
			assert(ind < SubR2[prev_sub].Dp.size());
			i = SubR2[prev_sub].Dp[ind];
			link.push_back(1);
			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR2[prev_sub].Dp[ind], Chain);
		}
		else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
			assert(prev_sub <  SubC1.StackSub.size());
			assert(prev_ind < SubC1[prev_sub].Ep.size());
			unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
			assert(ind < SubC1[prev_sub].Dp.size());
			i = SubC1[prev_sub].Dp[ind];
			link.push_back(0);
			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC1[prev_sub].Dp[ind], Chain);
		}		
		else { // The previous subproblem is SubC2
			assert(prev_sub < SubC2.StackSub.size());
			assert(prev_ind < SubC2[prev_sub].Ep.size());
			unsigned int ind = SubC2[prev_sub].Ep[prev_ind];
			assert(ind < SubC2[prev_sub].Dp.size());
			i = SubC2[prev_sub].Dp[ind];
			link.push_back(1);
			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC2[prev_sub].Dp[ind], Chain);
		}
		prev_sub = Value[i].prev_sub;
		prev_ind = Value[i].prev_ind;		
		Chain.push_back(i);
	}
}

//
// Compute the # of anchors on this chain
//
void 
ComputeNumOfAnchors (const vector<Cluster> & FragInput, const vector<unsigned int> & onechain, int &Num_Anchors) {
	for (int oc = 0; oc < onechain.size(); oc++) {
		Num_Anchors += FragInput[onechain[oc]].NumofAnchors0;
	}
}

//
// This function is for splitClusters;
//
void 
DecidePrimaryChains(const vector<Cluster> & FragInput, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
					const vector<Fragment_Info> & Value, vector<Primary_chain> &Primary_chains, Read & read, Options & opts) {

	vector<bool> used(Value.size(), 0);
	Fragment_valueOrder fragments_valueOrder(&Value);
	float value_thres = max(opts.alnthres * fragments_valueOrder[0], fragments_valueOrder[0] - 100*opts.globalK);//30 for 50kb
	//float value_thres = opts.alnthres*fragments_valueOrder[0];
	// cerr << "value_thres: " << value_thres << endl;
	// cerr << "fragments_valueOrder[0]: " << fragments_valueOrder[0] << " fragments_valueOrder[1]: " << 
	// 			fragments_valueOrder[1] << endl;
	int fv = 0;
	while (fv < fragments_valueOrder.size() and fragments_valueOrder[fv] >= value_thres) {
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
			// If this chain overlap with read greater than 0.5%, insert it to chains
			//
			assert(qEnd - qStart>10);
			if (((float)(qEnd - qStart)/read.length) > 0.005) {
				//
				// Compute the # of anchors on this chain
				//
				int Num_Anchors = 0;
				ComputeNumOfAnchors (FragInput, onechain, Num_Anchors);
				//
				// Compare onechain to all the primary chains we've found. 
				// If onechain overlaps with one primary chain over 50% ---> onechain is a secondary chain 
				// If onechain overlaps with all the primary chains less than 50% ---> onechain is another primary chain
				//
				if (Primary_chains.size() == 0) {
					Primary_chain Pc(CHain(qStart, qEnd, tStart, tEnd, onechain, link, fragments_valueOrder[fv], Num_Anchors));
					Primary_chains.push_back(Pc);
				} 
				else if (Primary_chains[0].chains.size() < opts.NumAln) {
					Primary_chains[0].chains.push_back(CHain(qStart, qEnd, tStart, tEnd, onechain, link, fragments_valueOrder[fv], Num_Anchors));					
				}
				else {
					break;
					// bool newpr = 1, inserted = 0;
					// int p = 0;
					// while (p < Primary_chains.size()) {
					// 	if (Primary_chains[p].chains[0].OverlapsOnQ(qStart, qEnd, 0.5)) {
					// 		if (!Primary_chains[p].chains[0].OverlapsOnT(tStart, tEnd, 0.3)) {
					// 			if (Primary_chains[p].chains.size() < opts.NumAln) {
					// 				Primary_chains[p].chains.push_back(CHain(qStart, qEnd, tStart, tEnd, onechain, 
					// 															link, fragments_valueOrder[fv], Num_Anchors));
					// 				inserted = 1;								
					// 			}
					// 			break;
					// 		}
					// 	}
					// 	else{
					// 		newpr = 0;
					// 	}
					// 	++p;
					// }			
					// if (p == Primary_chains.size() - 1 and inserted == 0 and newpr == 0) {		
					// 	if (Primary_chains.size() < 2) { // TODO(Jingwen): how to decide the number of Primary alignments
					// 		Primary_chain Pc(CHain(qStart, qEnd, tStart, tEnd, onechain, link, fragments_valueOrder[fv], Num_Anchors));
					// 		Primary_chains.push_back(Pc);
					// 	}	
					// 	else break;		
					// }	
				}
			}
			else break;				
		}
		onechain.clear();
		link.clear();
		fv++;
	}
}

//
// This function is for pure matches
//
void 
DecidePrimaryChains(vector<Cluster> & FragInput, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, StackOfSubProblems & SubR2, StackOfSubProblems & SubC2,
					const vector<Fragment_Info> & Value, vector<UltimateChain> &chains, Read & read, Options & opts, vector<int> &MatchStart) {

	Fragment_valueOrder fragments_valueOrder(&Value);
	float value_thres = max(0.95f * fragments_valueOrder[0], fragments_valueOrder[0] - 100*opts.globalK);//30 for 50kb
	//float value_thres = opts.alnthres*fragments_valueOrder[0];
	// cerr << "value_thres: " << value_thres << endl;
	// cerr << "fragments_valueOrder[0]: " << fragments_valueOrder[0] << " fragments_valueOrder[1]: " << 
	// 			fragments_valueOrder[1] << endl;
	int fv = 0;
	while (fv < fragments_valueOrder.size() and fragments_valueOrder[fv] >= value_thres) {
		unsigned int d = fragments_valueOrder.index[fv];
		if (fv < opts.NumAln) chains.push_back(UltimateChain(&FragInput));
		else break;
		TraceBack(SubR1, SubC1, SubR2, SubC2, Value, d, chains.back().chain, chains.back().link);

		if (chains.back().chain.size() != 0) {
			// Note: onechain store index from the last one to the first one
			int f = chains.back().chain[0]; int l = chains.back().chain.back();
			int fi = Value[f].clusterNum; int li = Value[l].clusterNum;
			chains.back().QEnd = FragInput[fi].matches[f - MatchStart[fi]].first.pos + FragInput[fi].matchesLengths[f - MatchStart[fi]];
			chains.back().QStart = FragInput[li].matches[l - MatchStart[li]].first.pos;
			chains.back().TEnd = FragInput[fi].matches[f - MatchStart[fi]].second.pos + FragInput[fi].matchesLengths[f - MatchStart[fi]];
			chains.back().TStart = FragInput[li].matches[l - MatchStart[li]].second.pos;

			for (int c = 0; c < chains.back().chain.size(); c++) {
				f = chains.back().chain[c]; fi = Value[f].clusterNum;
				chains.back().QEnd = max(chains.back().QEnd , FragInput[fi].matches[f - MatchStart[fi]].first.pos + FragInput[fi].matchesLengths[f - MatchStart[fi]]);
				chains.back().TEnd = max(chains.back().TEnd, FragInput[fi].matches[f - MatchStart[fi]].second.pos + FragInput[fi].matchesLengths[f - MatchStart[fi]]);
				chains.back().QStart = min(chains.back().QStart, FragInput[fi].matches[f - MatchStart[fi]].first.pos);
				chains.back().TStart = min(chains.back().TStart, FragInput[fi].matches[f - MatchStart[fi]].second.pos);
			}
			//
			// If this chain overlap with read greater than 0.5%, insert it to chains
			//
			// assert(qEnd - qStart > 10);
			if (chains.back().size() >= 5 and chains.back().QEnd > chains.back().QStart 
				and ((float)(chains.back().QEnd - chains.back().QStart)/read.length) > 0.005 
				and chains.back().QEnd - chains.back().QStart >= 1000) {
				if (chains.size() > 1 and chains.back().OverlapsOnT(chains[0].TStart, chains[0].TEnd, 0.1f)) {
					chains.back().NumOfAnchors0 = chains.back().chain.size();
					chains.back().FirstSDPValue = fragments_valueOrder[fv];				
				}
				else if (chains.size() > 1) {chains.pop_back();}
				else {
					chains.back().NumOfAnchors0 = chains.back().chain.size();
					chains.back().FirstSDPValue = fragments_valueOrder[fv];						
				}
			}
			else {chains.pop_back();}			
		}
		fv++;
	}
}
//
// The input for this function is SplitChain and vector<Cluster_SameDiag *>  ExtendClusters
// Tackle with fragments of different lengths;
// This SDP needs to insert 4 points for only anchors in the overlapping region between Clusters;
// This SDP needs to increase the cost for linking 2 anchors of different directions;
//
// Only insert s1, e1 for forward matches and s2, e2 for reverse matches
//
int SparseDP (SplitChain &inputChain, vector<Cluster_SameDiag *> &FragInput, FinalChain &finalchain, Options &opts, 
			 const vector<float> &LookUpTable, Read &read) {
	if (read.unaligned) return 0;
	if (inputChain.size() == 0) return 0;
	//
	// Compute the matches start in each Cluster;
	//
	vector<int> MatchStart(inputChain.size(), 0);
	int totalMatch = 0;
	ComputeMatchStart<Cluster_SameDiag>(MatchStart, totalMatch, FragInput, inputChain);
	// float rate = 2;//2
	//if ((float)totalMatch / (float) read.length <= 0.005) rate = 4; //2 0.01
	//cerr << "totalMatch/read.length: " << (float)totalMatch / (float) read.length << " rate: " << rate << endl;

	//
	// get points from FragInput and store them in H1;
	//
	vector<Point> H1;		
	int next, prev;
	int totalanchors = 0;
	for (int c = 0; c < inputChain.size(); c++) {

		int cm = inputChain[c];
		FragInput[cm]->matchStart = MatchStart[c]; // set the matchStart for each Cluster
		totalanchors += FragInput[cm]->size();
		prev = c - 1;
		next = c + 1;
		int curstrand = FragInput[cm]->strand;
		int prevstrand, nextstrand;
		//
		// Insert points into H1;
		//
		for (unsigned int i = 0; i < FragInput[cm]->size(); i++) {
			int add = 0;
			// If the current anchor is inside the non-overlapping region, then insert two points for it;
			// If it is forward stranded, then insert s1 and s2 for it;
			if (curstrand == 0) { 
				insertPointsPair(H1, i + MatchStart[c], FragInput[cm]->GetqStart(i),  FragInput[cm]->GettStart(i), FragInput[cm]->length(i), cm, 0, 1);
			}
			else { 
				insertPointsPair(H1, i + MatchStart[c], FragInput[cm]->GetqStart(i), FragInput[cm]->GettStart(i), FragInput[cm]->length(i), cm, 1, 0);	
			}
		}
	}
	//if (totalanchors >= 100000) cerr << "totalanchors: " << totalanchors << " read.name: " << read.name << endl;

	sort(H1.begin(), H1.end(), SortByRowOp<Point>()); // with same q and t coordinates, end point < start point
	vector<unsigned int> H2(H1.size());
	iota(H2.begin(), H2.end(), 0);
	sort(H2.begin(), H2.end(), SortByColOp<Point, unsigned int>(H1));
	
	vector<info> Row;
	vector<info> Col;
	GetRowInfo(H1, Row);
	GetColInfo(H1, H2, Col);

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

	DivideSubProbByRow1(H1, Row, 0, Row.size(), n1, SubR1, eeR1);
	DivideSubProbByCol1(H1, H2, Col, 0, Col.size(), m1, SubC1, eeC1);
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
			// int mi = H1[tt].matchstartNum; // inputChain index

			if (H1[tt].ind == 1 and H1[tt].inv == 1) { //H1[tt] is a start point (s1)
				Value[ii].SS_B_R1 = Row[t].SS_B1;
				Value[ii].counter_B_R1 = Row[t].SS_B1.size();
				Value[ii].val = FragInput[fi]->length(ii - FragInput[fi]->matchStart) * opts.second_anchor_rate;//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]]*rate;
				Value[ii].clusterNum = fi;
				// Value[ii].matchstartNum = mi;
				//Value[ii].val = min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 0 and H1[tt].inv == 1) { // H1[tt] is an end point (e1)
				Value[ii].SS_A_R1 = Row[t].SS_A1;
				Value[ii].counter_A_R1 = Row[t].SS_A1.size();
				//Value[ii].val = min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 1 and H1[tt].inv == 0) { //H1[tt] is a start point (s2)
				Value[ii].SS_B_R2 = Row[t].SS_B2;
				Value[ii].counter_B_R2 = Row[t].SS_B2.size();
				Value[ii].val = FragInput[fi]->length(ii - FragInput[fi]->matchStart) * opts.second_anchor_rate;//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]]*rate;
				Value[ii].clusterNum = fi;
				// Value[ii].matchstartNum = mi;
				//Value[ii].val = min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[tt].orient;
			}
			else { // H1[tt] is an end point (e2)
				Value[ii].SS_A_R2 = Row[t].SS_A2;
				Value[ii].counter_A_R2 = Row[t].SS_A2.size();
				//Value[ii].val = min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
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
			// int mi = H1[H2[tt]].matchstartNum;
			if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 1) { //H1[H2[tt]] a start point (s1)
				Value[ii].SS_B_C1 = Col[t].SS_B1;
				Value[ii].counter_B_C1 = Col[t].SS_B1.size();
				Value[ii].val = FragInput[fi]->length(ii - FragInput[fi]->matchStart) * opts.second_anchor_rate;//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]]*rate;
				Value[ii].clusterNum = fi;
				// Value[ii].matchstartNum = mi;
				//Value[ii].val = min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 0 and H1[H2[tt]].inv == 1) { // H1[H2[tt]] is an end point (e1)
				Value[ii].SS_A_C1 = Col[t].SS_A1;
				Value[ii].counter_A_C1 = Col[t].SS_A1.size();
				//Value[ii].val = min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 0) { //H1[H2[tt]] a start point (s2)
				Value[ii].SS_B_C2 = Col[t].SS_B2;
				Value[ii].counter_B_C2 = Col[t].SS_B2.size();
				Value[ii].val = FragInput[fi]->length(ii - FragInput[fi]->matchStart) * opts.second_anchor_rate;//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]]*rate;
				Value[ii].clusterNum = fi;
				// Value[ii].matchstartNum = mi;
				//Value[ii].val = min(FragInput[ii].qEnd - FragInput[ii].qStart, FragInput[ii].tEnd - FragInput[ii].tStart) * rate;
				Value[ii].orient = H1[H2[tt]].orient;				
			}
			else { // H1[H2[tt]] is an end point (e2)
				Value[ii].SS_A_C2 = Col[t].SS_A2;
				Value[ii].counter_A_C2 = Col[t].SS_A2.size();
				//Value[ii].val = min(FragInput[ii].qEnd - FragInput[ii].qStart + 1, FragInput[ii].tEnd - FragInput[ii].tStart + 1) * rate;
				//Value[ii].orient = H1[H2[tt]].orient;			
			}
		}
	}

	//cerr << "Value: " << Value << endl;
	//cerr << "ProcessPoint\n";
	// finalchain.InitializeOtherParts (MatchStart, totalMatch, Value);
	ProcessPoint<Cluster_SameDiag>(H1, Row, SubR1, SubC1, SubR2, SubC2, Value, opts, LookUpTable, FragInput, opts.second_anchor_rate); 
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
	vector<unsigned int> chain;
	finalchain.SecondSDPValue = max_value;
	vector<bool> link;
	TraceBack(SubR1, SubC1, SubR2, SubC2, Value, max_pos, chain, link); // NOTICE: This chain is from the last anchors to the first anchor;
	finalchain.Initialize(chain, Value);

	// Clear SubR and SubC
	chain.clear();
	clear(H1, H2, SubR1, SubC1, SubR2, SubC2, Value, Row, Col);
	return 0;
}

//
// The input for this function is splitclusters
// fragments of different lengths
// This SDP needs to insert 4 points for any anchors
//
int SparseDP (vector<Cluster> & FragInput, vector<Primary_chain> & Primary_chains, Options & opts, const vector<float> & LookUpTable, Read & read, float & rate) {
	if (FragInput.size() == 0) return 0;
	// if (Primary_chains.size() != 0 and Primary_chains[0].chains.size() == opts.NumAln) return 0;
	vector<Point>  H1;
	// FragInput is vector<Cluster>
	// get points from FragInput and store them in H1		
	for (unsigned int i = 0; i < FragInput.size(); i++) {
	
		// insert start point s1 into H1
		Point s1;
		H1.push_back(s1);
		H1.back().ind = 1; // start
		H1.back().inv = 1; // forward direction
		H1.back().frag_num = i;
		H1.back().se.first = FragInput[i].qStart + 1; 
		H1.back().se.second = FragInput[i].tStart + 1;	
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
		H1.back().se.first = FragInput[i].qEnd - 1; // -1 for chaining overlapped anchors
		H1.back().se.second = FragInput[i].tEnd - 1; // -1 for chaining overlapped anchors
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
		H1.back().se.first = FragInput[i].qStart + 1; 
		H1.back().se.second = FragInput[i].tEnd - 1; // -1 for chaining overlapped anchors	 
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
		H1.back().se.first = FragInput[i].qEnd - 1; // -1 for chaining overlapped anchors
		H1.back().se.second = FragInput[i].tStart + 1;	
		if (FragInput[i].strand == 0) {
			H1.back().orient = 1; 
		}
		else {
			H1.back().orient = 0; 				
		}
	}

	//clock_t begin = clock();

	//Sort the point by row
	if (H1.size() == 0) return 0;
	sort(H1.begin(), H1.end(), SortByRowOp<Point>()); // with same q and t coordinates, end point < start point
	vector<unsigned int> H2(H1.size());
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

	DivideSubProbByRow1(H1, Row, 0, Row.size(), n1, SubR1, eeR1);
	DivideSubProbByCol1(H1, H2, Col, 0, Col.size(), m1, SubC1, eeC1);
	DivideSubProbByRow2(H1, Row, 0, Row.size(), n2, SubR2, eeR2);	
	DivideSubProbByCol2(H1, H2, Col, 0, Col.size(), m2, SubC2, eeC2);

	//
	// Get SS_A_R1, SS_B_R1, SS_A_R2 and SS_B_R2 for each fragment
	//
	vector<Fragment_Info> Value(FragInput.size());
	for (unsigned int t = 0; t < Row.size(); ++t) {

		for (unsigned int tt = Row[t].pstart; tt < Row[t].pend; ++tt) {

			unsigned int ii = H1[tt].frag_num;

			if (H1[tt].ind == 1 and H1[tt].inv == 1) { //H1[tt] is a start point (s1)
				Value[ii].SS_B_R1 = Row[t].SS_B1;
				Value[ii].counter_B_R1 = Row[t].SS_B1.size();
				Value[ii].val = FragInput[ii].Val*rate;
				Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 0 and H1[tt].inv == 1) { // H1[tt] is an end point (e1)
				Value[ii].SS_A_R1 = Row[t].SS_A1;
				Value[ii].counter_A_R1 = Row[t].SS_A1.size();
			}
			else if (H1[tt].ind == 1 and H1[tt].inv == 0) { //H1[tt] is a start point (s2)
				Value[ii].SS_B_R2 = Row[t].SS_B2;
				Value[ii].counter_B_R2 = Row[t].SS_B2.size();
				Value[ii].val = FragInput[ii].Val*rate;
				Value[ii].orient = H1[tt].orient;
			}
			else { // H1[tt] is an end point (e2)
				Value[ii].SS_A_R2 = Row[t].SS_A2;
				Value[ii].counter_A_R2 = Row[t].SS_A2.size();
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
				Value[ii].val = FragInput[ii].Val*rate;
				Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 0 and H1[H2[tt]].inv == 1) { // H1[H2[tt]] is an end point (e1)
				Value[ii].SS_A_C1 = Col[t].SS_A1;
				Value[ii].counter_A_C1 = Col[t].SS_A1.size();
			}
			else if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 0) { //H1[H2[tt]] a start point (s2)
				Value[ii].SS_B_C2 = Col[t].SS_B2;
				Value[ii].counter_B_C2 = Col[t].SS_B2.size();
				Value[ii].val = FragInput[ii].Val*rate;
				Value[ii].orient = H1[H2[tt]].orient;				
			}
			else { // H1[H2[tt]] is an end point (e2)
				Value[ii].SS_A_C2 = Col[t].SS_A2;
				Value[ii].counter_A_C2 = Col[t].SS_A2.size();
			}

		}
	}
	//cerr << "Value: " << Value << endl;
	//cerr << "ProcessPoint\n";
	ProcessPoint<Cluster>(H1, Row, SubR1, SubC1, SubR2, SubC2, Value, opts, LookUpTable, FragInput, rate);	
	DecidePrimaryChains(FragInput, SubR1, SubC1, SubR2, SubC2, Value, Primary_chains, read, opts);

	// Clear SubR and SubC
	clear(H1, H2, SubR1, SubC1, SubR2, SubC2, Value, Row, Col);
	return 0;
}

//
// The input to this SparseDP is pure matches -- insert 4 points for each anchor
//
int SparseDP (vector<Cluster> &FragInput, vector<UltimateChain> &chains, Options &opts, const vector<float> &LookUpTable, Read &read, float rate) {
	if (read.unaligned) return 0;
	//
	// Compute the matches start in each Cluster;
	//
	vector<int> MatchStart(FragInput.size(), 0);
	for (int c = 1; c < FragInput.size(); c++){
		MatchStart[c] = FragInput[c - 1].matches.size();
	}
	for (int c = 1; c < FragInput.size(); c++){
		MatchStart[c] += MatchStart[c - 1];
	}	
	vector<Point>  H1;
	int totalMatch = 0;
	for (int cm = 0; cm < FragInput.size(); cm++) {
		//
		// Insert points into H1;
		//
		for (unsigned int i = 0; i < FragInput[cm].matches.size(); i++) {
			if (FragInput[cm].strand == 0) { 				
				insertPointsPair(H1, i + MatchStart[cm], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, FragInput[cm].matchesLengths[i], cm, 0, 1);
				if (i == 0 or i == FragInput[cm].matches.size() - 1) insertPointsPair(H1, i + MatchStart[cm], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, FragInput[cm].matchesLengths[i], cm, 1, 1);
			}
			else { 
				// insert start point s2 into H1
				insertPointsPair(H1, i + MatchStart[cm], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, FragInput[cm].matchesLengths[i], cm, 1, 0);
				if (i == 0 or i == FragInput[cm].matches.size() - 1) insertPointsPair(H1, i + MatchStart[cm], FragInput[cm].matches[i].first.pos, FragInput[cm].matches[i].second.pos, FragInput[cm].matchesLengths[i], cm, 0, 0);	
			}
		}
		totalMatch += FragInput[cm].matches.size();
	}
	//if (totalanchors >= 100000) cerr << "totalanchors: " << totalanchors << " read.name: " << read.name << endl;

	sort(H1.begin(), H1.end(), SortByRowOp<Point>()); // with same q and t coordinates, end point < start point
	vector<unsigned int> H2(H1.size());
	iota(H2.begin(), H2.end(), 0);
	sort(H2.begin(), H2.end(), SortByColOp<Point, unsigned int>(H1));

	
	vector<info> Row;
	vector<info> Col;
	GetRowInfo(H1, Row);
	GetColInfo(H1, H2, Col);

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

	DivideSubProbByRow1(H1, Row, 0, Row.size(), n1, SubR1, eeR1);
	DivideSubProbByCol1(H1, H2, Col, 0, Col.size(), m1, SubC1, eeC1);
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
			// int mi = H1[tt].matchstartNum; // inputChain index
			if (H1[tt].ind == 1 and H1[tt].inv == 1) { //H1[tt] is a start point (s1)
				Value[ii].SS_B_R1 = Row[t].SS_B1;
				Value[ii].counter_B_R1 = Row[t].SS_B1.size();
				Value[ii].val = FragInput[fi].matchesLengths[ii - MatchStart[fi]] * rate;
				Value[ii].clusterNum = fi;
				Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 0 and H1[tt].inv == 1) { // H1[tt] is an end point (e1)
				Value[ii].SS_A_R1 = Row[t].SS_A1;
				Value[ii].counter_A_R1 = Row[t].SS_A1.size();
			}
			else if (H1[tt].ind == 1 and H1[tt].inv == 0) { //H1[tt] is a start point (s2)
				Value[ii].SS_B_R2 = Row[t].SS_B2;
				Value[ii].counter_B_R2 = Row[t].SS_B2.size();
				Value[ii].val = FragInput[fi].matchesLengths[ii - MatchStart[fi]] * rate;//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]]*rate;
				Value[ii].clusterNum = fi;
				Value[ii].orient = H1[tt].orient;
			}
			else { // H1[tt] is an end point (e2)
				Value[ii].SS_A_R2 = Row[t].SS_A2;
				Value[ii].counter_A_R2 = Row[t].SS_A2.size();
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
			// int mi = H1[H2[tt]].matchstartNum;
			if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 1) { //H1[H2[tt]] a start point (s1)
				Value[ii].SS_B_C1 = Col[t].SS_B1;
				Value[ii].counter_B_C1 = Col[t].SS_B1.size();
				Value[ii].val = FragInput[fi].matchesLengths[ii - MatchStart[fi]] * rate;//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]]*rate;
				Value[ii].clusterNum = fi;
				Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 0 and H1[H2[tt]].inv == 1) { // H1[H2[tt]] is an end point (e1)
				Value[ii].SS_A_C1 = Col[t].SS_A1;
				Value[ii].counter_A_C1 = Col[t].SS_A1.size();
			}
			else if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 0) { //H1[H2[tt]] a start point (s2)
				Value[ii].SS_B_C2 = Col[t].SS_B2;
				Value[ii].counter_B_C2 = Col[t].SS_B2.size();
				Value[ii].val = FragInput[fi].matchesLengths[ii - MatchStart[fi]] * rate;//FragInput[inputChain[mi]]->matchesLengths[ii - MatchStart[mi]]*rate;
				Value[ii].clusterNum = fi;
				Value[ii].orient = H1[H2[tt]].orient;				
			}
			else { // H1[H2[tt]] is an end point (e2)
				Value[ii].SS_A_C2 = Col[t].SS_A2;
				Value[ii].counter_A_C2 = Col[t].SS_A2.size();
			}

		}
	}

	ProcessPoint<Cluster>(H1, Row, SubR1, SubC1, SubR2, SubC2, Value, opts, LookUpTable, FragInput, MatchStart, rate); 
	DecidePrimaryChains(FragInput, SubR1, SubC1, SubR2, SubC2, Value, chains, read, opts, MatchStart);
	for (int c = 0; c < chains.size(); c++) {
		chains[c].ClusterIndex.resize(chains[c].chain.size());
		for (int s = 0; s < chains[c].chain.size(); s++) {
			int ii = chains[c].chain[s];
			chains[c].chain[s] -= MatchStart[Value[ii].clusterNum];
			chains[c].ClusterIndex[s] = Value[ii].clusterNum;
		}
	}
	// Clear SubR and SubC
	clear(H1, H2, SubR1, SubC1, SubR2, SubC2, Value, Row, Col);
	return 0;
}

//
// The input to this SparseDP is for one cluster
//
int SparseDP (int ClusterIndex, vector<Cluster> &FragInput, UltimateChain &ultimatechain, Options &opts, const vector<float> &LookUpTable, Read &read) {
	if (read.unaligned) return 0;
	if (FragInput[ClusterIndex].matches.size() == 0) return 0;
	//
	// get points from FragInput and store them in H1;
	//
	vector<Point> H1;		
	int next, prev;
	FragInput[ClusterIndex].matchStart = 0;
	vector<int> MatchStart(FragInput.size(), 0);
	for (int i = 0; i < FragInput[ClusterIndex].matches.size(); i++) {
		if (FragInput[ClusterIndex].strand == 0) { // Insert s1, e1 into H1;
			insertPointsPair(H1, i, FragInput[ClusterIndex].matches[i].first.pos, FragInput[ClusterIndex].matches[i].second.pos, 
				FragInput[ClusterIndex].matchesLengths[i], ClusterIndex, 0, 1);
		}
		else { // insert start point s2, e2 into H1
			insertPointsPair(H1, i, FragInput[ClusterIndex].matches[i].first.pos, FragInput[ClusterIndex].matches[i].second.pos, 
				FragInput[ClusterIndex].matchesLengths[i], ClusterIndex, 1, 0);
		}
	}
	int totalMatch = FragInput[ClusterIndex].matches.size();
	//if (totalanchors >= 100000) cerr << "totalanchors: " << totalanchors << " read.name: " << read.name << endl;

	//Sort the point by row
	sort(H1.begin(), H1.end(), SortByRowOp<Point>()); // with same q and t coordinates, end point < start point
	vector<unsigned int> H2(H1.size());
	iota(H2.begin(), H2.end(), 0);
	sort(H2.begin(), H2.end(), SortByColOp<Point, unsigned int>(H1));

	
	vector<info> Row;
	vector<info> Col;
	GetRowInfo(H1, Row);
	GetColInfo(H1, H2, Col);

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

	DivideSubProbByRow1(H1, Row, 0, Row.size(), n1, SubR1, eeR1);
	DivideSubProbByCol1(H1, H2, Col, 0, Col.size(), m1, SubC1, eeC1);
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
			// int mi = H1[tt].matchstartNum; // inputChain index

			if (H1[tt].ind == 1 and H1[tt].inv == 1) { //H1[tt] is a start point (s1)
				Value[ii].SS_B_R1 = Row[t].SS_B1;
				Value[ii].counter_B_R1 = Row[t].SS_B1.size();
				Value[ii].val = FragInput[fi].matchesLengths[ii - MatchStart[fi]] * opts.second_anchor_rate;
				Value[ii].clusterNum = fi;
				Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 0 and H1[tt].inv == 1) { // H1[tt] is an end point (e1)
				Value[ii].SS_A_R1 = Row[t].SS_A1;
				Value[ii].counter_A_R1 = Row[t].SS_A1.size();
			}
			else if (H1[tt].ind == 1 and H1[tt].inv == 0) { //H1[tt] is a start point (s2)
				Value[ii].SS_B_R2 = Row[t].SS_B2;
				Value[ii].counter_B_R2 = Row[t].SS_B2.size();
				Value[ii].val = FragInput[fi].matchesLengths[ii - MatchStart[fi]] * opts.second_anchor_rate;
				Value[ii].clusterNum = fi;
				Value[ii].orient = H1[tt].orient;
			}
			else { // H1[tt] is an end point (e2)
				Value[ii].SS_A_R2 = Row[t].SS_A2;
				Value[ii].counter_A_R2 = Row[t].SS_A2.size();
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
			// int mi = H1[H2[tt]].matchstartNum;

			if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 1) { //H1[H2[tt]] a start point (s1)
				Value[ii].SS_B_C1 = Col[t].SS_B1;
				Value[ii].counter_B_C1 = Col[t].SS_B1.size();
				Value[ii].val = FragInput[fi].matchesLengths[ii - MatchStart[fi]] * opts.second_anchor_rate;
				Value[ii].clusterNum = fi;
				Value[ii].orient = H1[H2[tt]].orient;
			}
			else if (H1[H2[tt]].ind == 0 and H1[H2[tt]].inv == 1) { // H1[H2[tt]] is an end point (e1)
				Value[ii].SS_A_C1 = Col[t].SS_A1;
				Value[ii].counter_A_C1 = Col[t].SS_A1.size();
			}
			else if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 0) { //H1[H2[tt]] a start point (s2)
				Value[ii].SS_B_C2 = Col[t].SS_B2;
				Value[ii].counter_B_C2 = Col[t].SS_B2.size();
				Value[ii].val = FragInput[fi].matchesLengths[ii - MatchStart[fi]] * opts.second_anchor_rate;
				Value[ii].clusterNum = fi;
				Value[ii].orient = H1[H2[tt]].orient;				
			}
			else { // H1[H2[tt]] is an end point (e2)
				Value[ii].SS_A_C2 = Col[t].SS_A2;
				Value[ii].counter_A_C2 = Col[t].SS_A2.size();
			}

		}
	}

	//cerr << "Value: " << Value << endl;
	//cerr << "ProcessPoint\n";
	ProcessPoint<Cluster>(H1, Row, SubR1, SubC1, SubR2, SubC2, Value, opts, LookUpTable, FragInput, MatchStart, opts.second_anchor_rate); 
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

	vector<unsigned int> chain; vector<bool> link;
	ultimatechain.FirstSDPValue = max_value;
	TraceBack(SubR1, SubC1, SubR2, SubC2, Value, max_pos, chain, link); // NOTICE: This chain is from the last anchors to the first anchor;
	ultimatechain.Initialize(chain, Value, link);
	ultimatechain.NumOfAnchors0 = chain.size();
	chain.clear(); link.clear();
	clear(H1, H2, SubR1, SubC1, SubR2, SubC2, Value, Row, Col);
	return 0;
}

#endif
