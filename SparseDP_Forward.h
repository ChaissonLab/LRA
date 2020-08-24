// This program is implemented Sparse Dynamic Programming algorithm described in David Eppstein paper
// Author: Jingwen Ren

#ifndef SPARSE_DP_FORWARD_
#define SPARSE_DP_FORWARD_

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


// Note: Each fragment has the same length
void 
ProcessPoint_ForwardOnly (const std::vector<Point> & H1, const vector<int> &MatchLengths, std::vector<info> & V, StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, 
							std::vector<Fragment_Info> & Value, Options & opts, const std::vector<float> & LookUpTable, int rate) {

//ProcessPoint (const std::vector<Point> & H1, const std::vector<unsigned int> & H3, std::vector<info> & V, StackOfSubProblems & SubR, StackOfSubProblems & SubC,
//				  std::vector<Fragment_Info> & Value, Options & opts, const std::vector<float> & LookUpTable, int rate) {

	bool step_sdp = 1;
	for (unsigned int i = 0; i < H1.size(); ++i) { // process points by row

		long int ForwardDiag = static_cast<long int>(H1[i].se.second) - static_cast<long int>(H1[i].se.first);

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
						Maximization (SubR1[j].now, SubR1[j].last, SubR1[j].Di, SubR1[j].Ei, SubR1[j].Dv, SubR1[j].Db, SubR1[j].Block, SubR1[j].S_1, LookUpTable, opts, step_sdp); // TODO(Jingwen) anything change for SubC????
						SubR1[j].last = SubR1[j].Eb[*t];

						//cerr << "retrieve the value from the maximization structure\n";
						unsigned int i1 = *t; // i1 stores the index in Ei that ForwardDiag is in
						unsigned int i2; // i2 stores the index in Di which is the best candidate for ForwardDiag
						FindValueInBlock(ForwardDiag, SubR1[j].S_1, SubR1[j].Ei, SubR1[j].Block, i1, i2);
						//cerr << "the index in Ei that ForwardDiag is in----i1: " << i1 << "\n";
						//cerr << "the index in Di which is the best candidate for ForwardDiag ---- i2: " << i2 << "\n";
						

						SubR1[j].Ev[i1] = SubR1[j].Dv[i2] + w(SubR1[j].Di[i2], SubR1[j].Ei[i1], LookUpTable, opts, step_sdp) + MatchLengths[ii] * rate; //opts.globalK * rate; 
						SubR1[j].Ep[i1] = i2;							

						//cerr << "SubR1[" << j << "].Ev[" << i1 << "]: " << SubR1[j].Ev[i1] << ", SubR1[" << j << "].Ep[" << i1 << "]: " << SubR1[j].Ep[i1] << "\n"; 


						// Update the value of this point
						//TODO(Jingwen): if this point is a s1 of a reverse orientated anchor, only update the value when  SubR1[j].Dp[Ep[i1]] points to a forward orientated anchor

						//TODO(Jingwen): only for debug
						assert(Value[ii].orient == H1[i].orient);
						int p = SubR1[j].Dp[SubR1[j].Ep[i1]];
						// if the current anchor is reverse oriented and the previous anchor is also reverse oriented, then do not update Value[ii]
						// Otherwise update Value[ii]
						if (! (Value[ii].orient == 0 and Value[p].orient == 0)) { 
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

						SubC1[j].Ev[i1] = SubC1[j].Dv[i2] + w(SubC1[j].Di[i2], SubC1[j].Ei[i1], LookUpTable, opts, step_sdp) + MatchLengths[ii] * rate; //opts.globalK * rate; 
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
						if (! (Value[ii].orient == 0 and Value[p].orient == 0)) {
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
						}
						--Value[ii].counter_B_C1;	
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
	}
}


// void 
// TraceBack_ForwardOnly (StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, const vector<Fragment_Info> & Value, 
// 							unsigned int i, vector<unsigned int> & Chain) {

// 	long int prev_sub = Value[i].prev_sub;
// 	long int prev_ind = Value[i].prev_ind;
// 	Chain.push_back(i);

// 	if (prev_sub != -1 and prev_ind != -1) {
		
// 		// if (Value[i].prev == 1 and Value[i].inv == 1) { // The previous subproblem is SubR1
// 		// 	unsigned int ind = SubR1[prev_sub].Ep[prev_ind];
// 		// 	TraceBack_ForwardOnly(SubR1, SubC1, Value, SubR1[prev_sub].Dp[ind], FinalChain);
// 		// }
// 		if (Value[i].prev == 1 and Value[i].inv == 1) { // The previous subproblem is SubR1
// 			assert(prev_sub <  SubR1.StackSub.size());
// 			assert(prev_ind < SubR1[prev_sub].Ep.size());
// 			unsigned int ind = SubR1[prev_sub].Ep[prev_ind];
// 			assert(ind < SubR1[prev_sub].Dp.size());
// 			i = SubR1[prev_sub].Dp[ind];
// 			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubR1[prev_sub].Dp[ind], Chain);
// 		}
// 		else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
// 			assert(prev_sub <  SubC1.StackSub.size());
// 			assert(prev_ind < SubC1[prev_sub].Ep.size());
// 			unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
// 			assert(ind < SubC1[prev_sub].Dp.size());
// 			i = SubC1[prev_sub].Dp[ind];
// 			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC1[prev_sub].Dp[ind], Chain);
// 		}	
// 		// else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
// 		// 	unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
// 		// 	TraceBack_ForwardOnly(SubR1, SubC1, Value, SubC1[prev_sub].Dp[ind], FinalChain);
// 		// }	
// 		prev_sub = Value[i].prev_sub;
// 		prev_ind = Value[i].prev_ind;		
// 		Chain.push_back(i);
// 	}
// }


//
// This function is for tracing back a chain;
//
void 
TraceBack_ForwardOnly (StackOfSubProblems & SubR1, StackOfSubProblems & SubC1, const vector<Fragment_Info> & Value, 
							unsigned int i, vector<unsigned int> & Chain) {

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
		else if (Value[i].prev == 0 and Value[i].inv == 1) { // The previous subproblem is SubC1
			assert(prev_sub <  SubC1.StackSub.size());
			assert(prev_ind < SubC1[prev_sub].Ep.size());
			unsigned int ind = SubC1[prev_sub].Ep[prev_ind];
			assert(ind < SubC1[prev_sub].Dp.size());
			i = SubC1[prev_sub].Dp[ind];
			//TraceBack(SubR1, SubC1, SubR2, SubC2, Value, SubC1[prev_sub].Dp[ind], Chain);
		}		
		prev_sub = Value[i].prev_sub;
		prev_ind = Value[i].prev_ind;		
		Chain.push_back(i);
	}
}


// The input for this function is GenomePairs which is from gapPairs (from snd SDP)
// Each fragment has the same length
//
int SparseDP_ForwardOnly (const GenomePairs &FragInput, const vector<int> &MatchLengths, std::vector<unsigned int> &chain, Options &opts, 
							const std::vector<float> &LookUpTable, float &inv_value, int &inv_NumOfAnchors, int rate = 5) {
	
	if (FragInput.size() == 0) return 0;
	std::vector<Point> H1;
	// FragInput is vector<GenomePair>
	// get points from FragInput and store them in H1
	for (unsigned int i = 0; i < FragInput.size(); i++) { 

			// insert start point s1 into H1
			Point s1;
			H1.push_back(s1);
			H1.back().ind = 1; // start
			H1.back().inv = 1; // forward direction
			H1.back().frag_num = i;
			H1.back().se.first = FragInput[i].first.pos; 
			H1.back().se.second = FragInput[i].second.pos;	
			H1.back().orient = 1; // the point comes from a forward oriented anchor

			// insert end point e1 into H1
			Point e1;
			H1.push_back(e1);
			H1.back().ind = 0; // end
			H1.back().inv = 1; // forward direction		
			H1.back().frag_num = i;
			H1.back().se.first = FragInput[i].first.pos + MatchLengths[i];
			H1.back().se.second = FragInput[i].second.pos + MatchLengths[i];					
		
	}
	
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

	//std::vector<Subproblem> SubR;
	//std::vector<Subproblem> SubC;
	StackOfSubProblems SubR1;
	StackOfSubProblems SubC1;
	int eeR1 = 0, eeC1 = 0;

	//cerr << "DivideSubByRow\n";
	DivideSubProbByRow1(H1, Row, 0, Row.size(), n1, SubR1, eeR1);
	//cerr << "SubR: " << SubR << endl;

	//cerr << "DivideSubByCol\n";
	DivideSubProbByCol1(H1, H2, Col, 0, Col.size(), m1, SubC1, eeC1);
	//cerr << "SubC: " << SubC << endl;

	// Get SS_A_R1, SS_B_R1 for each fragment
	std::vector<Fragment_Info> Value(FragInput.size());
	for (unsigned int t = 0; t < Row.size(); ++t) {
		for (unsigned int tt = Row[t].pstart; tt < Row[t].pend; ++tt) {

			unsigned int ii = H1[tt].frag_num;

			if (H1[tt].ind == 1 and H1[tt].inv == 1) { //H1[tt] is a start point (s1)
				Value[ii].SS_B_R1 = Row[t].SS_B1;
				Value[ii].counter_B_R1 = Row[t].SS_B1.size();
				Value[ii].val = MatchLengths[ii]*rate; //opts.globalK * rate;
				Value[ii].orient = H1[tt].orient;
			}
			else if (H1[tt].ind == 0 and H1[tt].inv == 1) { // H1[tt] is an end point (e1)
				Value[ii].SS_A_R1 = Row[t].SS_A1;
				Value[ii].counter_A_R1 = Row[t].SS_A1.size();
				//Value[ii].val = (opts.globalK); // TODO(Jingwen): make sure that this is redundant
			}
		}
	}


	// Get SS_A_C2 and SS_B_C2 for each fragment
	for (unsigned int t = 0; t < Col.size(); ++t) {
		for (unsigned int tt = Col[t].pstart; tt < Col[t].pend; ++tt) { 

			unsigned int ii = H1[H2[tt]].frag_num;

			if (H1[H2[tt]].ind == 1 and H1[H2[tt]].inv == 1) { //H1[H2[tt]] a start point (s1)
				Value[ii].SS_B_C1 = Col[t].SS_B1;
				Value[ii].counter_B_C1 = Col[t].SS_B1.size();
				Value[ii].val = MatchLengths[ii]*rate; // opts.globalK * rate;
			}
			else if (H1[H2[tt]].ind == 0 and H1[H2[tt]].inv == 1) { // H1[H2[tt]] is an end point (e1)
				Value[ii].SS_A_C1 = Col[t].SS_A1;
				Value[ii].counter_A_C1 = Col[t].SS_A1.size();
				//Value[ii].val = (opts.globalK);
			}
		}
	}


	//cerr << "Value: " << Value << endl;
	

	//cerr << "ProcessPoint\n";

	ProcessPoint_ForwardOnly(H1, MatchLengths, Row, SubR1, SubC1, Value, opts, LookUpTable, rate);

	
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
	
	inv_value = max_value;
	inv_NumOfAnchors = chain.size();
	//cerr << "TraceBack\n";
	// Trace back to get the FinalChain
	// store the index of points
	chain.clear();
	TraceBack_ForwardOnly(SubR1, SubC1, Value, max_pos, chain);
	//std::reverse(chain.begin(), chain.end());

	// Clear SubR and SubC
	SubR1.Clear(eeR1);
	SubC1.Clear(eeC1);

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
