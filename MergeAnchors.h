// Merge Anchors 
#ifndef MERGE_ANCHORS_H_
#define MERGE_ANCHORS_H_

#include <iostream> //std::cout 
#include <fstream>   
#include <cstdlib>   // std::labs, std::EXIT FAILURE, std::EXIT SUCCESS
#include <string> 
#include <cmath>        // std::labs
#include <cstdio>    // std::FILE std::perror
#include <vector>
#include <numeric>
#include <utility>
#include <set>
#include <list>
#include "Sorting.h"
#include "Clustering.h"
#include "Options.h"


int 
gapdifference (int &e, GenomePairs &matches, int &strand) {
	int prevDiag, curDiag;
	if (strand == 0) {
		prevDiag = matches[e-1].second.pos - matches[e-1].first.pos;
		curDiag = matches[e].second.pos - matches[e].first.pos;		
	}
	else {
		prevDiag = matches[e-1].second.pos + matches[e-1].first.pos;
		curDiag = matches[e].second.pos + matches[e].first.pos;	
		//cerr << "prevDiag: " << prevDiag << endl;
		//cerr << "curDiag: " << curDiag << endl;
		//cerr << "prevDiag - curDiag: " << prevDiag - curDiag << endl;
	}
	return std::abs(prevDiag - curDiag);
}


// TODO(Jingwen): check all the inequility
void 
MergeAnchors (Options & opts, vector<Cluster> &refinedClusters, vector<LogCluster> &refinedLogClusters, 
					int r, vector<ClusterCoordinates> &mergedAnchors) {

	GenomePos qBoundary_s = 0, qBoundary_e = 0, tBoundary_s = 0, tBoundary_e = 0, next_qBoundary_s = 0, next_tBoundary_s = 0;


	for (int l = refinedLogClusters[r].SubCluster.size() - 1; l >= 0; --l) {
		//cerr << "l " << l << " refinedLogClusters[r].SubCluster[l].strand: " << refinedLogClusters[r].SubCluster[l].strand << endl;
		if (refinedLogClusters[r].SubCluster.size() != 1) {
			GenomePos cur_qStart = refinedLogClusters[r].SubCluster[l].qStart;
			GenomePos cur_qEnd = refinedLogClusters[r].SubCluster[l].qEnd;
			GenomePos cur_tStart = refinedLogClusters[r].SubCluster[l].tStart;
			GenomePos cur_tEnd = refinedLogClusters[r].SubCluster[l].tEnd;

			if ( l > 0) {
				GenomePos next_qStart = refinedLogClusters[r].SubCluster[l-1].qStart;
				GenomePos next_qEnd = refinedLogClusters[r].SubCluster[l-1].qEnd;
				GenomePos next_tStart = refinedLogClusters[r].SubCluster[l-1].tStart;
				GenomePos next_tEnd = refinedLogClusters[r].SubCluster[l-1].tEnd;


				if (next_qStart < cur_qEnd) {
					qBoundary_e = next_qStart; // qBoundary_e is inclusive
					next_qBoundary_s = cur_qEnd; // qBoundary_s is exclusive
				}
				else {
					qBoundary_e = 0;
					next_qBoundary_s = 0;
				}

				if (next_tStart < cur_tEnd) {
					tBoundary_e = next_tStart;
					next_tBoundary_s = cur_tEnd;
				}
				else {
					tBoundary_e = 0;
					next_tBoundary_s = 0;
				}				
			}
			else {
				qBoundary_e = 0;
				tBoundary_e	= 0;			
			}
		}
		//cerr << "qBoundary_s: " << qBoundary_s << "  tBoundary_s: " << tBoundary_s << "     qBoundary_e: " << qBoundary_e << "  tBoundary_e: " << tBoundary_e << endl;
		//cerr << "next_qBoundary_s: " << next_qBoundary_s << "  next_tBoundary_s: " << next_tBoundary_s << endl;

		//if (qBoundary_e != 0 or tBoundary_e != 0 or qBoundary_s != 0 or tBoundary_s != 0) {

			//DiagonalSort<GenomeTuple>(refinedClusters[r].matches.begin() + refinedLogClusters[r].SubCluster[l].start, 
			//							refinedClusters[r].matches.begin()+refinedLogClusters[r].SubCluster[l].end);
			//CleanOffDiagonal(refinedClusters[r].matches, refinedLogClusters[r].SubCluster[l].start,refinedLogClusters[r].SubCluster[l].end, 
			//						opts, refinedLogClusters[r].SubCluster[l].strand);

			// there is boudary

		GenomePos lastQ = 0, lastT = 0; // lastQ, lastT means the qEnd, tEnd of the last mergedAnchors; 
										// The reason why we use lastQ, lastT is that we do not want two long merged Anchors overlapped with each other
										// SDP cannot choose overlapped merged Anchors
		bool str = 0;
		CartesianSort<GenomeTuple>(refinedClusters[r].matches, refinedLogClusters[r].SubCluster[l].start, refinedLogClusters[r].SubCluster[l].end); 
		int s = refinedLogClusters[r].SubCluster[l].start;
		while (s < refinedLogClusters[r].SubCluster[l].end) {
		//for (int s = refinedLogClusters[r].SubCluster[l].start; s < refinedLogClusters[r].SubCluster[l].end; ++s) {
			str = refinedLogClusters[r].SubCluster[l].strand;
			int e = s;
			GenomePos qStart = refinedClusters[r].matches[s].first.pos, 
					  qEnd = refinedClusters[r].matches[s].first.pos + opts.globalK, 
					  tStart = refinedClusters[r].matches[s].second.pos, 
					  tEnd = refinedClusters[r].matches[s].second.pos + opts.globalK;

			while (	(qBoundary_e == 0 or refinedClusters[r].matches[s].first.pos + opts.globalK <= qBoundary_e) // qBoundary_e is inclusive
					and (tBoundary_e == 0 or refinedClusters[r].matches[s].second.pos + opts.globalK <= tBoundary_e)
					and (qBoundary_s == 0 or refinedClusters[r].matches[s].first.pos >= qBoundary_s) // qBoundary_s is exclusive
					and (tBoundary_s == 0 or refinedClusters[r].matches[s].second.pos >= tBoundary_s)
					and (lastQ == 0 or refinedClusters[r].matches[s].first.pos >= lastQ) // lastQ and lastT are exclusive
					and (lastT == 0 or (str == 1 or refinedClusters[r].matches[s].second.pos >= lastT)) // If this is forward stranded, then the anchor s.second.pos >= lastT 
					and (lastT == 0 or (str == 0 or refinedClusters[r].matches[s].second.pos + opts.globalK - 1 <= lastT))	) { // If this is rev stranded, then the anchor s.second.pos + opts.globalK <= lastT

				assert(s < refinedLogClusters[r].SubCluster[l].end);
				++e;
				int rgap, ggap;
				if (refinedLogClusters[r].SubCluster[l].strand == 0) {
					rgap = refinedClusters[r].matches[e].first.pos - refinedClusters[r].matches[e-1].first.pos - opts.globalK;
					ggap = refinedClusters[r].matches[e].second.pos - refinedClusters[r].matches[e-1].second.pos - opts.globalK;
				}
				else {
					rgap = refinedClusters[r].matches[e].first.pos - refinedClusters[r].matches[e-1].first.pos - opts.globalK;
					ggap = refinedClusters[r].matches[e - 1].second.pos - refinedClusters[r].matches[e].second.pos - opts.globalK;								
				}

				while (	e < refinedLogClusters[r].SubCluster[l].end 
					and (qBoundary_e == 0 or refinedClusters[r].matches[e].first.pos + opts.globalK <= qBoundary_e) // qBoundary_e is inclusive
					and (tBoundary_e == 0 or refinedClusters[r].matches[e].second.pos + opts.globalK <= tBoundary_e)
					and (qBoundary_s == 0 or refinedClusters[r].matches[e].first.pos >= qBoundary_s) // qBoundary_s is exclusive
					and (tBoundary_s == 0 or refinedClusters[r].matches[e].second.pos >= tBoundary_s)
					and gapdifference(e, refinedClusters[r].matches, refinedLogClusters[r].SubCluster[l].strand) < opts.maxDiag 
					and std::max(std::abs(rgap), std::abs(ggap)) <= opts.maxGap	) {

					qStart = min(qStart, refinedClusters[r].matches[e].first.pos);
					qEnd = max(qEnd, refinedClusters[r].matches[e].first.pos + opts.globalK);
					tStart = min(tStart, refinedClusters[r].matches[e].second.pos);
					tEnd = max(tEnd, refinedClusters[r].matches[e].second.pos + opts.globalK);
					++e;
				}

				if (e == refinedLogClusters[r].SubCluster[l].end) {
					// insert
					mergedAnchors.push_back(ClusterCoordinates(s, e, qStart, qEnd, tStart, tEnd, refinedLogClusters[r].SubCluster[l].strand));
					//cerr << "1   s: " << s << " e: " << e << " qStart: " << qStart << "  qEnd: " << qEnd << "  tStart: " << tStart << "  tEnd: " << tEnd << endl;
					s = e;
					break; 
				}
				else {
					// insert
					mergedAnchors.push_back(ClusterCoordinates(s, e, qStart, qEnd, tStart, tEnd, refinedLogClusters[r].SubCluster[l].strand));
					//cerr << "2   s: " << s << " e: " << e << " qStart: " << qStart << "  qEnd: " << qEnd << "  tStart: " << tStart << "  tEnd: " << tEnd << endl;

					if (refinedLogClusters[r].SubCluster[l].strand == 0) {
						lastQ = qEnd;
						lastT = tEnd;
					}
					else {
						lastQ = qEnd;
						lastT = tStart - 1; // make lastT exclusive
					}

					s = e;
					qStart = refinedClusters[r].matches[s].first.pos, 
					qEnd = refinedClusters[r].matches[s].first.pos + opts.globalK, 
					tStart = refinedClusters[r].matches[s].second.pos, 
					tEnd = refinedClusters[r].matches[s].second.pos + opts.globalK;
				}
			}

		/*	if (s < refinedLogClusters[r].SubCluster[l].end and !((refinedClusters[r].matches[s].first.pos + opts.globalK <= qBoundary_e or qBoundary_e == 0) 
					and (refinedClusters[r].matches[s].second.pos + opts.globalK <= tBoundary_e or tBoundary_e == 0)
					and (refinedClusters[r].matches[s].first.pos > qBoundary_s or qBoundary_s == 0)
					and (refinedClusters[r].matches[s].second.pos > tBoundary_s or tBoundary_s == 0)))
		*/
			if (s < refinedLogClusters[r].SubCluster[l].end) { // This IF condition is aimed at anchor s that fails at the previous while condition
				// insert
				mergedAnchors.push_back(ClusterCoordinates(s, s+1, qStart, qEnd, tStart, tEnd, refinedLogClusters[r].SubCluster[l].strand));
				//cerr << "3   s: " << s << " s+1: " << s+1 << " qStart: " << qStart << "  qEnd: " << qEnd << "  tStart: " << tStart << "  tEnd: " << tEnd << endl;
				
				++s;
			}
		}
		//}
		qBoundary_s = next_qBoundary_s;
		tBoundary_s = next_tBoundary_s;
	}		

}

#endif