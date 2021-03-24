#ifndef CHAIN_REFINE_H_
#define CHAIN_REFINE_H_
#include <math.h>
#include "Genome.h"
#include "Read.h"
#include "Options.h"
#include "CompareLists.h"
#include "Sorting.h"
#include "TupleOps.h"
#include "Clustering.h"
#include "TupleOps.h"
#include "Chain.h"
#include "ClusterRefine.h"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cmath>	// std::log 
#include <sstream>
#include <thread>
#include <climits>

void append_to_closetcluster(GenomePairs &matches, int start, int end, Cluster *cluster, Cluster *prevcluster, Options &opts, bool st, int K) {
	int dist_cur = 0, dist_prev = 0;
	GenomePos qStart = matches[start].first.pos,  qEnd = qStart + K,
			  tStart = matches[start].second.pos, tEnd = tStart + K;
	for (int i = start + 1; i < end; i++) {
		tEnd = max(tEnd, matches[i].second.pos + K);
		tStart = min(tStart, matches[i].second.pos);
		qEnd = max(qEnd, matches[i].first.pos + K);
		qStart = min(qStart, matches[i].first.pos);
	}
	int qdist = 0, tdist = 0;
	qdist = (qStart >= cluster->qEnd) ? qStart - cluster->qEnd : 0;
	if (st == 0) tdist = (tStart >= cluster->tEnd) ? tStart - cluster->tEnd : 0;
	else tdist = (cluster->tStart >= tEnd) ? cluster->tStart - tEnd : 0;
	dist_cur = max(qdist, tdist);

	qdist = (prevcluster->qStart >= qEnd) ? prevcluster->qStart - qEnd : 0;
	if (st == 0) tdist = (prevcluster->tStart >= tEnd) ? prevcluster->tStart - tEnd : 0;
	else tdist = (tStart >= prevcluster->tEnd) ? tStart - prevcluster->tEnd : 0;
	dist_prev = max(qdist, tdist);
	if (dist_cur <= dist_prev) {
		cluster->matches.insert(cluster->matches.end(), matches.begin() + start, matches.begin() + end); 
		cluster->SetClusterBoundariesFromMatches(opts);
		// cluster->refinespace = 1;			
	}
	else {
		prevcluster->matches.insert(prevcluster->matches.end(), matches.begin() + start, matches.begin() + end); 
		prevcluster->SetClusterBoundariesFromMatches(opts);
		// prevcluster->refinespace = 1;			
	}
}
//
// This function find anchors btwn two adjacent Clusters;
//
int 			
RefineBtwnSpace_AppendCloseCluster (vector<Cluster> &RevBtwnCluster, bool twoblocks, Cluster *cluster, Cluster *prevcluster, Options &opts, Genome &genome, Read &read, char *strands[2], GenomePos qe, GenomePos qs, 
				GenomePos te, GenomePos ts, bool st, GenomePos lrts=0, GenomePos lrlength=0) {

	int ChromIndex = cluster->chromIndex;
	if (st == 1) { 	// If st == 1, then we need to flip this Cluster, since the following code of fining matches requiers that;
		GenomePos t = qs;
		qs = read.length - qe;
		qe = read.length - t;
	}
	//
	// Find matches in read and reference 
	//
	GenomePairs EndPairs;
	RefineSpace(opts.globalK, opts.globalW, 100, 1, EndPairs, opts, genome, read, strands, ChromIndex, qe, qs, te, ts, st, lrts, lrlength);
	float eff = ((float) EndPairs.size()) / min(qe - qs, te - ts);
	// if (twoblocks) cerr << "refineEffiency: " << eff << " original: " << cluster->refineEffiency << endl;
	//
	// If highest pairwise dist between anchors is too large, then divide EndPairs into parts and assign to close clusters; 
	// Else assign the whole EndPairs to close cluster
	//
	if (EndPairs.size() == 0) return 0;
	if (twoblocks and eff >= opts.anchorstoosparse * 2) { // two block happends, just insert if anchors are dense
		cluster->matches.insert(cluster->matches.end(), EndPairs.begin(), EndPairs.end()); 
		cluster->SetClusterBoundariesFromMatches(opts);
		cluster->refinespace = 1;
		return 0;
	}
	if (twoblocks) return 0;
	CartesianSort(EndPairs);
	GenomePos max_pairdist = 0;
	for (int e = 1; e < EndPairs.size(); e++) { max_pairdist = max(max_pairdist, EndPairs[e].first.pos - (EndPairs[e - 1].first.pos + opts.globalK));}
	if (max_pairdist <= 100 and eff >= opts.anchorstoosparse * 2) {
		append_to_closetcluster(EndPairs, 0, EndPairs.size(), cluster, prevcluster, opts, st, opts.globalK);
		// cluster->matches.insert(cluster->matches.end(), EndPairs.begin(), EndPairs.end()); 
		// cluster->SetClusterBoundariesFromMatches(opts);
		// cluster->refinespace = 1;
		return 0;		
	}
	else {
		int start = 0, end = 1;
		while (start < EndPairs.size()) {
			end = start + 1;
			while (end < EndPairs.size() and minGapDifference(EndPairs[end - 1], EndPairs[end]) <= 200) {
				end++;
			}
			if (end - start >= 4) {
				append_to_closetcluster(EndPairs, start, end, cluster, prevcluster, opts, st, opts.globalK);
			}
			start = end;
		}
		return 0;
	}
	//
	// Find inversion
	//
	// bool rst = (st == 1? 0 : 1);
	// GenomePos t = qs;
	// qs = read.length - qe;
	// qe = read.length - t;
	// GenomePairs revEndPairs;
	// RefineSpace(300, 1, revEndPairs, opts, genome, read, strands, ChromIndex, qe, qs, te, ts, rst, lrts, lrlength);		
	// float reff = ((float) revEndPairs.size()) / min(qe - qs, te - ts);
	
	// // cerr << "refineEffiency: " << eff << "  reff: " <<  reff << " original: " << cluster->refineEffiency << endl;
	// if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
	// 	ofstream clust("RevBtwnCluster.tab", std::ofstream::app);
	// 	for (int t = 0; t < revEndPairs.size(); t++) {
	// 			if (rst == 0) {
	// 				clust << revEndPairs[t].first.pos << "\t"
	// 					  << revEndPairs[t].second.pos << "\t"
	// 					  << revEndPairs[t].first.pos + opts.globalK << "\t"
	// 					  << revEndPairs[t].second.pos + opts.globalK << "\t"
	// 					  << rst << endl;
	// 			}
	// 			else {
	// 				clust << revEndPairs[t].first.pos << "\t"
	// 					  << revEndPairs[t].second.pos + opts.globalK << "\t"
	// 					  << revEndPairs[t].first.pos + opts.globalK << "\t"
	// 					  << revEndPairs[t].second.pos<< "\t"
	// 					  << rst << endl;					
	// 			}
	// 	}
	// 	for (int t = 0; t < EndPairs.size(); t++) {
	// 			if (st == 0) {
	// 				clust << EndPairs[t].first.pos << "\t"
	// 					  << EndPairs[t].second.pos << "\t"
	// 					  << EndPairs[t].first.pos + opts.globalK << "\t"
	// 					  << EndPairs[t].second.pos + opts.globalK << "\t"
	// 					  << st << endl;
	// 			}
	// 			else {
	// 				clust << EndPairs[t].first.pos << "\t"
	// 					  << EndPairs[t].second.pos + opts.globalK << "\t"
	// 					  << EndPairs[t].first.pos + opts.globalK << "\t"
	// 					  << EndPairs[t].second.pos<< "\t"
	// 					  << st << endl;					
	// 			}
	// 	}
	// 	clust.close();
	// }	
	// if (eff >= reff) {
	// 	cluster->matches.insert(cluster->matches.end(), EndPairs.begin(), EndPairs.end()); 
	// 	cluster->SetClusterBoundariesFromMatches(opts);
	// 	cluster->refinespace = 1;
	// 	cluster->anchorfreq = 1.0f;
	// 	return 0;
	// }
	// else{
	// 	cerr << "rev happens, refineEffiency: " << eff << "  reff: " <<  reff << " original: " << cluster->refineEffiency << endl;
	// 	cerr << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
	// 	RevBtwnCluster.push_back(Cluster(0, 0, rst));
	// 	RevBtwnCluster.back().matches.insert(RevBtwnCluster.back().matches.end(), revEndPairs.begin(), revEndPairs.end()); 
	// 	RevBtwnCluster.back().SetClusterBoundariesFromMatches(opts);
	// 	RevBtwnCluster.back().refinespace = 1;
	// 	RevBtwnCluster.back().anchorfreq = 1.0f;
	// 	return 1;
	// }
}

//
// This function refines the Clusters in chain and store refined anchors in refinedClusters
// NOTICE: Inside this function, we need to flip reversed Cluster into forward direction to find refined matches;
// And flip them back after the refining step;
//
int 
Refine_splitchain(vector<SplitChain> &splitchains, UltimateChain &chain, vector<Cluster> & refinedclusters, vector<Cluster> &clusters, Genome & genome, Read & read,  
				LocalIndex & glIndex, LocalIndex *localIndexes[2], Options & smallOpts, Options & opts) {
	if (read.unaligned) return 0;
	for (int ph = 0; ph < splitchains.size(); ph++) {
		//
		// Get the boundaries of the cluster in genome sequence.
		//
		if (splitchains[ph].size() == 0) continue;
		refinedclusters[ph].chromIndex = splitchains[ph].chromIndex;
		refinedclusters[ph].refined = 1;
		//
		// Make the anchors reference this chromosome for easier bookkeeping 
		// NOTICE: Remember to add chromOffset back in refinedclusters
		//
		GenomePos chromOffset = genome.header.pos[splitchains[ph].chromIndex];
		for (int c = 0; c < splitchains[ph].ClusterIndex.size(); c++) {
			int cI = splitchains[ph].ClusterIndex[c];
			if (clusters[cI].flip == 0) {
				for (int m = 0; m < clusters[cI].matches.size(); m++) {
					clusters[cI].matches[m].second.pos -= chromOffset; // comment!
				}	
				if (clusters[cI].strand == 1) {
					SwapStrand(read, opts, clusters[cI], opts.globalK);
				}	
				clusters[cI].flip = 1;
			}	
		}

		GenomePos GenomeClusterEnd = splitchains[ph].TEnd;
		GenomePos chromEndOffset = genome.header.GetNextOffset(GenomeClusterEnd);

		int64_t maxDN, minDN;
		maxDN = (int64_t) splitchains[ph].tStart(0) - (int64_t) splitchains[ph].qStart(0); // trans_ takes care of reverse strand and offset
		minDN = maxDN;
		for (int db = 0; db < splitchains[ph].size() - 1; db++) {
			maxDN = max(maxDN, (int64_t) splitchains[ph].tStart(db) - (int64_t) splitchains[ph].qStart(db));
			minDN = min(minDN, (int64_t) splitchains[ph].tStart(db) - (int64_t) splitchains[ph].qStart(db));
		}						
		refinedclusters[ph].maxDiagNum = maxDN + 50; //20
		refinedclusters[ph].minDiagNum = minDN - 50;//20
		//
		// Get shorthand access to alignment boundaries.
		//
		// sorted by second.pos and then first.pos
		GenomePos genomeClusterSegStart, genomeClusterSegEnd;
		// genomeClusterSegStart = splitchains[ph].TStart - chromOffset;
		// genomeClusterSegEnd = splitchains[ph].TEnd - chromOffset;
		genomeClusterSegStart = splitchains[ph].TStart;
		genomeClusterSegEnd = splitchains[ph].TEnd;
		//
		// Search region starts in window, or beginning of chromosome
		//
		int ls, le;
		GenomePos wts = (genomeClusterSegStart >= chromOffset + smallOpts.window)? genomeClusterSegStart - smallOpts.window : chromOffset;
		GenomePos wte = (genomeClusterSegEnd + smallOpts.window < chromEndOffset) ? genomeClusterSegEnd + smallOpts.window : chromEndOffset;

		ls = glIndex.LookupIndex(wts);
		le = glIndex.LookupIndex(wte);
		// 
		// Get quick access to the local index
		//
		LocalIndex *readIndex;
		readIndex = localIndexes[splitchains[ph].Strand];

		for (int lsi = ls; lsi <= le; lsi++) {
			//
			// Find the coordinates in the cluster fragment that start in this local index.
			//
			GenomePos genomeLocalIndexStart = glIndex.seqOffsets[lsi]  - chromOffset;
			GenomePos genomeLocalIndexEnd   = glIndex.seqOffsets[lsi + 1] - 1 - chromOffset;

			// assert(clusters.back().matches.size() == 1); // dummy matches
			// clusters.back().matches[0].first.pos = 0; clusters.back().matches[0].second.pos = genomeLocalIndexStart;
			// int matchStart = splitchains[ph].CartesianTargetLowerBound(0, 1, splitchains[ph].sptc.size() - 1);

			// if (matchStart >= splitchains[ph].size() - 1) continue;
			// clusters.back().matches[0].first.pos = 0; clusters.back().matches[0].second.pos = genomeLocalIndexEnd;
			// int matchEnd = splitchains[ph].CartesianTargetUpperBound(matchStart, 1, splitchains[ph].sptc.size() - 1);

			// // matchEnd += matchStart;
			// assert(matchEnd >= matchStart);
			// if (matchEnd == splitchains[ph].size() - 1) matchEnd--;
			// if (matchStart >= splitchains[ph].size() - 1 or matchStart >= matchEnd) continue; // If there is no overlap with this cluster

			GenomePos prev_readEnd = 0;
			GenomePos prev_readStart = read.length;
			// GenomePos readStart = splitchains[ph].qStart(matchStart);
			// GenomePos readEnd = splitchains[ph].qStart(matchEnd);
			GenomePos readStart = splitchains[ph].qStart(0);
			GenomePos readEnd = splitchains[ph].qEnd(splitchains[ph].size()-1);
			if (readStart == readEnd) { // there is a gap
				if (lsi > ls and readStart > prev_readEnd) {readStart = prev_readEnd;} 
			}
			//
			// Expand boundaries of read to match.
			//
			if (lsi == ls) {readStart = (readStart < smallOpts.window)? 0 : readStart - smallOpts.window;}
			if (lsi == le) { readEnd = (readEnd + smallOpts.window > read.length) ? read.length : readEnd + smallOpts.window;}			
			if (readStart > readEnd) continue; // tandem repear -- get picked up by 3rd SDP;
			//cerr << "readStart: " << readStart << " readEnd: " << readEnd << " ph: " << ph << endl;
			//
			// Find the boundaries where in the query the matches should be added.
			//
			int queryIndexStart = readIndex->LookupIndex(readStart);
			int queryIndexEnd = readIndex->LookupIndex(min(readEnd, (GenomePos)read.length - 1));
			assert(queryIndexEnd < readIndex->seqOffsets.size() + 1);
			GenomePos qStart, qEnd;
			for (int qi = queryIndexStart; qi <= queryIndexEnd; ++qi){ 
				LocalPairs smallMatches;
				GenomePos qStartBoundary = readIndex->tupleBoundaries[qi];
				GenomePos qEndBoundary   = readIndex->tupleBoundaries[qi+1];
				GenomePos readSegmentStart= readIndex->seqOffsets[qi];
				GenomePos readSegmentEnd  = readIndex->seqOffsets[qi+1];

				CompareLists<LocalTuple, SmallTuple>(readIndex->minimizers.begin() + qStartBoundary, readIndex->minimizers.begin() + qEndBoundary, 
									glIndex.minimizers.begin() + glIndex.tupleBoundaries[lsi], 
									glIndex.minimizers.begin()+ glIndex.tupleBoundaries[lsi+1], 
									smallMatches, smallOpts, false, 0, 0, false);

				//
				// Add refined anchors if they fall into the diagonal band and cluster box
				//
				if (splitchains[ph].Strand == 0) {qStart = splitchains[ph].QStart; qEnd = splitchains[ph].QEnd;}
				else {qStart = read.length - splitchains[ph].QEnd; qEnd = read.length - splitchains[ph].QStart;}
				// AppendValues<LocalPairs>(refinedclusters[ph].matches, smallMatches.begin(), smallMatches.end(), readSegmentStart, 
				// 			genomeLocalIndexStart, refinedclusters[ph].maxDiagNum, refinedclusters[ph].minDiagNum, qStart, 
				// 			qEnd, splitchains[ph].TStart - chromOffset, splitchains[ph].TEnd - chromOffset, prev_readStart, prev_readEnd);	
				AppendValues<LocalPairs>(refinedclusters[ph].matches, smallMatches.begin(), smallMatches.end(), readSegmentStart, 
 							genomeLocalIndexStart, refinedclusters[ph].maxDiagNum, refinedclusters[ph].minDiagNum, qStart, 
 							qEnd, splitchains[ph].TStart - chromOffset, splitchains[ph].TEnd - chromOffset, prev_readStart, prev_readEnd);					
			}
		}
		splitchains[ph].clusterIndex = ph;
		if (refinedclusters[ph].matches.size() == 0) continue;
		if (splitchains[ph].Strand == 1) SwapStrand(read, smallOpts, refinedclusters[ph], smallOpts.globalK);
		refinedclusters[ph].SetClusterBoundariesFromMatches(smallOpts);
		refinedclusters[ph].strand = splitchains[ph].Strand;
		refinedclusters[ph].coarse = ph;
		refinedclusters[ph].refinespace = 0;
		refinedclusters[ph].refineEffiency = ((float) refinedclusters[ph].matches.size()) / min(refinedclusters[ph].qEnd - 
											refinedclusters[ph].qStart, refinedclusters[ph].tEnd - refinedclusters[ph].tStart);
		for (int c = 0; c < splitchains[ph].ClusterIndex.size(); c++) {
			int cI = splitchains[ph].ClusterIndex[c];
			if (clusters[cI].flip == 1) {
				for (int m = 0; m < clusters[cI].matches.size(); m++) {
					clusters[cI].matches[m].second.pos += chromOffset; // comment!
				}	
				if (clusters[cI].strand == 1) {
					SwapStrand(read, opts, clusters[cI], opts.globalK);
				}	
				clusters[cI].flip = 0;
			}	
		}
	}
	return 0;
}

void
Refine_Btwnsplitchain(vector<SplitChain> &splitchains, vector<Cluster> &RefinedClusters, vector<Cluster> &RevBtwnCluster, 
						vector<tuple<int, int, int> > &tracerev, Genome &genome, Read &read, Options &opts, char *strands[2], vector<bool> &spchain_link) {
	//
	// Find matches btwn every two adjacent Clusters;
	//
	int c = 1;
	GenomePos qe, qs, te1, ts1;
	GenomePos te2 = 0, ts2 = 0;
	bool st1, st2; 
	bool twoblocks = 0; // twoblocks = 1 when INV happens (need to refine the end of INV)
	GenomePos SpaceLength;
	while (c < splitchains.size()) {
		int cur = c; 
		int prev = c - 1;
		//
		// Decide the boudaries of space and strand direction btwn RefinedClusters[cur] and RefinedClusters[prev]
		//
		if (RefinedClusters[cur].matches.size() == 0 or RefinedClusters[prev].matches.size() == 0) {c++; continue;}
		qs = RefinedClusters[cur].qEnd; qe = RefinedClusters[prev].qStart;
		if (qe <= qs) {c++; continue;}
		if (RefinedClusters[cur].strand == RefinedClusters[prev].strand and spchain_link[c - 1] == 0) { // Missing TRA or INV
			twoblocks = 0;
			st1 = RefinedClusters[cur].strand;
			if (RefinedClusters[cur].tEnd <= RefinedClusters[prev].tStart) {
				ts1 = RefinedClusters[cur].tEnd;
				te1 = RefinedClusters[prev].tStart;
			}
			else if (RefinedClusters[cur].tStart > RefinedClusters[prev].tEnd) {
				ts1 = RefinedClusters[prev].tEnd;
				te1 = RefinedClusters[cur].tStart;
			}
			else {c++; continue;}// No need to refine the space! (DUP)
			ts2 = 0; te2 = 0;
		}
		else if (RefinedClusters[cur].strand != RefinedClusters[prev].strand and spchain_link[c - 1] == 1) { // INV
			st1 = RefinedClusters[cur].strand; 
			st2 = RefinedClusters[prev].strand; 
			twoblocks = 1;
			if (RefinedClusters[cur].tEnd <= RefinedClusters[prev].tStart) {
				if (st1 == 0) {
					ts1 = RefinedClusters[cur].tEnd; te1 = ts1 + qe - qs;
					ts2 = RefinedClusters[prev].tEnd; te2 = ts2 + qe - qs;
				}
				else {
					te1 = RefinedClusters[cur].tStart; ts1 = (te1 > (qe - qs)? te1 - (qe - qs) : 0);
					te2 = RefinedClusters[prev].tStart; ts2 = (te2 > (qe - qs)? te2 - (qe - qs) : 0);
				}
			}
			else if (RefinedClusters[cur].tStart > RefinedClusters[prev].tEnd) {
				if (st1 == 0) {
					ts1 = RefinedClusters[cur].tEnd; te1 = ts1 + qe - qs;
					te2 = RefinedClusters[cur].tStart; ts2 = (te2 > (qe - qs)? te2 - (qe - qs) : 0); 
				}
				else {
					te1 = RefinedClusters[cur].tStart; ts1 = (te1 > (qe - qs)? te1 - (qe - qs) : 0);
					te2 = RefinedClusters[prev].tStart; ts2 = (te2 > (qe - qs)? te2 - (qe - qs) : 0);
				}
			}
			else {c++; continue;}// No need to refine the space!
		}
		else if (RefinedClusters[cur].strand == RefinedClusters[prev].strand and spchain_link[c - 1] == 1) { // DUP
			st1 = RefinedClusters[cur].strand; 
			st2 = st1; 
			twoblocks = 1;
			if (st1 == 0 and RefinedClusters[cur].tEnd > RefinedClusters[prev].tStart) {
				ts1 = RefinedClusters[cur].tEnd; te1 = ts1 + qe - qs;
				te2 = RefinedClusters[prev].tStart; ts2 = (te2 > (qe - qs) ? te2 - (qe - qs) : 0);
			}
			else if (st1 == 1 and RefinedClusters[cur].tStart < RefinedClusters[prev].tEnd) {
				te1 = RefinedClusters[cur].tStart; ts1 = (te1 > (qe - qs) ? te1 - (qe - qs) : 0);
				ts2 = RefinedClusters[prev].tEnd; te2 = ts2 + (qe - qs);
			}
			else {c++; continue;}// No need to refine the space!			
		}

		// //cerr << "btwn  p: " << p << " h: " << h << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
		// if (qe > qs and te1 > ts1) {
		// 	SpaceLength = max(qe - qs, te1 - ts1); 
		// 	//cerr << "SpaceLength: " << SpaceLength << "st: " << st << endl; 
		// 	//used to be 100000; mapping contigs requires larger threshold;
		// 	if (SpaceLength <= 100000 and RefinedClusters[cur]->chromIndex == RefinedClusters[prev]->chromIndex) {
		// 		// btwnClusters have GenomePos, st, matches, coarse
		// 		// This function also set the "coarse" flag for RefinedClusters[cur]
		// 		RefineBtwnSpace(RefinedClusters[cur], smallOpts, genome, read, strands, qe, qs, te1, ts1, st);
		// 	}
		// }
		if (te1 <= ts1) {c++; continue;}
		SpaceLength = min(qe - qs, te1 - ts1); 
		// if (twoblocks) {cerr << "refinetwoblocks " << read.name << " qs: " << qs << " qe: " << qe << " ts1: " << ts1 << " te1: " << te1 << endl;}
		if (SpaceLength <= opts.refineSpaceDist and RefinedClusters[cur].chromIndex == RefinedClusters[prev].chromIndex) {//used to be 100000; mapping contigs requires larger threshold;
			if (RefineBtwnSpace_AppendCloseCluster(RevBtwnCluster, twoblocks, &RefinedClusters[cur], &RefinedClusters[prev], opts, genome, read, strands, qe, qs, te1, ts1, st1)) {
				tracerev.push_back(make_tuple(0, c, RevBtwnCluster.size() - 1)); // Insert a rev cluster
			}
		}		
		if (twoblocks) {
			if (te2 <= ts2) {c++; continue;}
			SpaceLength = min(qe - qs, te2 - ts2); 
			// if (twoblocks) {cerr << "refinetwoblocks " << read.name << " qs: " << qs << " qe: " << qe << " ts1: " << ts1 << " te1: " << te1 << endl;}
			if (SpaceLength <= opts.refineSpaceDist and RefinedClusters[cur].chromIndex == RefinedClusters[prev].chromIndex) {//used to be 100000; mapping contigs requires larger threshold;
				RefineBtwnSpace(opts.globalK, opts.globalW, RevBtwnCluster, twoblocks, &RefinedClusters[prev], opts, genome, read, strands, qe, qs, te2, ts2, st2);
			}				
		}
		c++;
	}
	//
	// Find matches at the right end;
	//
	GenomePos te, ts;
	bool st;
	int rh = splitchains[0].clusterIndex;
	if (RefinedClusters[rh].matches.size() > 0) {
		st = RefinedClusters[rh].strand;
		qs = RefinedClusters[rh].qEnd;
		qe = read.length;
		if (st == 0) {
			ts = RefinedClusters[rh].tEnd;
			te = ts + qe - qs;				
		}
		else {
			te = RefinedClusters[rh].tStart;
			if (te > qe - qs) ts = te - (qe - qs);
			else te = 0;
		}
		// cerr << "rigt: " << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
		if (qe > qs and te > ts) {
			SpaceLength = min(qe - qs, te - ts); 
			if (SpaceLength < opts.refineSpaceDist and te + 500 < genome.lengths[RefinedClusters[rh].chromIndex]) { // used (1000, 6000)
				GenomePos lrts=0, lrlength=0;
				if (st==0) {
					lrts=0;
					lrlength=500;				
				}
				else {
					if (ts>500) lrts=500;
					lrlength=lrts;						
				}
				//				cerr << "right: " << ts-lrts << " length: " << te-ts+lrlength<< endl;
				RefineBtwnSpace(opts.globalK, opts.globalW, RevBtwnCluster, 1, &RefinedClusters[rh], opts, genome, read, strands, qe, qs, te, ts, st, lrts, lrlength);
			}			
		}		
	}

	//
	// Find matches at the left end
	//		
	int lh = splitchains.back().clusterIndex;
	if (RefinedClusters[lh].matches.size() > 0) {
		qs = 0;
		qe = RefinedClusters[lh].qStart;
		st = RefinedClusters[lh].strand;
		if (st == 0) {
			te = RefinedClusters[lh].tStart;
			if (te > qe - qs) ts = te - (qe - qs);
			else ts = 0;
		}
		else {
			ts = RefinedClusters[lh].tEnd;
			te = ts + (qe - qs);
		}
		// cerr << "left qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
		if (qe > qs and te > ts) {
			SpaceLength = min(qe - qs, te - ts);
			if (SpaceLength < opts.refineSpaceDist and te+500 < genome.lengths[RefinedClusters[lh].chromIndex]) { // used (1000, 6000)
				GenomePos lrts=0, lrlength=0;
				if (st==0) { 
					if (ts>500) lrts=500;
					lrlength=lrts;	
				}
				else {
					lrts=0;
					lrlength=500;						
				}
				//cerr << "left: " << ts-lrts << " length: " << te-ts+lrlength<< endl;
				RefineBtwnSpace(opts.globalK, opts.globalW, RevBtwnCluster, 1, &RefinedClusters[lh], opts, genome, read, strands, qe, qs, te, ts, st, lrts, lrlength);
			}			
		}		
	}
}

//
// Merge adjacent cluster when very close on q-coordinates or t-coordinates
// This merge cluster either on a linear chain or INS or DEL
//
void MergeChain(vector<Cluster *> &clusters, vector<Merge_SplitChain> &mergeinfo, SplitChain &merge_sp, SplitChain &sp) {
	vector<int> onec;
	onec.push_back(0);
	GenomePos qe, qs, te, ts;	
	int t = 1, cur = 0, prev = 0;
	int qdist = 0; int tdist = 0;
	while (t < sp.size()) {
		cur = sp[t]; prev = sp[t - 1];
		qdist = 9999; tdist = 9999;
		if (clusters[prev]->strand == clusters[cur]->strand) {
			qdist = (clusters[prev]->qStart > clusters[cur]->qEnd) ? (clusters[prev]->qStart - clusters[cur]->qEnd) : 0;
			if (clusters[prev]->strand == 0) {
				tdist = (clusters[prev]->tStart >= clusters[cur]->tEnd)? clusters[prev]->tStart - clusters[cur]->tEnd : 9999;
			}
			else if (clusters[prev]->strand == 1) {
				tdist = (clusters[prev]->tEnd <= clusters[cur]->tStart)? clusters[cur]->tStart - clusters[prev]->tEnd : 9999;
			}
		}
		if (qdist <= 500 and tdist <= 500) {
			onec.push_back(cur);
		}
		else {
			mergeinfo.push_back(Merge_SplitChain(onec, &clusters));
			onec.clear();
			onec.push_back(cur);
		}
		t++;
	}
	if (!onec.empty()) {
		mergeinfo.push_back(Merge_SplitChain(onec, &clusters));
	}
	merge_sp.sptc.resize(mergeinfo.size());
	for (int t = 0; t < mergeinfo.size(); t++) {
		merge_sp.sptc[t] = t;
	}
}

#endif