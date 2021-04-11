#ifndef CLUSTER_REFINE_H_
#define CLUSTER_REFINE_H_
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

#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cmath>	// std::log 
#include <sstream>
#include <thread>
#include <climits>

void SwapStrand (Read & read, const Options & opts, Cluster & cluster, int K) {
	for (int m = 0; m < cluster.matches.size(); m++) {
		cluster.matches[m].first.pos = read.length - (cluster.matches[m].first.pos + K);
	}
	GenomePos r = cluster.qStart;
	cluster.qStart = read.length - cluster.qEnd;
	cluster.qEnd = read.length - r;
}

void SwapStrand(Read &read, const Options &opts, GenomePairs &matches, int K) {
	for (int m=0; m < matches.size(); m++) {
		matches[m].first.pos = read.length - (matches[m].first.pos + K);
	}
}

void SwapStrand(Read &read, const Options &opts, GenomePairs &matches, int start, int end, int K) {
	for (int m=start; m<end; m++) {
		matches[m].first.pos = read.length - (matches[m].first.pos + K);
	}
}

//
// This function refines the Clusters in chain and store refined anchors in refinedClusters
// NOTICE: Inside this function, we need to flip reversed Cluster into forward direction to find refined matches;
// And flip them back after the refining step;
//
int 
REFINEclusters(vector<Cluster> & clusters, vector<Cluster> & refinedclusters, Genome & genome, Read & read,  
				LocalIndex & glIndex, LocalIndex *localIndexes[2], const Options & smallOpts, const Options & opts) {
	if (read.unaligned) return 0;
	for (int ph = 0; ph < clusters.size(); ph++) {
		//
		// Get the boundaries of the cluster in genome sequence.
		//
		if (clusters[ph].matches.size() == 0) continue;
		if (clusters[ph].refined == 1) continue; // this Cluster has been refined;

		int pass = clusters[ph].CHROMIndex(genome);
		if (pass == 1) {
			clusters[ph].matches.clear();
			continue;			
		}
		refinedclusters[ph].chromIndex = clusters[ph].chromIndex;
		clusters[ph].refined = 1;
		//
		// Make the anchors reference this chromosome for easier bookkeeping 
		// NOTICE: Remember to add chromOffset back in refinedclusters
		//
		GenomePos chromOffset = genome.header.pos[clusters[ph].chromIndex];
		for (int m = 0; m < clusters[ph].matches.size(); m++) {
			clusters[ph].matches[m].second.pos -= chromOffset;
		}
		GenomePos GenomeClusterEnd = clusters[ph].tEnd;
		GenomePos chromEndOffset = genome.header.GetNextOffset(GenomeClusterEnd);
		//
		// If the current Cluster is reversed stranded, swap the anchors and GenomePos to forward strand; 
		// This is for finding refined anchors;
		// NOTICE: Need to flip such Cluster back into reversed strand;
		//
		if (clusters[ph].strand == 1) SwapStrand(read, opts, clusters[ph], opts.globalK);
		// 
		// Decide the diagonal band for each clusters[ph]
		// Find the digonal band that each clusters[ph] is in; 
		// NOTICE: here every diagnoal have already subtracted chromOffset, so it's in the same scale with local matches
		// 
		int64_t maxDN, minDN;
		maxDN = (int64_t) clusters[ph].matches[0].second.pos - (int64_t) clusters[ph].matches[0].first.pos;
		minDN = maxDN;
		for (int db = 0; db < clusters[ph].matches.size(); db++) {
			maxDN = max(maxDN, (int64_t)clusters[ph].matches[db].second.pos - (int64_t)clusters[ph].matches[db].first.pos);
			minDN = min(minDN, (int64_t)clusters[ph].matches[db].second.pos - (int64_t)clusters[ph].matches[db].first.pos);
		}						
		clusters[ph].maxDiagNum = maxDN + (int64_t)100; //20
		clusters[ph].minDiagNum = minDN - (int64_t)100;//20
		//
		// Get shorthand access to alignment boundaries.
		//
		// sorted by second.pos and then first.pos
		CartesianTargetSort<GenomeTuple>(clusters[ph].matches.begin(), clusters[ph].matches.end()); 
		GenomePos genomeClusterSegStart, genomeClusterSegEnd;
		genomeClusterSegStart = clusters[ph].tStart;
		genomeClusterSegEnd = clusters[ph].tEnd;
		//
		// Search region starts in window, or beginning of chromosome
		//
		int ls, le;
		GenomePos wts, wte;
		if (chromOffset + smallOpts.window > genomeClusterSegStart) {
			wts = chromOffset;
		}
		else {
			wts = genomeClusterSegStart - smallOpts.window;
		}
				
		if (genomeClusterSegEnd + smallOpts.window > chromEndOffset) {
			wte = chromEndOffset-1;
		}
		else {
			wte = genomeClusterSegEnd + smallOpts.window;
		}
			
		ls = glIndex.LookupIndex(wts);
		le = glIndex.LookupIndex(wte);
		// 
		// Get quick access to the local index
		//
		LocalIndex *readIndex;
		readIndex = localIndexes[clusters[ph].strand];

		for (int lsi = ls; lsi <= le; lsi++) {
			//
			// Find the coordinates in the cluster fragment that start in this local index.
			//
			if (glIndex.seqOffsets[lsi] < chromOffset or glIndex.seqOffsets[lsi + 1] < chromOffset) continue; 
			GenomePos genomeLocalIndexStart = glIndex.seqOffsets[lsi]  - chromOffset;
			GenomePos genomeLocalIndexEnd = glIndex.seqOffsets[lsi+1] - 1 - chromOffset;
			if (genomeLocalIndexStart >= genomeLocalIndexEnd) continue;

			int matchStart = CartesianTargetLowerBound<GenomeTuple>(clusters[ph].matches.begin(), clusters[ph].matches.end(), genomeLocalIndexStart);

			int matchEnd = CartesianTargetUpperBound<GenomeTuple>(clusters[ph].matches.begin()+matchStart, clusters[ph].matches.end(), genomeLocalIndexEnd);
			matchEnd += matchStart;
			//assert(matchEnd > matchStart);
			if (matchEnd == clusters[ph].matches.size()) matchEnd--;
			//
			// If there is no overlap with this cluster
			//
			if (matchStart >= clusters[ph].matches.size()) {
				continue;
			}
			GenomePos prev_readEnd = 0;
			GenomePos prev_readStart = read.length;
			GenomePos readStart = clusters[ph].matches[matchStart].first.pos;
			GenomePos readEnd=clusters[ph].matches[matchEnd].first.pos;
			// long startDiag=clusters[ph].matches[matchStart].second.pos - clusters[ph].matches[matchStart].first.pos;
			// long endDiag=clusters[ph].matches[matchEnd].second.pos - clusters[ph].matches[matchEnd].first.pos;
			// long minDiag=min(startDiag, endDiag);
			// long maxDiag=max(startDiag, endDiag);
			if (readStart == readEnd) { // there is a gap
				if (lsi > ls and readStart > prev_readEnd) {
					readStart = prev_readEnd;
				} 
			}
			//
			// Expand boundaries of read to match.
			//
			if (lsi == ls) {
				if (readStart < smallOpts.window) {
					readStart = 0;
				}
				else {
					readStart -= smallOpts.window;
				}
			}
			if (lsi == le) {
				if (readEnd + smallOpts.window > read.length) {
					readEnd = read.length; 
				}
				else { 
					readEnd += smallOpts.window;	
				}
			}			
			if (readStart > readEnd) continue; // tandem repear -- get picked up by 3rd SDP;
			//cerr << "readStart: " << readStart << " readEnd: " << readEnd << " ph: " << ph << endl;
			//
			// when readStart > readEnd, it means that the refine strategy does not work well due to a tandem repeat;
			// Re-find readStart and readEnd based on maxdiag and mindiag
			//
			// if (readStart > readEnd) {
			// 	if (genomeLocalIndexStart <= clusters[ph].maxDiagNum) readStart = 0;
			// 	else readStart = genomeLocalIndexStart - clusters[ph].maxDiagNum;
			// 	if (genomeLocalIndexEnd <= clusters[ph].minDiagNum) readEnd = 0;
			// 	else readEnd = genomeLocalIndexEnd - clusters[ph].minDiagNum;
			// }

			//
			// Find the boundaries where in the query the matches should be added.
			//
			int queryIndexStart = readIndex->LookupIndex(readStart);
			int queryIndexEnd = readIndex->LookupIndex(min(readEnd, (GenomePos) read.length-1));
			assert(queryIndexEnd < readIndex->seqOffsets.size()+1);
			//			cout << "ph\t"<< ph << "\tlsi " << lsi << "\tle" << le << "\tqis " << queryIndexStart << "\t" << queryIndexEnd << "\t" << readIndex->tupleBoundaries.size() << endl;
			for (int qi = queryIndexStart; qi <= queryIndexEnd; ++qi){ 
				LocalPairs smallMatches;
				GenomePos qStartBoundary = readIndex->tupleBoundaries[qi];
				GenomePos qEndBoundary   = readIndex->tupleBoundaries[qi+1];
				GenomePos readSegmentStart= readIndex->seqOffsets[qi];
				GenomePos readSegmentEnd  = readIndex->seqOffsets[qi+1];

				CompareLists<LocalTuple, SmallTuple>(readIndex->minimizers.begin()+qStartBoundary, 
									readIndex->minimizers.begin()+qEndBoundary, 
									glIndex.minimizers.begin()+ glIndex.tupleBoundaries[lsi], 
									glIndex.minimizers.begin()+ glIndex.tupleBoundaries[lsi+1], smallMatches, smallOpts, false, 0, 0, false);
				//
				// Add refined anchors if they fall into the diagonal band and cluster box

				//				AppendValues<LocalPairs>(refinedclusters[ph].matches, smallMatches.begin(), smallMatches.end(), readSegmentStart, genomeLocalIndexStart);

				AppendValues<LocalPairs>(refinedclusters[ph].matches, smallMatches.begin(), smallMatches.end(), readSegmentStart, 
										genomeLocalIndexStart, clusters[ph].maxDiagNum, clusters[ph].minDiagNum, clusters[ph].qStart, 
										clusters[ph].qEnd, clusters[ph].tStart-chromOffset, clusters[ph].tEnd-chromOffset, 
										prev_readStart, prev_readEnd);

			
			}
		}
		if (refinedclusters[ph].matches.size() == 0) continue;
		if (clusters[ph].strand == 1) SwapStrand(read, smallOpts, refinedclusters[ph], smallOpts.globalK);
		refinedclusters[ph].SetClusterBoundariesFromMatches(smallOpts);
		refinedclusters[ph].strand = clusters[ph].strand;
		refinedclusters[ph].coarse = -1;
		refinedclusters[ph].refinespace = 0;
		refinedclusters[ph].refineEffiency = ((float) refinedclusters[ph].matches.size()) / min(refinedclusters[ph].qEnd - refinedclusters[ph].qStart, 
																							refinedclusters[ph].tEnd - refinedclusters[ph].tStart);
	}
	return 0;
}

int RefineSpace(int K, int W, int refineSpaceDiag, bool consider_str, GenomePairs &EndPairs, const Options & opts, Genome & genome, Read & read, char *strands[2], int &ChromIndex, GenomePos qe, 
				GenomePos qs, GenomePos te, GenomePos ts, bool st, GenomePos lrts=0, GenomePos lrlength=0) {
	//
	// Decide the diagonal band for this space
	//
	int64_t minDiagNum, maxDiagNum; 
	int64_t diag1, diag2;
	diag1 = 0;
	diag2 = (int64_t) (te - (ts - lrts)) - (int64_t) (qe - qs); // scale diag1 and diag2 to the local coordinates
	minDiagNum = min(diag1, diag2) - refineSpaceDiag; 
	maxDiagNum = max(diag1, diag2) + refineSpaceDiag; 

	vector<GenomeTuple> EndReadTup, EndGenomeTup;

	string refSeq(genome.seqs[ChromIndex]+(ts-lrts), te-ts+lrlength);
	string querySeq(strands[st] + qs, qe - qs);
	  
	StoreMinimizers_noncanonical<GenomeTuple, Tuple>(genome.seqs[ChromIndex]+(ts-lrts), te-ts+lrlength, K, W, EndGenomeTup, false); // local minimizer
	sort(EndGenomeTup.begin(), EndGenomeTup.end());
	StoreMinimizers_noncanonical<GenomeTuple, Tuple>(strands[st] + qs, qe - qs, K, W, EndReadTup, false);
	sort(EndReadTup.begin(), EndReadTup.end());
	CompareLists<GenomeTuple, Tuple>(EndReadTup.begin(), EndReadTup.end(), EndGenomeTup.begin(), EndGenomeTup.end(), EndPairs, opts, false, maxDiagNum, minDiagNum, false); // By passing maxDiagNum and minDiagNum, this function
										// filters out anchors that are outside the diagonal band;

	for (int rm = 0; rm < EndPairs.size(); rm++) {
		EndPairs[rm].first.pos += qs;
		EndPairs[rm].second.pos += ts-lrts;
		assert(EndPairs[rm].first.pos + K <= read.length);
		assert(EndPairs[rm].second.pos + K <= genome.lengths[ChromIndex]);
		if (consider_str == true and st == 1) EndPairs[rm].first.pos = read.length - EndPairs[rm].first.pos - K;
		assert(EndPairs[rm].first.pos + K <= read.length);
	}	
	return 0;
}

//
// This function find anchors btwn two adjacent Clusters;
//
int 			
RefineBtwnSpace(int K, int W, vector<Cluster> &RevBtwnCluster, bool twoblocks, Cluster *cluster, const Options &opts, Genome &genome, Read &read, char *strands[2], GenomePos qe, GenomePos qs, 
				GenomePos te, GenomePos ts, bool st, GenomePos lrts=0, GenomePos lrlength=0) {

	int ChromIndex = cluster->chromIndex;
	if (st == 1) { 	// If st == 1, then we need to flip this Cluster, since the following code of fining matches requiers that;
		GenomePos t = qs;
		qs = read.length - qe;
		qe = read.length - t;
	}
	int refineSpaceDiag = 0;
	if (opts.readType == Options::contig or opts.readType == Options::ccs ) {
	  refineSpaceDiag = min((int) floor(max(100.f,0.01f * (qe - qs))), 100);
	}
	else if (opts.readType == Options::ccs) {
		refineSpaceDiag = min((int) floor(0.02f * (qe - qs)), 300);
	}
	else if (opts.readType == Options::clr or opts.readType == Options::ont) {
	  refineSpaceDiag = min((int) floor(max(100.f,0.15f * (qe - qs))), 1000);	
	}

	// cerr << refineSpaceDiag << endl;
	//
	// Find matches in read and reference 
	//
	GenomePairs EndPairs;
	RefineSpace(K, W, refineSpaceDiag, 1, EndPairs, opts, genome, read, strands, ChromIndex, qe, qs, te, ts, st, lrts, lrlength);
	float eff = ((float) EndPairs.size()) / min(qe - qs, te - ts);
	// if (twoblocks) cerr << "refineEffiency: " << eff << " original: " << cluster->refineEffiency << endl;
	if ((EndPairs.size() > 0 and twoblocks) or (EndPairs.size() > 0 and eff >= opts.anchorstoosparse * 2)) {
		cluster->matches.insert(cluster->matches.end(), EndPairs.begin(), EndPairs.end()); 
		cluster->SetClusterBoundariesFromMatches(opts);
		cluster->refinespace = 1;
		return 0;
	}
	// cerr << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;

	if (twoblocks) return 0;
	bool rst = (st == 1? 0 : 1);
	GenomePos t = qs;
	qs = read.length - qe;
	qe = read.length - t;
	GenomePairs revEndPairs;
	RefineSpace(K, W, refineSpaceDiag, 1, revEndPairs, opts, genome, read, strands, ChromIndex, qe, qs, te, ts, rst, lrts, lrlength);		
	float reff = ((float) revEndPairs.size()) / min(qe - qs, te - ts);
	
	// cerr << "refineEffiency: " << eff << "  reff: " <<  reff << " original: " << cluster->refineEffiency << endl;
	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
		ofstream clust("RevBtwnCluster.tab", std::ofstream::app);
		for (int t = 0; t < revEndPairs.size(); t++) {
				if (rst == 0) {
					clust << revEndPairs[t].first.pos << "\t"
						  << revEndPairs[t].second.pos << "\t"
						  << revEndPairs[t].first.pos + K << "\t"
						  << revEndPairs[t].second.pos + K << "\t"
						  << rst << endl;
				}
				else {
					clust << revEndPairs[t].first.pos << "\t"
						  << revEndPairs[t].second.pos + K << "\t"
						  << revEndPairs[t].first.pos + K << "\t"
						  << revEndPairs[t].second.pos<< "\t"
						  << rst << endl;					
				}
		}
		for (int t = 0; t < EndPairs.size(); t++) {
				if (st == 0) {
					clust << EndPairs[t].first.pos << "\t"
						  << EndPairs[t].second.pos << "\t"
						  << EndPairs[t].first.pos + K << "\t"
						  << EndPairs[t].second.pos + K << "\t"
						  << st << endl;
				}
				else {
					clust << EndPairs[t].first.pos << "\t"
						  << EndPairs[t].second.pos + K << "\t"
						  << EndPairs[t].first.pos + K << "\t"
						  << EndPairs[t].second.pos<< "\t"
						  << st << endl;					
				}
		}
		clust.close();
	}	
	if (eff >= reff) {
		cluster->matches.insert(cluster->matches.end(), EndPairs.begin(), EndPairs.end()); 
		cluster->SetClusterBoundariesFromMatches(opts);
		cluster->refinespace = 1;
		cluster->anchorfreq = 1.0f;
		return 0;
	}
	else{
		//		cerr << "rev happens, refineEffiency: " << eff << "  reff: " <<  reff << " original: " << cluster->refineEffiency << endl;
		//		cerr << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
		RevBtwnCluster.push_back(Cluster(0, 0, rst));
		RevBtwnCluster.back().matches.insert(RevBtwnCluster.back().matches.end(), revEndPairs.begin(), revEndPairs.end()); 
		RevBtwnCluster.back().SetClusterBoundariesFromMatches(opts);
		RevBtwnCluster.back().refinespace = 1;
		RevBtwnCluster.back().anchorfreq = 1.0f;
		return 1;
	}
}

void
RefineBtwnClusters_chain(int K, int W, vector<Primary_chain> &Primary_chains, vector<Cluster*> &RefinedClusters, vector<Cluster> &RevBtwnCluster, 
						vector<tuple<int, int, int> > &tracerev, Genome &genome, Read &read, const Options &smallOpts, int &p, int &h, char *strands[2]) {
	//
	// Find matches btwn every two adjacent Clusters;
	//
	int c = 1;
	GenomePos qe, qs, te1, ts1;
	GenomePos te2 = 0, ts2 = 0;
	bool st1, st2; 
	bool twoblocks = 0;
	int SpaceLength;
	int low_b = 20;
	int SpaceLength_upper;
	if (smallOpts.readType == Options::contig) {
		low_b = 1000;
		SpaceLength_upper = 100000;
	} 
	else {
		SpaceLength_upper = 50000;
	}
	while (c < Primary_chains[p].chains[h].ch.size()) {

		int cur = Primary_chains[p].chains[h].ch[c];
		int prev = Primary_chains[p].chains[h].ch[c - 1];
		//
		// Decide the boudaries of space and strand direction btwn RefinedClusters[cur] and RefinedClusters[prev]
		//
		qs = RefinedClusters[cur]->qEnd; 
		qe = RefinedClusters[prev]->qStart;
		te1 = 0; ts1 = 0; te2 = 0; ts2 = 0;
		if (qe <= qs) {c++; continue;}
		if (smallOpts.readType == Options::contig) twoblocks = 0; // Do not refine end for INV and DUP when aligning contig
		if (RefinedClusters[cur]->strand == RefinedClusters[prev]->strand) {
			twoblocks = 0;
			st1 = RefinedClusters[cur]->strand;
			if (RefinedClusters[cur]->tEnd <= RefinedClusters[prev]->tStart) {
				ts1 = RefinedClusters[cur]->tEnd;
				te1 = RefinedClusters[prev]->tStart;
			}
			else if (RefinedClusters[cur]->tStart > RefinedClusters[prev]->tEnd) {
				ts1 = RefinedClusters[prev]->tEnd;
				te1 = RefinedClusters[cur]->tStart;
			}
			else {c++; continue;}// No need to refine the space!
			ts2 = 0; te2 = 0;
		}
		else if (smallOpts.readType != Options::contig) {
			st1 = RefinedClusters[cur]->strand; 
			st2 = RefinedClusters[prev]->strand; 
			twoblocks = 1;
			if (RefinedClusters[cur]->tEnd <= RefinedClusters[prev]->tStart) {
				if (st1 == 0) {
					ts1 = RefinedClusters[cur]->tEnd; te1 = ts1 + qe - qs;
					ts2 = RefinedClusters[prev]->tEnd; te2 = ts2 + qe - qs;
				}
				else {
					te1 = RefinedClusters[cur]->tStart; ts1 = (te1 > (qe - qs)? te1 - (qe - qs) : 0);
					te2 = RefinedClusters[prev]->tStart; ts2 = (te2 > (qe - qs)? te2 - (qe - qs) : 0);
				}
			}
			else if (RefinedClusters[cur]->tStart > RefinedClusters[prev]->tEnd) {
				if (st1 == 0) {
					ts1 = RefinedClusters[cur]->tEnd; te1 = ts1 + qe - qs;
					te2 = RefinedClusters[cur]->tStart; ts2 = (te2 > (qe - qs)? te2 - (qe - qs) : 0); 
				}
				else {
					te1 = RefinedClusters[cur]->tStart; ts1 = (te1 > (qe - qs)? te1 - (qe - qs) : 0);
					te2 = RefinedClusters[prev]->tStart; ts2 = (te2 > (qe - qs)? te2 - (qe - qs) : 0);
				}
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
		SpaceLength = max(qe - qs, te1 - ts1); 
		//		if (twoblocks) cerr << read.name << " qs: " << qs << " qe: " << qe << " ts1: " << ts1 << " te1: " << te1 << endl;
		if (SpaceLength >= low_b and SpaceLength <= SpaceLength_upper and RefinedClusters[cur]->chromIndex == RefinedClusters[prev]->chromIndex) {//used to be 100000; mapping contigs requires larger threshold;
			// btwnClusters have GenomePos, st, matches, coarse
			// This function also set the "coarse" flag for RefinedClusters[cur]
			if (RefineBtwnSpace(K, W, RevBtwnCluster, twoblocks, RefinedClusters[cur], smallOpts, genome, read, strands, qe, qs, te1, ts1, st1)) {
				// 
				// Insert a rev cluster
				//
				tracerev.push_back(make_tuple(h, c, RevBtwnCluster.size() - 1));
			}
		}		
		if (te2 <= ts2) {c++; continue;}
		SpaceLength = max(qe - qs, te2 - ts2); 
		//		if (twoblocks) cerr << read.name << " qs: " << qs << " qe: " << qe << " ts2: " << ts2 << " te2: " << te2 << endl;
		if (SpaceLength >= low_b and SpaceLength <= SpaceLength_upper and RefinedClusters[cur]->chromIndex == RefinedClusters[prev]->chromIndex) {//used to be 100000; mapping contigs requires larger threshold;
			if (smallOpts.readType == Options::contig) {
				cerr << "Shouldn't run this " << read.name << endl;
			}
			// btwnClusters have GenomePos, st, matches, coarse
			// This function also set the "coarse" flag for RefinedClusters[cur]
			RefineBtwnSpace(K, W, RevBtwnCluster, twoblocks, RefinedClusters[prev], smallOpts, genome, read, strands, qe, qs, te2, ts2, st2);
		}			
		c++;
	}
	//
	// Find matches at the right end;
	//
	GenomePos te, ts;
	bool st;
	int rh = Primary_chains[p].chains[h].ch[0];
	st = RefinedClusters[rh]->strand;
	qs = RefinedClusters[rh]->qEnd;
	qe = read.length;
	if (st == 0) {
		ts = RefinedClusters[rh]->tEnd;
		te = ts + qe - qs;				
	}
	else {
		te = RefinedClusters[rh]->tStart;
		if (te > qe - qs) ts = te - (qe - qs);
		else te = 0;
	}
	if (qe > qs and te > ts) {
		SpaceLength = max(qe - qs, te - ts); 
		if (SpaceLength >= low_b and SpaceLength < SpaceLength_upper and te+500 < genome.lengths[RefinedClusters[rh]->chromIndex]) { // used (1000, 6000)
			GenomePos lrts=0, lrlength=0;
			if (st==0) {
				lrts=0;
				lrlength=500;				
			}
			else {
				if (ts>500) lrts=500;
				lrlength=lrts;						
			}
			//cerr << "right: " << ts-lrts << " length: " << te-ts+lrlength<< endl;
			RefineBtwnSpace(K, W, RevBtwnCluster, 1, RefinedClusters[rh], smallOpts, genome, read, strands, qe, qs, te, ts, st, lrts, lrlength);
		}			
	}
	//
	// Find matches at the left end
	//		
	int lh = Primary_chains[p].chains[h].ch.back();
	qs = 0;
	qe = RefinedClusters[lh]->qStart;
	st = RefinedClusters[lh]->strand;
	if (st == 0) {
		te = RefinedClusters[lh]->tStart;
		if (te > qe - qs) ts = te - (qe - qs);
		else ts = 0;
	}
	else {
		ts = RefinedClusters[lh]->tEnd;
		te = ts + (qe - qs);
	}
	// cerr << "left  p: " << p << " h: " << h << "chrom: " << RefinedClusters[rh]->chromIndex << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
	if (qe > qs and te > ts) {
		SpaceLength = max(qe - qs, te - ts);
		if (SpaceLength >= low_b and SpaceLength < SpaceLength_upper and te+500 < genome.lengths[RefinedClusters[lh]->chromIndex]) { // used (1000, 6000)
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
			RefineBtwnSpace(K, W, RevBtwnCluster, 1, RefinedClusters[lh], smallOpts, genome, read, strands, qe, qs, te, ts, st, lrts, lrlength);
		}			
	}
	// timing.Tick("refine_btwnclusters");
}
#endif
