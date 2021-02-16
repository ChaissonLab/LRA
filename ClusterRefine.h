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
#include "overload.h"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cmath>	// std::log 
#include <sstream>
#include <thread>
#include <climits>

void SwapStrand (Read & read, Options & opts, Cluster & cluster) {
	for (int m = 0; m < cluster.matches.size(); m++) {
		cluster.matches[m].first.pos = read.length - (cluster.matches[m].first.pos + opts.globalK);
	}
	GenomePos r = cluster.qStart;
	cluster.qStart = read.length - cluster.qEnd;
	cluster.qEnd = read.length - r;
}

void SwapStrand(Read &read, Options &opts, GenomePairs &matches) {
	for (int m=0; m < matches.size(); m++) {
		matches[m].first.pos = read.length - (matches[m].first.pos + opts.globalK);
	}
}

void SwapStrand(Read &read, Options &opts, GenomePairs &matches, int start, int end) {
	for (int m=start; m<end; m++) {
		matches[m].first.pos = read.length - (matches[m].first.pos + opts.globalK);
	}
}

//
// This function decide the chromIndex
//
int
CHROMIndex(Cluster & cluster, Genome & genome) {

	if (cluster.matches.size() == 0) return 1;
	GenomePos tPos = cluster.tStart;
	int firstChromIndex = genome.header.Find(tPos);
	int lastChromIndex;
	tPos = cluster.tEnd;
	lastChromIndex = genome.header.Find(tPos);
	if (firstChromIndex != lastChromIndex ) {
		return 1;
	}
	cluster.chromIndex = firstChromIndex;  
	return 0;
}

//
// This function refines the Clusters in chain and store refined anchors in refinedClusters
// NOTICE: Inside this function, we need to flip reversed Cluster into forward direction to find refined matches;
// And flip them back after the refining step;
//
int 
REFINEclusters(vector<Cluster> & clusters, vector<Cluster> & refinedclusters, Genome & genome, Read & read,  
				LocalIndex & glIndex, LocalIndex *localIndexes[2], Options & smallOpts, Options & opts) {

	for (int ph = 0; ph < clusters.size(); ph++) {
		//
		// Get the boundaries of the cluster in genome sequence.
		//
		if (clusters[ph].matches.size() == 0) continue;
		if (clusters[ph].refined == 1) continue; // this Cluster has been refined;

		int pass = CHROMIndex(clusters[ph], genome);
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
		if (clusters[ph].strand == 1) SwapStrand(read, opts, clusters[ph]);
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
		clusters[ph].maxDiagNum = maxDN + 20; //20
		clusters[ph].minDiagNum = minDN - 20;//20
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
			GenomePos genomeLocalIndexStart = glIndex.seqOffsets[lsi]  - chromOffset;
			GenomePos genomeLocalIndexEnd   = glIndex.seqOffsets[lsi+1] - 1 - chromOffset;

			int matchStart = CartesianTargetLowerBound<GenomeTuple>(clusters[ph].matches.begin(), clusters[ph].matches.end(), 
																		genomeLocalIndexStart);

			int matchEnd = CartesianTargetUpperBound<GenomeTuple>(clusters[ph].matches.begin()+matchStart, clusters[ph].matches.end(), 
																		genomeLocalIndexEnd);
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
				//
				AppendValues<LocalPairs>(refinedclusters[ph].matches, smallMatches.begin(), smallMatches.end(), readSegmentStart, 
										genomeLocalIndexStart,clusters[ph].maxDiagNum, clusters[ph].minDiagNum, clusters[ph].qStart, 
										clusters[ph].qEnd, clusters[ph].tStart-chromOffset, clusters[ph].tEnd-chromOffset, 
										prev_readStart, prev_readEnd);
			}
		}
		if (refinedclusters[ph].matches.size() == 0) continue;
		if (clusters[ph].strand == 1) SwapStrand(read, smallOpts, refinedclusters[ph]);
		refinedclusters[ph].SetClusterBoundariesFromMatches(smallOpts);
		refinedclusters[ph].strand = clusters[ph].strand;
		//refinedclusters[ph].minDiagNum = clusters[ph].minDiagNum;
		//refinedclusters[ph].maxDiagNum = clusters[ph].maxDiagNum;
		refinedclusters[ph].coarse = -1;
		refinedclusters[ph].refinespace = 0;
	}
	return 0;
}

int RefineSpace(bool consider_str, GenomePairs &EndPairs, Options & opts, Genome & genome, Read & read, char *strands[2], int &ChromIndex, GenomePos qe, 
				GenomePos qs, GenomePos te, GenomePos ts, bool st, GenomePos lrts=0, GenomePos lrlength=0) {
	//
	// Decide the diagonal band for this space
	//
	int64_t minDiagNum, maxDiagNum; 
	int64_t diag1, diag2;
	diag1 = 0;
	diag2 = (int64_t) (te - (ts - lrts)) - (int64_t) (qe - qs); // scale diag1 and diag2 to the local coordinates
	minDiagNum = min(diag1, diag2) - opts.refineSpaceDiag; 
	maxDiagNum = max(diag1, diag2) + opts.refineSpaceDiag; 

	vector<GenomeTuple> EndReadTup, EndGenomeTup;
	StoreMinimizers<GenomeTuple, Tuple>(genome.seqs[ChromIndex]+(ts-lrts), te-ts+lrlength, opts.globalK, opts.globalW, 
											EndGenomeTup, true, false);
	sort(EndGenomeTup.begin(), EndGenomeTup.end());
	StoreMinimizers<GenomeTuple, Tuple>(strands[st] + qs, qe - qs, opts.globalK, opts.globalW, EndReadTup, true, false);
	sort(EndReadTup.begin(), EndReadTup.end());
	CompareLists<GenomeTuple, Tuple>(EndReadTup.begin(), EndReadTup.end(), EndGenomeTup.begin(), EndGenomeTup.end(), EndPairs, 
										opts, true, maxDiagNum, minDiagNum, false); // By passing maxDiagNum and minDiagNum, this function
										// filters out anchors that are outside the diagonal band;

	for (int rm = 0; rm < EndPairs.size(); rm++) {
		EndPairs[rm].first.pos += qs;
		EndPairs[rm].second.pos += ts-lrts;
		assert(EndPairs[rm].first.pos < read.length);
		if (consider_str == true and st == 1) EndPairs[rm].first.pos = read.length - EndPairs[rm].first.pos - opts.globalK;
		
	}	
	return 0;
}

//
// This function find anchors btwn two adjacent Clusters;
//
int 			
RefineBtwnSpace(Cluster *cluster, Options &opts, Genome &genome, Read &read, char *strands[2], GenomePos qe, GenomePos qs, 
							GenomePos te, GenomePos ts, bool st, GenomePos lrts=0, GenomePos lrlength=0) {

	int ChromIndex = cluster->chromIndex;
	// 
	// If st == 1, then we need to flip this Cluster, since the following code of fining matches requiers that;
	//
	if (st == 1) {
		GenomePos t = qs;
		qs = read.length - qe;
		qe = read.length - t;
	}
	//
	// Find matches in read and reference 
	//
	GenomePairs EndPairs;
	opts.refineSpaceDiag = 300;
	RefineSpace(1, EndPairs, opts, genome, read, strands, ChromIndex, qe, qs, te, ts, st, lrts, lrlength);

	if (EndPairs.size() > 0) {
		cluster->matches.insert(cluster->matches.end(), EndPairs.begin(), EndPairs.end());  // TODO(Jingwen): Timing consuming???????
		cluster->SetClusterBoundariesFromMatches(opts);
		cluster->refinespace = 1;
	}
	return 0;
}

void
RefineBtwenClusters_chain(vector<Primary_chain> &Primary_chains, vector<SegAlignmentGroup> &alignments, vector<Cluster*> &RefinedClusters, Genome &genome, 
	Read &read, Options &smallOpts, int &p, int &h, char *strands[2]) {
	// for (int p = 0; p < Primary_chains.size(); p++) {
	
	// 	for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			
			//
			// Find matches btwn every two adjacent Clusters;
			//
			int c = 1;
			GenomePos qe, qs, te, ts;
			bool st; GenomePos SpaceLength;
			while (c < Primary_chains[p].chains[h].ch.size()) {

				int cur = Primary_chains[p].chains[h].ch[c];
				int prev = Primary_chains[p].chains[h].ch[c - 1];
				//
				// Decide the boudaries of space and strand direction btwn RefinedClusters[cur] and RefinedClusters[prev]
				//
				if (RefinedClusters[cur]->strand == RefinedClusters[prev]->strand) st = RefinedClusters[cur]->strand;
				else st = RefinedClusters[cur]->strand;
				qs = RefinedClusters[cur]->qEnd; 
				qe = RefinedClusters[prev]->qStart;

				if (RefinedClusters[cur]->tEnd <= RefinedClusters[prev]->tStart) {
					ts = RefinedClusters[cur]->tEnd;
					te = RefinedClusters[prev]->tStart;
					//st = 0;
					//st = 1;
				}
				else if (RefinedClusters[cur]->tStart > RefinedClusters[prev]->tEnd) {
					ts = RefinedClusters[prev]->tEnd;
					te = RefinedClusters[cur]->tStart;
					//st = 1;
					//st = 0;
				}
				else {
					c++;
					continue; // No need to refine the space!					
				}
				//cerr << "btwn  p: " << p << " h: " << h << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
				if (qe > qs and te > ts) {
					SpaceLength = max(qe - qs, te - ts); 
					//cerr << "SpaceLength: " << SpaceLength << "st: " << st << endl; 
					//used to be 100000; mapping contigs requires larger threshold;
					if (SpaceLength <= 1500000 and RefinedClusters[cur]->chromIndex == RefinedClusters[prev]->chromIndex) {
						// btwnClusters have GenomePos, st, matches, coarse
						// This function also set the "coarse" flag for RefinedClusters[cur]
						RefineBtwnSpace(RefinedClusters[cur], smallOpts, genome, read, strands, qe, qs, te, ts, st);
					}
				}
				c++;
			}
			//
			// Find matches at the right end;
			//
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
			//cerr << "right  p: " << p << " h: " << h << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
			if (qe > qs and te > ts) {
				SpaceLength = max(qe - qs, te - ts); 
				if (SpaceLength < 50000 and te+500 < genome.lengths[RefinedClusters[rh]->chromIndex]) { // used (1000, 6000)
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
					RefineBtwnSpace(RefinedClusters[rh], smallOpts, genome, read, strands, qe, qs, te, ts, st, lrts, lrlength);
				}
				// else {
				// 	RefineBtwnSpace(RefinedClusters[rh], smallOpts, genome, read, strands, qe, qe-1000, te, te-1000, st);
				// }				
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
			//cerr << "left  p: " << p << " h: " << h << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
			if (qe > qs and te > ts) {
				SpaceLength = max(qe - qs, te - ts);
				if (SpaceLength < 50000 and te+500 < genome.lengths[RefinedClusters[lh]->chromIndex]) { // used (1000, 6000)
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
					RefineBtwnSpace(RefinedClusters[lh], smallOpts, genome, read, strands, qe, qs, te, ts, st, lrts, lrlength);
				}	
				// else {
				// 	RefineBtwnSpace(RefinedClusters[lh], smallOpts, genome, read, strands, qe, qe-1000, te, te-1000, st);						
				// }			
			}
			// timing.Tick("refine_btwnclusters");
	// 	}
	// }
}
#endif