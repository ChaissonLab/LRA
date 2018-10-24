#ifndef MAP_READ_H_
#define MAP_READ_H_
#include "MMIndex.h"
#include "Genome.h"
#include "Read.h"
#include "Options.h"
#include "CompareLists.h"
#include "Sorting.h"
#include "TupleOps.h"
#include <iostream>
#include <algorithm>

#include <sstream>
#include "IndexedSeed.h"
#include "MergeSplit.h"

#include "Clustering.h"
#include "seqan/seeds.h"
#include "seqan/align.h"
#include "AffineOneGapAlign.h"

// Print results to stdout.
typedef seqan::String<seqan::Dna> TSequence;                 // sequence type
//typedef seqan::Infix<seqan::String<seqan::Dna> >::Type  Substring;
typedef seqan::Infix<char* >::Type  Substring;

typedef seqan::Align<Substring, seqan::ArrayGaps> TAlign;     // align type
typedef seqan::Row<TAlign>::Type TRow;                 // gapped sequence typ

typedef seqan::Align<TSequence, seqan::ArrayGaps> TSeqAlign;   
typedef seqan::Row<TSeqAlign>::Type TSeqRow;                 // gapped sequence typ


void SwapStrand(Read &read, Options &opts, GenomePairs &matches) {
	for (int m=0; m < matches.size(); m++) {
		matches[m].first.pos = read.length - (matches[m].first.pos + opts.globalK - 1);
	}
}



template <typename T>
void StoreSubset(vector<T> &a, vector<int> &idx) {
	int c, i;
	vector<T> t;
	for (c=0,i=0; i < idx.size(); i++, c++) {
		t.push_back(a[idx[i]]);
	}
	a=t;
}

int Matched(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te) {
	return min(qe-qs, te-ts);
}

void SetMatchAndGaps(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te, int &m, int &qg, int &tg) {
	m=Matched(qs, qe, ts, te);
	qg=qe-qs-m;
	tg=te-ts-m;
}


void RemoveOverlappingClusters(vector<Cluster> &clusters, Options &opts) {
	int a=0;
	int ovp=a;
	if (clusters.size() == 0) {
		return;
	}
	while (a < clusters.size() - 1) {
		ovp=a+1;
		float num=1.0;
		float denom=1.0;

		while ( ovp < clusters.size() and
						clusters[a].Overlaps(clusters[ovp], 0.8 ) ) {
			ovp++;
		}
		if (ovp - a > opts.maxCandidates) {
			for (int i=a+opts.maxCandidates; i < ovp; i++) {
				clusters[i].matches.clear();
			}
		}
		a=ovp;
	}
	int c=0;
	for (int i=0; i < clusters.size(); i++) {
		if (clusters[i].matches.size() > 0) {
			clusters[c] = clusters[i];
			c++;
		}
	}
	clusters.resize(c);
}

void SimpleMapQV(vector<Alignment*> &alignments) {
	int a=0;
	int ovp=a;
	if (alignments.size() == 0) {
		return;
	}
	if (alignments.size() == 1) {
		alignments[a]->mapqv=255;
		return;
	}
	while (a < alignments.size() - 1) {
		ovp=a+1;
		float num=1.0;
		float denom=1.0;

		while (ovp < alignments.size()  and
					 alignments[a]->Overlaps(*alignments[ovp], 0.9) ) {
			int nmDiff = alignments[a]->nm - alignments[ovp]->nm;
			if (nmDiff < 10) {
				denom+= pow(0.5,nmDiff);
			}
			ovp++;
		}
		if (ovp == a+1){
			alignments[a]->mapqv=255;
		}
		else {
			// comparing float ok because it is potentially not set.
			if (denom == 1.0) {
				alignments[a]->mapqv=255;
			}
			else {
				alignments[a]->mapqv=-10*log10(1-num/denom);
			}
		}
		a=ovp;
	}			
}

typedef seqan::Iterator<TSeqRow>::Type TRowIterator;
typedef char TChar;                             // character type

int AlignSubstrings(char *qSeq, GenomePos &qStart, GenomePos &qEnd, 
										char *tSeq, GenomePos &tStart, GenomePos &tEnd, 
										vector<int> &scoreMat, vector<Arrow> &pathMat, Alignment &aln) {
	
	int qLen = qEnd-qStart;
	int tLen = tEnd-tStart;
	int drift = abs(qLen - tLen);
	int k = max(7, drift+1);
	

	
	/*
	int score = KBandAlign(&qSeq[qStart], qEnd-qStart, &tSeq[tStart], tEnd-tStart, 
												 -5,3,2,2, k, // make these smart later.
												 scoreMat, pathMat, aln);*/

	string readSeq(&qSeq[qStart], qEnd-qStart);
	string chromSeq(&tSeq[tStart],tEnd-tStart);
	TSeqAlign align;
	/*
	cout << "Aligning " << endl;
	cout << readSeq << endl;
	cout << chromSeq << endl;
	*/
	int score = AffineOneGapAlign(readSeq, chromSeq, 4, -4, -3, 15, aln);
	/*
	seqan::resize(seqan::rows(align), 2);
	seqan::assignSource(seqan::row(align, 0), readSeq.c_str());
	seqan::assignSource(seqan::row(align, 1), chromSeq.c_str());
	int score = seqan::globalAlignment(align, seqan::Score<int, seqan::Simple>(4, -4, -13, -0), -k, k);

	cout << qStart << "\t" << tStart << "\t" << score << "\t" << qEnd-qStart << "\t" << tEnd - tStart << "\t" << k << endl;
	cout << readSeq << endl;
	cout << chromSeq << endl;
	
	cout << align << endl;

	TSeqRow & row1 = seqan::row(align, 0);
	TSeqRow & row2 = seqan::row(align, 1);
	int q=0, t=0;
	int rowLen = seqan::length(align);
	int i=0;


	TRowIterator qit = begin(row1),
		qitEnd = end(row1),
		tit = begin(row2);
	while ( qit != qitEnd ) {
		while(qit  < qitEnd and (seqan::isGap(qit) == true or seqan::isGap(tit) == true)) {
			if (seqan::isGap(qit)) { t+=1;}
			if (seqan::isGap(tit)) { q+=1;}
			qit++;
			tit++;
		}
		if (qit != qitEnd and seqan::isGap(qit) == false and seqan::isGap(tit) == false) {
			int bqs = q;
			int bts = t;
			int is=i;
			while(qit != qitEnd and seqan::isGap(qit) == false and seqan::isGap(tit) == false) {
					qit++;
					tit++;
					q++;
					t++;
			}
			aln.blocks.push_back(Block(bqs, bts, q-bqs));
		}
	}
	*/
	return score;
}

class SortClusterBySize {
public:
	bool operator()(Cluster &a, Cluster &b) {
		return a.matches.size() > b.matches.size();
	}
};

class SortAlignmentsByMatches {
public:
	bool operator()(const Alignment *a, const Alignment *b) const {
		return a->nm > b->nm;
	}
};

void RankClustersByScore(vector<Cluster> &clusters) {
	sort(clusters.begin(), clusters.end(), SortClusterBySize());
}

int SetStrand(Read &read, Genome &genome, Options &opts, GenomePairs &matches) { 
	int nSame=0;
	int nDifferent=0;
	for (int m=0; m< matches.size(); m++) {
		int chromIndex = genome.header.Find(matches[m].second.pos);
		char *chrom=genome.seqs[chromIndex];
		int chromPos = matches[m].second.pos - genome.header.pos[chromIndex];
		GenomeTuple readTup, genomeTup;
		StoreTuple(read.seq, matches[m].first.pos, opts.globalK, readTup);
		StoreTuple(chrom, chromPos, opts.globalK, genomeTup);
		if (readTup.t == genomeTup.t) {
			nSame++;
		}
		else {
			nDifferent++;
		}
	}
	if (nSame > nDifferent) {
		return 0;
	}
	else {
		return 1;
	}
}
template<typename T>
void SwapReadCoordinates(vector<T> &matches,
												 GenomePos readLength, GenomePos kmer){

	for (int i=0; i < matches.size(); i++) {
		matches[i].first.pos = readLength - (matches[i].first.pos+ kmer);
	}
}

void ReverseClusterStrand(Read &read, Genome &genome, Options &opts, 
											vector<Cluster> &clusters) {
	for (int c = 0; c < clusters.size(); c++) {
			SwapStrand(read, opts, clusters[c].matches);
			clusters[c].strand = 1;
	}
}


void SetClusterStrand(Read &read, Genome &genome, Options &opts, 
											vector<Cluster> &clusters) {
	for (int c = 0; c < clusters.size(); c++) {
		clusters[c].strand = SetStrand(read, genome, opts, clusters[c].matches);
		if (clusters[c].strand == 1) {
			SwapStrand(read, opts, clusters[c].matches);
		}
	}
}
template<typename T>
void UpdateBoundaries(T &matches, 
											GenomePos &qStart, GenomePos &qEnd, 
											GenomePos &tStart, GenomePos &tEnd) {
	for (int i =0; i< matches.size(); i++) {
		qStart=min(qStart, matches[i].first.pos);
		qEnd=max(qEnd, matches[i].first.pos);
		tStart=min(tStart, matches[i].second.pos);
		tEnd=max(tEnd, matches[i].second.pos);
	}
}
void RemoveEmptyClusters(vector<Cluster> &clusters, int minSize=1) {
	int cCur=0;
	for(int c=0; c<clusters.size(); c++) {
		if (clusters[c].tEnd== 0 or clusters[c].matches.size() < minSize ) {
			continue;
		}
		else {
					clusters[cCur] = clusters[c];
					cCur++;
		}
	}
	if (cCur < clusters.size() ) {
		clusters.resize(cCur);
	}
}

void MergeAdjacentClusters(ClusterOrder &order, Genome &genome, Options &opts) {
	int c=0;
	int cn=0;
	c=0;
	while(c< order.size()) {
		cn=c+1;
		int curEndChrom = genome.header.Find(order[c].tEnd);
		while (cn < order.size()) {
			int nextStartChrom = genome.header.Find(order[cn].tStart);
			int gap;
			gap = abs((int)((int)(order[cn].tStart - order[cn].qStart) - (int)(order[c].tEnd-order[c].qEnd)));

			if (nextStartChrom == curEndChrom and
					gap < opts.maxGap and
					order[c].strand == order[cn].strand ) {
				order[c].matches.insert(order[c].matches.end(), 
																order[cn].matches.begin(),
																order[cn].matches.end());
				order[c].qEnd = order[cn].qEnd;
				order[c].tEnd = order[cn].tEnd;
				order[cn].tEnd=0;
				cn++;
			}
			else {

				int cn2=cn;
				int MAX_AHEAD=10;
				while (cn2 < order.size() and 
							 cn2-cn < MAX_AHEAD and 
							 nextStartChrom == curEndChrom and 
							 order[c].strand == order[cn2].strand ) {
					gap = abs((int)((int)(order[cn2].tStart - order[cn2].qStart) - (int)(order[c].tEnd-order[c].qEnd)));
					nextStartChrom = genome.header.Find(order[cn2].tStart);
					if (gap < opts.maxGap) {
						break;
					}
					cn2++;
				}
				if (cn2 < order.size() and cn2 - cn < MAX_AHEAD and cn2 > cn and gap < opts.maxGap) {
					cn=cn2;
					order[c].matches.insert(order[c].matches.end(), 
																	order[cn].matches.begin(),
																	order[cn].matches.end());
					order[c].qEnd = order[cn].qEnd;
					order[c].tEnd = order[cn].tEnd;
					order[cn].tEnd=0;
					cn++;
				}
				else {
					break;
				}
			}
		}
		c=cn;
	}
}

void MergeOverlappingClusters(ClusterOrder &order) {
	int cCur = 0;
	while(cCur < order.size()){
		int cNext;
		
		cNext = cCur + 1;
		while ( cNext < order.size() and
						order[cNext].OverlapsPrevious(order[cCur])) {
			order[cCur].matches.insert(order[cCur].matches.end(),
															order[cNext].matches.begin(),
															order[cNext].matches.end());
			order[cCur].UpdateBoundaries(order[cNext]);
			//
			// Signal to remove cm;
			//
			order[cNext].start=0;
			order[cNext].end=0;
			cNext+=1;
		}
		cCur=cNext;
	}
	//
	// Remove merged clusters.
	//
	RemoveEmptyClusters(*order.clusters);
}

void RefineSubstrings(char *read,   GenomePos readSubStart, GenomePos readSubEnd, 
											char *genome, GenomePos genomeSubStart, GenomePos genomeSubEnd, 
											vector<int> &scoreMat, vector<Arrow> &pathMat, 
											Alignment &aln) {

	aln.blocks.clear();
	AlignSubstrings( read, readSubStart, readSubEnd, genome, genomeSubStart, genomeSubEnd, scoreMat, pathMat, aln);
	for (int b = 0; b < aln.blocks.size(); b++) {
		aln.blocks[b].qPos += readSubStart;
		aln.blocks[b].tPos += genomeSubStart;
	}
	
}


void MapRead(Read &read, 
						Genome &genome,
						vector<GenomeTuple> &genomemm,
						LocalIndex &glIndex,
						 Options &opts,
						 ostream *output,
						 pthread_mutex_t *semaphore=NULL) {
						 

	vector<GenomeTuple> readmm;
	vector<pair<GenomeTuple, GenomeTuple> > matches;
								
	if (opts.storeAll) {
		Options allOpts = opts;
		allOpts.globalW=1;
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length,
																				allOpts.globalK, allOpts.globalW, readmm);			
	}
	else {
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, opts.globalK, opts.globalW, readmm);
	}
	sort(readmm.begin(), readmm.end());
	CompareLists(readmm, genomemm, matches, opts);

	DiagonalSort<GenomeTuple>(matches);

	CleanOffDiagonal(matches, opts);

	vector<Cluster> clusters;

	StoreDiagonalClusters(matches, clusters, opts, false);

	AntiDiagonalSort<GenomeTuple>(matches, genome.GetSize());

	SwapStrand(read, opts, matches);
	
	vector<Cluster> revClusters;

	StoreDiagonalClusters(matches, revClusters, opts, false, 1);
	//
	// Add pointers to seq that make code more readable.
	//
	char *readRC;

	CreateRC(read.seq, read.length, readRC);
	char *strands[2] = { read.seq, readRC };

	clusters.insert(clusters.end(), revClusters.begin(), revClusters.end());

	ClusterOrder clusterOrder(&clusters);
	MergeAdjacentClusters(clusterOrder, genome, opts);
	int cl=0;
	int cn;
	for(cl=0, cn=0; cn < clusters.size(); cn++) {
		if (clusters[cn].tEnd > 0) {
			clusters[cl] = clusters[cn];
			cl++;
		}
	}
	clusters.resize(cl);
	
	ClusterOrder reducedOrder(&clusters);
	vector<Cluster> reducedClusters;
	for (int i = 0; i < reducedOrder.size(); i++) {
		reducedClusters.push_back(reducedOrder[i]);
	}
	clusters= reducedClusters;
	sort(clusters.begin(), clusters.end(), OrderClusterBySize());	
	/*	ofstream clust("clusters.tab");
	for (int c =0; c < clusters.size(); c++) {
		for (int m=0; m < clusters[c].matches.size(); m++) {
			clust << clusters[c].matches[m].second.pos << "\t" << clusters[c].matches[m].first.pos << "\t" << c << "\t" << clusters[c].strand << endl;
		}
	}
	clust.close();
	*/


  RemoveOverlappingClusters(clusters, opts);

	


	//
	// Build local index for refining alignments.
	//
	LocalIndex forwardIndex(glIndex);
	LocalIndex reverseIndex(glIndex);

	LocalIndex *localIndexes[2] = {&forwardIndex, &reverseIndex};
	forwardIndex.IndexSeq(read.seq, read.length);
	reverseIndex.IndexSeq(readRC, read.length);

	Options smallOpts = opts;
	smallOpts.globalK=glIndex.k;
	smallOpts.globalW=glIndex.w;
	smallOpts.globalMaxFreq=6;
	smallOpts.maxDiag=25;
	smallOpts.minDiagCluster=3;

	Options tinyOpts = opts;
	tinyOpts.globalMaxFreq=3;
	tinyOpts.maxDiag=5;
	tinyOpts.minDiagCluster=2;
	tinyOpts.globalK=smallOpts.globalK-3;


	//
	// Merge overlapping clusters
	//
	RemoveEmptyClusters(clusters, opts.minClusterSize);
	


	if (opts.mergeGapped) {
		ClusterOrder clusterOrder(&clusters);
		clusterOrder.Sort();
		MergeOverlappingClusters(clusterOrder);
	}
	
	vector<Cluster> refinedClusters(clusters.size());

	vector<Alignment*> alignments;

	for (int c = 0; c < clusters.size(); c++) {

		if (clusters[c].start == clusters[c].end) {
			continue;
		}			
		//
		// Get the boundaries of the cluster in both sequences.
		//
		CartesianTargetSort<GenomeTuple>(clusters[c].matches.begin(), 
																		 clusters[c].matches.end());
		int nMatch = clusters[c].matches.size();
		GenomePos tPos=clusters[c].matches[0].second.pos;
		int firstChromIndex   = genome.header.Find(tPos);
		int lastChromIndex;
		if (nMatch > 1 ) {
			tPos = clusters[c].matches[nMatch-1].second.pos;
			lastChromIndex = genome.header.Find(tPos);
		} else { lastChromIndex = firstChromIndex; }
		clusters[c].chromIndex = firstChromIndex;
		if (firstChromIndex != lastChromIndex ) {
			clusters[c].matches.clear();
			continue;
		}
					

		//
		// Make the anchors reference this chromosome for easier bookkeeping 
		//
		GenomePos chromOffset = genome.header.pos[firstChromIndex];
		for (int m=0; m < clusters[c].matches.size(); m++) {
			clusters[c].matches[m].second.pos-=chromOffset;
		}
		//
		// Get shorthand access to alignment boundaries.
		//

		GenomePos readClusterStart, readClusterEnd, genomeClusterStart, genomeClusterEnd;
		readClusterStart = clusters[c].matches[0].first.pos;
		genomeClusterStart = clusters[c].matches[0].second.pos + chromOffset;

		int cl = clusters[c].matches.size();
		readClusterEnd = clusters[c].matches[cl-1].first.pos + opts.globalK;
		genomeClusterEnd = clusters[c].matches[cl-1].second.pos + opts.globalK + chromOffset;

		int ls, le;

		GenomePos chromEndOffset   = genome.header.GetNextOffset(genomeClusterEnd);
		// Search region starts in window, or beginning of chromosome
		GenomePos wts, wte;
		if ( chromOffset + opts.window > genomeClusterStart ) {
			wts = chromOffset;
		}
		else {
			wts = genomeClusterStart - opts.window;
		}
				
		if (genomeClusterEnd + opts.window > chromEndOffset) {
			wte = chromEndOffset-1;
		}
		else {
			wte = genomeClusterEnd + opts.window;
		}
			
		ls = glIndex.LookupIndex(wts);
		le = glIndex.LookupIndex(wte);
				
			
		// 
		// Get quick access to the local index
		//
		LocalIndex *readIndex;
		readIndex = localIndexes[clusters[c].strand];

		int lmIndex=0;
		GenomePos alnStart=-1;
		GenomePos alnEnd=0;

		for (int lsi=ls; lsi <= le; lsi++) {
			//
			// Find the coordinates in the cluster that start in this local index.
			//
			GenomePos genomeLocalIndexStart = 
				glIndex.seqOffsets[lsi]  - 
				chromOffset;
			GenomePos genomeLocalIndexEnd   = 
				glIndex.seqOffsets[lsi+1] - 1 - chromOffset;

			int matchStart = 
				CartesianTargetLowerBound<GenomeTuple>(clusters[c].matches.begin(),
																							 clusters[c].matches.end(),
																							 genomeLocalIndexStart);

			int matchEnd   = 
				CartesianTargetUpperBound<GenomeTuple>(clusters[c].matches.begin(),
																							 clusters[c].matches.end(),
																							 genomeLocalIndexEnd);

			//
			// If there is no overlap with this cluster
			if (matchStart >= clusters[c].matches.size()) {
				continue;
			}
			GenomePos readStart = clusters[c].matches[matchStart].first.pos;
			if (lsi == ls) {
				if (readStart < opts.window) {
					readStart = 0;
				}
				else {
					readStart -= opts.window;
				}
			}
			GenomePos readEnd;
			if (matchEnd > matchStart) {
				readEnd = clusters[c].matches[matchEnd-1].first.pos;
			}
			else {
				readEnd = clusters[c].matches[matchStart].first.pos + opts.globalK;
			}
			//
			// Expand boundaries of read to match.
			if (lsi == le) {
				if (readEnd + opts.window > read.length) {
					readEnd = read.length; 
				}
				else { 
					readEnd += opts.window;	
				}
			}			
				
			//
			// Find the boundaries where in the query the matches should be added.
			//

			int queryIndexStart = readIndex->LookupIndex(readStart);
			int queryIndexEnd   = 
				readIndex->LookupIndex(min(readEnd, 
																	 (GenomePos) read.length-1));
			assert(queryIndexEnd < readIndex->seqOffsets.size()+1);

			for (int qi = queryIndexStart; qi <= queryIndexEnd; ++qi){ 
				LocalPairs smallMatches;

				GenomePos qStartBoundary = readIndex->tupleBoundaries[qi];
				GenomePos qEndBoundary   = readIndex->tupleBoundaries[qi+1];
				GenomePos readSegmentStart= readIndex->seqOffsets[qi];
				GenomePos readSegmentEnd  = readIndex->seqOffsets[qi+1];

				CompareLists<LocalTuple>(readIndex->minimizers.begin() + 
																 qStartBoundary,

																 readIndex->minimizers.begin() +
																 qEndBoundary,

																 glIndex.minimizers.begin()+
																 glIndex.tupleBoundaries[lsi], 

																 glIndex.minimizers.begin()+
																 glIndex.tupleBoundaries[lsi+1], 

																 smallMatches, smallOpts);

				lmIndex+=smallMatches.size();


				//
				// Do local processing of matches to ensure the region that is searched returns reasonable anchors.
				//

				DiagonalSort<LocalTuple>(smallMatches);
				CleanOffDiagonal(smallMatches, smallOpts);					

				AppendValues<LocalPairs>(refinedClusters[c].matches, 
																 smallMatches.begin(), smallMatches.end(), 
																 readSegmentStart, genomeLocalIndexStart );
			}
				
		}

		refinedClusters[c].SetClusterBoundariesFromMatches();
		refinedClusters[c].strand = clusters[c].strand;
		refinedClusters[c].chromIndex = clusters[c].chromIndex;
		refinedClusters[c].coarse = c;
	}

	RemoveEmptyClusters(refinedClusters);

	for (int r =0; r < refinedClusters.size(); r++) {

		if (refinedClusters[r].matches.size() < opts.minRefinedClusterSize) {
			continue;
		}

		//
		// Clean local matches to reduce chaining burden.
		//
		DiagonalSort<GenomeTuple>(refinedClusters[r].matches);
		CleanOffDiagonal(refinedClusters[r].matches, smallOpts);

		if (refinedClusters[r].matches.size() == 0) {
			continue;
		}
		int k=glIndex.k;
		if (refinedClusters[r].matches.size() > read.length or 
				opts.refineLevel & REF_LOC == 0) {
			refinedClusters[r].matches= clusters[refinedClusters[r].coarse].matches;
			k=opts.globalK;
			continue;
		}

		//
		// At this point in the code, refinedClusters[r].matches is a list
		// of k-mer matches between the read and the genome. 
		
		//
		// This is where the code should go for merging colinear matches

		// Build SeedSet from refined (small) matches.
		seqan::SeedSet<IndSeed, seqan::Unordered> seedSet;


		// Jingwen: Instead of copying directly from refinedClusters[r].matches into the seed set, you can use your code to:
		//  1. Merge adjacent anchors (using "merge" from MergeSplit.h)
		//  2. Split overlapping anchors.
		// 
		// The container of matches is 
		//	vector<Cluster> refinedClusters(clusters.size());
		// Each of these has a vector 'matches', refinedClusters[r].matches
		// 'matches' is defined as 	GenomePairs matches;
		// GenomePairs is a list of GenomePair:
		// typedef vector<GenomePair > GenomePairs;
		// a GenomePair is a pair of indexes into a genome or a read (both called genomes)
		// it is defined as:
		// typedef pair<GenomeTuple, GenomeTuple> GenomePair;
		// You access the first using 'first', as in:
		// GenomePair p;
		// p.first.pos = 1000;
		// p.second.pos=2000;
		//   Each GenomeTuple contains a 'pos' , or the position of a k-mer in a genome (or read), and a tuple, which is just a k-mer. 
		//  When there is a GenomePair, the two tuples are the same!
		//

		// The k-mer that is used in the refined matches is store in 
		// smallOpts.globalK

		if (opts.mergeClusters) {

			// add merge split code here
			// The result should be to store the merged clusters in seedSet
			CartesianSort<GenomeTuple>(refinedClusters[r].matches);

			vector<vector<GenomePair*>> v;
			int h = 500; 
			merge (refinedClusters[r], v, h, k);

	//		vector<vector<GenomePair*>> v_prime;
	//		int g = 10;
	//		RemoveSomeSeed (v, v_prime, g, k);

			vector<vector<IndSeed>> splitHSeed;
			SplitH (v, splitHSeed, k);

			vector<vector<IndSeed>> splitVSeed;
			SplitV (splitHSeed, splitVSeed);

			AddInset (splitVSeed, seedSet);
		}
		else {
			for (int m=0; m< refinedClusters[r].matches.size(); m++) {
				seqan::addSeed(seedSet, 
											 seqan::Seed<IndexedSeed>(refinedClusters[r].matches[m].second.pos, 
																								refinedClusters[r].matches[m].first.pos, 
																								k),
											 seqan::Single());
			}
		}


		// Perform sparse chaining, uses time O(n log n).
		//
		// Merge anchors that are not overlapping

		//
		// The input is a set of seqan::Seed<seqan::Simple> > seeds
		//

		seqan::String<seqan::Seed<seqan::Simple> > chain;
		if (seqan::length(seedSet) > 0) {
			seqan::chainSeedsGlobally(chain, seedSet, seqan::SparseChaining());
		}
		
		//Debug code
		/*
		// Print results to stdout.
		const string filename4("/home/cmb-16/mjc/jingwenr/lra/lra_test/sparse_result.txt");  // file to store original_seed
		 // save the original seeds
    	FILE *fi = fopen(filename4.c_str(), "w");
   		SaveSparse (chain, fi);
   		fclose(fi);
	    */


		vector<GenomePair> tupChain;
		int qPrev=0, tPrev=0;
		int csg = seqan::length(chain);
		if (opts.mergeClusters) {
			//
			// Add small anchors to tupChain.
			//
			/*
				tupChain.push_back(MatchingPos(GenomeTuple(0, beginPositionV(chain[ch])),
																			 GenomeTuple(0, beginPositionH(chain[ch])), length_of_split_kmer));
			*/
		}
		else {
			for (int ch=0; ch < seqan::length(chain); ch++) {
				tupChain.push_back(GenomePair(GenomeTuple(0, beginPositionV(chain[ch])),
																			GenomeTuple(0, beginPositionH(chain[ch]))));
			}
		}
		if (tupChain.size() == 0) {
			refinedClusters[r].matches.clear();
			continue;
		}
		vector<Cluster> chainClust;
		Options diagOpts;
		diagOpts = smallOpts;
		diagOpts.maxDiag=15;
		diagOpts.minClusterSize=1;
		StoreDiagonalClusters(tupChain, chainClust, diagOpts, false, refinedClusters[r].strand);

		RemovePairedIndels(tupChain, chainClust, smallOpts);
			
		seqan::clear(chain);
		for (int m=0; m< tupChain.size(); m++) {
			seqan::append(chain, 
										seqan::Seed<seqan::Simple>(tupChain[m].second.pos,
																							 tupChain[m].first.pos, glIndex.k));

		}
		int prevq=0;
		int prevt=0;
		
		//		cout << "Chain is on " << chainClust.size() << " diagonals " << endl;
		/*
		for (int cc = 0; cc < chainClust.size(); cc++ ) {
			if (cc > 0) {
								cout << (int) (chainClust[cc].qStart - prevq - (chainClust[cc].tStart - prevt)) << "\t" 
						 << chainClust[cc].qStart - prevq << "\t" << chainClust[cc].qStart << "\t" << chainClust[cc].qEnd << "\tt " 
						 << chainClust[cc].tStart - prevt << "\t" << chainClust[cc].tStart << "\t" << chainClust[cc].tEnd  << endl;
			}

			prevq = chainClust[cc].qEnd;
			prevt = chainClust[cc].tEnd;
		}
		*/
		GenomePos chainGenomeStart = seqan::beginPositionH(chain[0]);
		GenomePos chainGenomeEnd   = seqan::endPositionH(chain[seqan::length(chain)-1]);

		GenomePos chainReadStart = seqan::beginPositionV(chain[0]);
		GenomePos chainReadEnd   = seqan::endPositionV(chain[seqan::length(chain)-1]);


		GenomePos globalStart = refinedClusters[r].tStart;
		int chromIndex  = refinedClusters[r].chromIndex;
			
		//
		// Create subsequences that will be used to generate the alignment.  Gaps should be inserted 
		// with respect to an offset from chainGenomeStart and chainReadStart
		//
		GenomePos genomeAlnOffset = chainGenomeStart;
		GenomePos readAlnOffset = chainReadStart;
		vector<GenomeTuple> gapReadTup, gapGenomeTup;
		GenomePairs gapPairs;
		Options gapOpts=opts;
		gapOpts.globalMaxFreq=5;
		gapOpts.globalK=7;
		vector< seqan::String<seqan::Seed<seqan::Simple> > > refinedChains(seqan::length(chain)-1);
		// Build SeedSet.
		seqan::SeedSet<seqan::Seed<seqan::Simple>, seqan::Unordered> gapSeedSet;


		vector<int> scoreMat;
		vector<Arrow> pathMat;

		int chainLength = seqan::length(chain);
		for (int c = 0; chainLength > 0 and c < seqan::length(chain)-1; c++) {
			GenomePos curGenomeEnd = seqan::endPositionH(chain[c]);
			GenomePos nextGenomeStart = seqan::beginPositionH(chain[c+1]);

			GenomePos curReadEnd = seqan::endPositionV(chain[c]);
			GenomePos nextReadStart = seqan::beginPositionV(chain[c+1]);
				
			int rg=nextReadStart-curReadEnd, gg=nextGenomeStart-curGenomeEnd;
			int mg=min(rg, gg);
			int rm=rg-mg;
			int gm=gg-mg;
			rg-=mg;
			gg-=mg;

			GenomePos subreadLength = nextReadStart-curReadEnd;
			GenomePos subgenomeLength = nextGenomeStart-curGenomeEnd;

			if (nextReadStart > curReadEnd and nextGenomeStart > curGenomeEnd) {

				if (subreadLength > 50 and 
						subgenomeLength > 50 and opts.refineLevel & REF_DYN ) {

					GenomePos maxLen = max(subreadLength, subgenomeLength);
					if (maxLen < 500) {
						tinyOpts.globalK=5;
					}
					else if (maxLen < 2000) {
						tinyOpts.globalK=5;
					}
					else {
						tinyOpts.globalK=9;
					}
					gapGenomeTup.clear();
					gapReadTup.clear();
					gapPairs.clear();

					//
					// Find matches between read and reference in the coordinate space of read and chromosome
					//
					assert(curGenomeEnd < genome.lengths[chromIndex]);
					assert(curGenomeEnd + subgenomeLength < genome.lengths[chromIndex]);
					StoreMinimizers<GenomeTuple, Tuple>( genome.seqs[chromIndex] + curGenomeEnd,
																							 subgenomeLength, tinyOpts.globalK, 1, gapGenomeTup, false);

					sort(gapGenomeTup.begin(), gapGenomeTup.end());
					StoreMinimizers<GenomeTuple, Tuple>( strands[refinedClusters[r].strand] + curReadEnd,
																							 subreadLength, tinyOpts.globalK, 1, gapReadTup, false);
					sort(gapReadTup.begin(), gapReadTup.end());
					CompareLists(gapReadTup.begin(),
											 gapReadTup.end(),
											 gapGenomeTup.begin(),
											 gapGenomeTup.end(),
											 gapPairs, gapOpts);
					//
					// Remove egregious off diagonal seeds
					//
					DiagonalSort<GenomeTuple>(gapPairs);
					
					int tinyDiagStart = curGenomeEnd - curReadEnd;
					int tinyDiagEnd   = nextGenomeStart - nextReadStart;
					int diagDiff = abs((int) tinyDiagStart - (int) tinyDiagEnd);

					CleanOffDiagonal(gapPairs, tinyOpts, 0, diagDiff);
					//					CartesianTargetSort<GenomeTuple>(gapPairs);

					for(int rm=0; rm < gapPairs.size(); rm++) {
						gapPairs[rm].first.pos  += curReadEnd;
						gapPairs[rm].second.pos += curGenomeEnd;
					}
					int gp=0;
					int nSaved =0;

					seqan::clear(gapSeedSet);
					int gpStart;
					int clusterIndex=0;

					seqan::SeedSet<seqan::Seed<seqan::Simple>, seqan::Unordered> gapSeedSet;
					for (int m=0; m< gapPairs.size(); m++) {
						seqan::addSeed(gapSeedSet, 
													 seqan::Seed<seqan::Simple>(gapPairs[m].second.pos, 
																											gapPairs[m].first.pos, tinyOpts.globalK),
													 seqan::Single());
					}
					seqan::String<seqan::Seed<seqan::Simple> > gapChain;
					if (seqan::length(gapSeedSet) > 0) {
						seqan::chainSeedsGlobally(gapChain, gapSeedSet, seqan::SparseChaining());
					}
					int pre=seqan::length(gapChain);
					RemovePairedIndels(curReadEnd, curGenomeEnd, nextReadStart, nextGenomeStart, gapChain, tinyOpts);

					refinedChains[c]=gapChain;						
				}
			}
		}
			
		//
		// Refine and store the alignment
		//
		//
		// The alignment is on a substring that starts at the beginning of the first chain.
		//
		GenomePos alnReadPos = seqan::beginPositionV(chain[0]);
		GenomePos alnRefPos  = seqan::beginPositionH(chain[0]);
		Alignment *alignment = new Alignment(strands[refinedClusters[r].strand], 
																				 read.seq,
																				 read.length, read.name,
																				 refinedClusters[r].strand,
																				 genome.seqs[chromIndex],  
																				 genome.lengths[chromIndex],
																				 genome.header.names[chromIndex],
																				 chromIndex) ;

		alignments.push_back(alignment);
		/*
		ofstream refc("refc.tab");
		ofstream reft("reft.tab");
		for (int c=0; c < refinedChains.size(); c++) {
			for (int m=0; m < seqan::length(refinedChains[c]); m++) {
				reft <<  seqan::beginPositionV(refinedChains[c][m]) << "\t" << 
					seqan::beginPositionH(refinedChains[c][m]) << endl;
			}
		}
		*/	
		for (int c = 0; chainLength> 0 and  c < chainLength-1; c++) {
				
			//
			// Chain is with respect to full sequence
			//
			GenomePos curGenomeEnd     = seqan::endPositionH(chain[c]);
			GenomePos nextGenomeStart  = seqan::beginPositionH(chain[c+1]);

			GenomePos curReadEnd       = seqan::endPositionV(chain[c]);
			GenomePos nextReadStart    = seqan::beginPositionV(chain[c+1]);
			int curRefinedReadEnd      = curReadEnd;
			int curRefinedGenomeEnd    = curGenomeEnd;
			int nextRefinedReadStart   = nextReadStart;
			int nextRefinedGenomeStart = nextGenomeStart;

			alignment->blocks.push_back(Block(seqan::beginPositionV(chain[c]),
																				seqan::beginPositionH(chain[c]), glIndex.k));
			//			string curAnchor = string(genome.seqs[chromIndex], seqan::beginPositionV(chain[c]), glIndex.k );

			for (int cs = 0; cs < seqan::length(refinedChains[c]); cs++) {
				//
				// Refined anchors are with respect to the chained sequence
				nextRefinedReadStart   = seqan::beginPositionV(refinedChains[c][cs]);
				nextRefinedGenomeStart = seqan::beginPositionH(refinedChains[c][cs]);
				

				int m, rg, gg;
				SetMatchAndGaps(curRefinedReadEnd, nextRefinedReadStart,
												curRefinedGenomeEnd, nextRefinedGenomeStart, m, rg, gg);

				if (m > 0) {
					Alignment betweenAnchorAlignment;
					if (opts.refineLevel & REF_DP) {						
						RefineSubstrings(strands[refinedClusters[r].strand], curRefinedReadEnd, nextRefinedReadStart,
														 genome.seqs[chromIndex], curRefinedGenomeEnd, nextRefinedGenomeStart,
														 scoreMat, pathMat, betweenAnchorAlignment);
						alignment->blocks.insert(alignment->blocks.end(), 
																		 betweenAnchorAlignment.blocks.begin(), 
																		 betweenAnchorAlignment.blocks.end());
						betweenAnchorAlignment.blocks.clear();
					}
				}

				curRefinedReadEnd   = seqan::endPositionV(refinedChains[c][cs]);
				curRefinedGenomeEnd = seqan::endPositionH(refinedChains[c][cs]);
				alignment->blocks.push_back(Block(nextRefinedReadStart, nextRefinedGenomeStart, 
																					curRefinedReadEnd - nextRefinedReadStart));
			}

			// Add the last gap, or the only one if no refinements happened here.
							 
			int match, readGap, genomeGap;
			SetMatchAndGaps(curRefinedReadEnd, nextReadStart,
											curRefinedGenomeEnd, nextGenomeStart, match, readGap, genomeGap);
				
			if (match > 0) {
				if (opts.refineLevel & REF_DP) {
					Alignment aln;
					assert(curRefinedReadEnd < read.length);
					assert(nextReadStart <read.length);
					assert(curRefinedGenomeEnd < genome.lengths[chromIndex]);
					assert(nextGenomeStart < genome.lengths[chromIndex]);
					RefineSubstrings(strands[refinedClusters[r].strand], curRefinedReadEnd, nextReadStart,
													 genome.seqs[chromIndex], curRefinedGenomeEnd, nextGenomeStart,
													 scoreMat, pathMat, aln);
					alignment->blocks.insert(alignment->blocks.end(), aln.blocks.begin(), aln.blocks.end());
					aln.blocks.clear();			
				}		
			}
		}
		//		refc.close();
		//			rdc.close();
			
		alignment->blocks.push_back(Block(seqan::beginPositionV(chain[chainLength-1]),
																			seqan::beginPositionH(chain[chainLength-1]), 
																			glIndex.k));

		int nm=0;
		for(int b=0; b< alignment->blocks.size(); b++) {
			nm+= alignment->blocks[b].length;
		}
		alignment->nblocks=seqan::length(chain);
	}
	for (int a=0; a < alignments.size(); a++) {
		alignments[a]->CalculateStatistics();
	}
	sort(alignments.begin(), alignments.end(), SortAlignmentsByMatches());

	SimpleMapQV(alignments);

	if (semaphore != NULL) {
		pthread_mutex_lock(semaphore);
	}
	for (int a=0; a < min(opts.bestn, (int) alignments.size()); a++ ){
		if (opts.printFormat == 'b') {
			alignments[a]->PrintBed(*output);
		}
		else if (opts.printFormat == 's') {
			alignments[a]->PrintSAM(*output, opts);
		}
		else if (opts.printFormat == 'p') {
			alignments[a]->PrintPairwise(*output);
		}
	}
	if (semaphore != NULL ) {
		pthread_mutex_init(semaphore, NULL);
	}

	//
	// Done with one read. Clean memory.
	//
	delete[] readRC;
	for (int a=0; a < alignments.size(); a++) {
		delete alignments[a];
	}
	read.Clear();
}


#endif
