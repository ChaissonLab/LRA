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
#include <thread>
#include "IndexedSeed.h"

#include "Clustering.h"
#include <thread>

#include "seqan/seeds.h"
#include "seqan/align.h"
using namespace std;
#include "AffineOneGapAlign.h"
#include "NaiveDP.h"
#include <iterator>

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
	int score = AffineOneGapAlign(readSeq, chromSeq, 4, -4, -3, 15, aln);
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
					order[c].strand == order[cn].strand and 
					order[c].Encompasses(order[cn],0.5) == false) {
				order[c].matches.insert(order[c].matches.end(), 
																order[cn].matches.begin(),
																order[cn].matches.end());
				order[c].qEnd = order[cn].qEnd;
				order[c].tEnd = order[cn].tEnd;
				order[cn].tEnd=0;
				order[cn].tStart=0;
				order[cn].matches.clear();
				cn++;
			}
			else {

				int cn2=cn;
				int MAX_AHEAD=10;
				while (cn2 < order.size() and 
							 cn2-cn < MAX_AHEAD and 
							 nextStartChrom == curEndChrom and 
							 order[c].strand == order[cn2].strand and 
							 order[c].Encompasses(order[cn2],0.5) == false) {
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
					order[cn].matches.clear();
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

void SeparateMatchesByStrand(Read &read, Genome &genome, int k,
														 vector<pair<GenomeTuple, GenomeTuple> > &allMatches, 
														 vector<pair<GenomeTuple, GenomeTuple> > &forMatches,
														 vector<pair<GenomeTuple, GenomeTuple> > &revMatches) {
	vector<bool> strand(allMatches.size());
	int nForward=0;
	for (int i=0; i < allMatches.size(); i++) {
		int readPos     = allMatches[i].first.pos;
		uint64_t refPos = allMatches[i].second.pos;
		char *genomePtr=genome.GlobalIndexToSeq(refPos);
		if (strncmp(&read.seq[readPos], genomePtr, k) == 0) {
			nForward++;
		}
		else {
			strand[i] = true;
		}
	}
	forMatches.resize(nForward);
	revMatches.resize(allMatches.size()-nForward);
	int i=0,r=0,f=0;
	for (i=0,r=0,f=0; i < allMatches.size(); i++) {
		if (strand[i] == 0) {
			forMatches[f] = allMatches[i];
			f++;
		}
		else {
			revMatches[r] = allMatches[i];
			r++;
		}
	}
}


//debug code Dont forget delete int &i
void MapRead(Read &read, 
						 Genome &genome,
						 vector<GenomeTuple> &genomemm,
						 LocalIndex &glIndex,
						 Options &opts,
						 ostream *output,
						 pthread_mutex_t *semaphore=NULL) {

	string baseName = read.name;


	for (int i=0; i < baseName.size(); i++) {	if (baseName[i] == '/') baseName[i] = '_';	}


	vector<GenomeTuple> readmm;
	vector<pair<GenomeTuple, GenomeTuple> > allMatches, forMatches, revMatches, matches;

	if (opts.storeAll) {
		Options allOpts = opts;
		allOpts.globalW=1;
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, allOpts.globalK, allOpts.globalW, readmm);			
	}
	else {
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, opts.globalK, opts.globalW, readmm);
	}
	sort(readmm.begin(), readmm.end());
	//
	// Add matches between the read and the genome.
	//
	CompareLists(readmm, genomemm, allMatches, opts);

	DiagonalSort<GenomeTuple>(allMatches);
	SeparateMatchesByStrand(read, genome, opts.globalK, allMatches, forMatches, revMatches);


	CleanOffDiagonal(forMatches, opts);

	if (opts.dotPlot) {
		ofstream clust("for-matches.dots");
		for (int m=0; m < forMatches.size(); m++) {
			clust << forMatches[m].first.pos << "\t" << forMatches[m].second.pos << "\t" << opts.globalK << "\t0\t0"<<endl;
		}
		clust.close();
		ofstream rclust("rev-matches.dots");
		for (int m=0; m < revMatches.size(); m++) {
			rclust << revMatches[m].first.pos << "\t" << revMatches[m].second.pos << "\t" << opts.globalK << "\t0\t0"<<endl;
		}
		rclust.close();

	}

	vector<Cluster> clusters;
	int forwardStrand=0;
	StoreDiagonalClusters(forMatches, clusters, opts, false, forwardStrand);


	AntiDiagonalSort<GenomeTuple>(revMatches, genome.GetSize());
	SwapStrand(read, opts, revMatches);
	CleanOffDiagonal(revMatches, opts);
	vector<Cluster> revClusters;
	
	int reverseStrand=1;
	StoreDiagonalClusters(revMatches, revClusters, opts, false, reverseStrand);
	//
	// Add pointers to seq that make code more readable.
	//
	char *readRC;

	CreateRC(read.seq, read.length, readRC);
	char *strands[2] = { read.seq, readRC };

	clusters.insert(clusters.end(), revClusters.begin(), revClusters.end());
	ClusterOrder clusterOrder(&clusters);
	if (opts.dotPlot) {
		ofstream clust("clusters-pre-merge.tab");
		for (int c =0; c < clusters.size(); c++) {
			for (int m=0; m < clusters[c].matches.size(); m++) {
				clust << clusters[c].matches[m].second.pos << "\t" << clusters[c].matches[m].first.pos << "\t" << opts.globalK << "\t" << c << "\t" << clusters[c].strand << endl;
			}
		}
		clust.close();
	
	}

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
	if (opts.dotPlot) {
		ofstream matchfile("long_matches.tab");
		for (int m =0; m < matches.size(); m++) {
			matchfile << matches[m].first.pos << "\t" << matches[m].second.pos << "\t" << opts.globalK << "\t0\t0" << endl;
		}
		matchfile.close();
		ofstream clust("clusters.tab");
		for (int c =0; c < clusters.size(); c++) {
			for (int m=0; m < clusters[c].matches.size(); m++) {
				clust << clusters[c].matches[m].second.pos << "\t" << clusters[c].matches[m].first.pos << "\t" << opts.globalK << "\t" << c << "\t" << clusters[c].strand << endl;
			}
		}
		clust.close();
	
	}

  RemoveOverlappingClusters(clusters, opts);

	
	//
	// Merge overlapping clusters
	//
	//RemoveEmptyClusters(clusters, opts.minClusterSize);
	
	if (opts.mergeGapped) {
		ClusterOrder clusterOrder(&clusters);
		clusterOrder.Sort();
		MergeOverlappingClusters(clusterOrder);
	}



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


	vector<Cluster> refinedClusters(clusters.size());

	vector<Alignment*> alignments;

	for (int c = 0; c < clusters.size(); c++) {

		if (clusters[c].start == clusters[c].end) {
			continue;
		}			
		if (opts.dotPlot) {
			stringstream outNameStrm;
			outNameStrm << baseName + "." << c << ".clust.dots";
			ofstream baseDots(outNameStrm.str().c_str());
			for (int m=0; m < clusters[c].matches.size(); m++) {
				baseDots << clusters[c].matches[m].second.pos << "\t" << clusters[c].matches[m].first.pos << "\t" << smallOpts.globalK << "\t" << c << "\t0" << endl;
			}
			baseDots.close();
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
		} else { 
			lastChromIndex = firstChromIndex; 
		}
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

		refinedClusters[c].SetClusterBoundariesFromMatches(opts);
		refinedClusters[c].strand = clusters[c].strand;
		refinedClusters[c].chromIndex = clusters[c].chromIndex;
		refinedClusters[c].coarse = c;
		if (opts.dotPlot) {
			stringstream outNameStrm;
			outNameStrm << baseName + "." << c << ".orig.dots";
			ofstream baseDots(outNameStrm.str().c_str());
			for (int m=0; m < refinedClusters[c].matches.size(); m++) {
				baseDots << refinedClusters[c].matches[m].second.pos << "\t" << refinedClusters[c].matches[m].first.pos << "\t" << smallOpts.globalK << "\t" << c << "\t0" << endl;
			}
			baseDots.close();
		}
			

	}
	//
	// Remove clusters under some minimal number of anchors. By default this is one. 
	//
	//RemoveEmptyClusters(refinedClusters);
	

	//cout << "refinedClusters.size(): " << refinedClusters.size() << endl;
	/*
	for (int r = 0; r < refinedClusters.size(); ++r) {
		cerr << "refinedClusters[" << r << "]: " << endl;
		for (int s = 0; s < refinedClusters[r].matches.size(); ++s) {
			cerr << refinedClusters[r].matches[s].first.pos << "\t" 
				 << refinedClusters[r].matches[s].second.pos << "\t"
				 << smallOpts.globalK << endl;
		}
	}
	*/


	for (int r = 0; r < refinedClusters.size(); ++r) {   
		//debug
		//cout <<"   refinedClusters[r]: " << r << "\t" << refinedClusters[r].matches.size() << endl;


		if (refinedClusters[r].matches.size() < opts.minRefinedClusterSize) {
			continue;
		}
		ofstream dotFile;
		if (opts.dotPlot) {
			stringstream outName;
			outName << baseName << "." << r << ".dots";
			dotFile.open(outName.str().c_str());
		}

		//
		// Clean local matches to reduce chaining burden.
		//
		DiagonalSort<GenomeTuple>(refinedClusters[r].matches);
		CleanOffDiagonal(refinedClusters[r].matches, smallOpts);

		if (opts.dotPlot) {
			stringstream outNameStrm;
			outNameStrm << baseName + "." << r << ".clean.dots";
			ofstream baseDots(outNameStrm.str().c_str());
			for (int m=0; m < refinedClusters[r].matches.size(); m++) {
				baseDots << refinedClusters[r].matches[m].second.pos << "\t" << refinedClusters[r].matches[m].first.pos << "\t" << smallOpts.globalK << "\t" << r << "\t0" << endl;
			}
			baseDots.close();
		}


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
		vector<Cluster> vt; // vt stores the split result of merged anchors


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

		/*
		// Debug code ----------- print out "seedSet"

		if (refinedClusters.size() == 9) {
		const string filename0("/home/cmb-16/mjc/jingwenr/lra/lra_test/test_lra1/Orignalseeds.txt");  
		FILE *frr = fopen(filename0.c_str(), "w");
		SaveOriginalSeed (refinedClusters[r], frr, opts.globalK); 
		fclose(frr);	
		}
		*/	

		/*
		// Debug code ----------- print out "seedSet"

			cerr << "1" << endl;
			const string filename("/home/cmb-16/mjc/jingwenr/lra/lra_test/TEST/Orignalseeds." + std::to_string(r) + ".txt");
			FILE *frr = fopen(filename.c_str(), "w");
			SaveOriginalSeed (refinedClusters[r], frr, opts.globalK); 
			fclose(frr);
			*/
		


		// debug code
		if (opts.mergeClusters) {
			vt.clear();
			// add merge split code here
			//------------------------------------------------------------------
			// Step 1: merge matches with diagnol difference smaller than maxDiag
			//------------------------------------------------------------------
			vector<Cluster> v;
			vector<Cluster> vq;
			Options mergeOpts=smallOpts;
			mergeOpts.maxGap = -1;
			mergeOpts.maxDiag = 20;

			StoreDiagonalClusters(refinedClusters[r].matches, v, mergeOpts, true); 
			if (v.size() != 0) {
				std::set<unsigned int> s1;     // s1 stores x boundary; s2 stores y boundary      
				std::set<unsigned int> s2;    //elements in s are arranged in strictly increasing order. and no duplicates
            
				/*
				// debug
				if (r == 2) {
				cout << "refinedClusters[r].matches.size(): " << refinedClusters[r].matches.size() << endl;
				}
				*/

				/*
				// Debug code ----------- print out "orignal"
			
				seqan::SeedSet<IndSeed, seqan::Unordered> seedSet1;
				for (unsigned int i = 0; i != v.size(); ++i) {
					seqan::addSeed(seedSet1, IndSeed(v[i].qStart, v[i].tStart, v[i].qEnd, v[i].tEnd, i), seqan::Single());
				} 
				for (TIterator tt = begin(seedSet1, seqan::Standard()); tt != end(seedSet1, seqan::Standard()); ++tt) {
					cerr << *tt << endl;
				}
			
				
				cout << "v.szie(): " << v.size() << endl;
				const string filename7("/home/cmb-16/mjc/jingwenr/lra/lra_test/TEST/Merge.txt");  
				FILE *fr = fopen(filename7.c_str(), "w");
				SaveseedSet (seedSet1, fr); 
				fclose(fr);
				*/	
			

				//---------------------------------------------------------------------------------
				// Step 2: sort matches by x first and then by y coordinates in each Cluster v[i](merged matches)
				//---------------------------------------------------------------------------------
				for (unsigned int i = 0; i < v.size(); ++i) {

					//insert x boundaries into s1 
					s1.insert(v[i].qStart);
					s1.insert(v[i].qEnd);

					// sort matches by x and y inside each merged matches cluster
					vector<pair<GenomeTuple, GenomeTuple>>::iterator begin = refinedClusters[r].matches.begin();
					if (v[i].start != 0) {
						std::advance(begin , v[i].start);
					}	 

					vector<pair<GenomeTuple, GenomeTuple>>::iterator end = refinedClusters[r].matches.begin();
					if (v[i].end != 0) {
						std::advance(end, v[i].end);
					}
					CartesianSort<GenomeTuple>(begin, end);            
				}

				//----------------------------------------------------------------------------
				// Step 3:  split based on x boundaries
				//----------------------------------------------------------------------------
				for (unsigned int i = 0; i < v.size(); ++i) {

					std::set<unsigned int>::iterator qs, qe; 
					qs = s1.lower_bound(v[i].qStart);
					qe = s1.lower_bound(v[i].qEnd);
					std::set<unsigned int>::iterator l = qs;
					std::advance(l, 1);
					unsigned int ii = v[i].start;
					unsigned int jj = ii;

					if (*l == *qe) {
						vq.push_back(Cluster(v[i].start, v[i].end));
					}
					else{
						while (ii != v[i].end && l != qe) { // end is end + 1!!!!!!
							while (refinedClusters[r].matches[ii].first.pos <= *l && ii != v[i].end) {
								if (refinedClusters[r].matches[ii].first.pos + smallOpts.globalK - 1 >= *l) {
									// split 
									if (jj + 1 <= ii) {
										assert(jj < refinedClusters[r].matches.size());
										vq.push_back(Cluster(jj, ii));
									}
									vq.push_back(Cluster(ii, ii + 1));
									jj = ii + 1;
								}
								++ii;
							}
							assert(ii - 1 < refinedClusters[r].matches.size());
							if (ii > jj && refinedClusters[r].matches[ii - 1].first.pos + smallOpts.globalK - 1 < *l) {
								vq.push_back(Cluster(jj, ii));
								jj = ii;
							}
							std::advance(l, 1);
						}	
						if (ii != v[i].end && l == qe) {
							assert(i < v.size());
							vq.push_back(Cluster(ii, v[i].end));
						}
					}
				}


				/*
				cerr << "length(seedSet2): " << length(seedSet2)<< endl;;
				for (TIterator tt = begin(seedSet2, seqan::Standard()); tt != end(seedSet2, seqan::Standard()); ++tt) {
					cerr << *tt << endl;
				}
				*/
				/*
				//cout << "vq.szie(): " << vq.size() << endl;
				const string filename2("/home/cmb-16/mjc/jingwenr/lra/lra_test/TEST/seedSet1.txt");  
				FILE *fm = fopen(filename2.c_str(), "w");
				SaveseedSet (seedSet2, fm); 
				fclose(fm);	
				*/		


				//----------------------------------------------------------------------------
				// Step 4: sort matches by y first and then by x coordinates in each Cluster vq[i](merged matches)  && split based on y boundaries
				//----------------------------------------------------------------------------
				for (unsigned int ij = 0; ij < vq.size(); ++ij) {
					int css = vq[ij].start; 
					int cee = vq[ij].start;
					cee = css + 1;
					GenomePos qStart = refinedClusters[r].matches[css].first.pos, 
					qEnd = refinedClusters[r].matches[css].first.pos + smallOpts.globalK, 
					tStart = refinedClusters[r].matches[css].second.pos, 
					tEnd = refinedClusters[r].matches[css].second.pos + smallOpts.globalK;	

					while(cee < vq[ij].end) {
						qStart = min(qStart, refinedClusters[r].matches[cee].first.pos);
						qEnd = max(qEnd, refinedClusters[r].matches[cee].first.pos + smallOpts.globalK);
						tStart = min(tStart, refinedClusters[r].matches[cee].second.pos);
						tEnd = max(tEnd, refinedClusters[r].matches[cee].second.pos + smallOpts.globalK);
						++cee;
					}
					vq[ij].qStart = qStart;
					vq[ij].qEnd = qEnd;
					vq[ij].tStart = tStart;
					vq[ij].tEnd = tEnd;

					//insert y boundaries into s2
					s2.insert(vq[ij].tStart);
					s2.insert(vq[ij].tEnd);

					// sort matches by y and then by x inside each vq[i]
					vector<pair<GenomeTuple, GenomeTuple>>::iterator begin = refinedClusters[r].matches.begin();
					if (vq[ij].start != 0) {
						std::advance(begin , vq[ij].start);
					}	 

					vector<pair<GenomeTuple, GenomeTuple>>::iterator end = refinedClusters[r].matches.begin();
					if (vq[ij].end != 0) {
						std::advance(end, vq[ij].end);
					}
					CartesianTargetSort<GenomeTuple>(begin, end);            
				}

				/*
				cerr << "finished sorting" << endl;
				// Debug code ----------- print out "seedSet"
				seqan::SeedSet<IndSeed, seqan::Unordered> seedSet2;
				for (unsigned int i = 0; i != vq.size(); ++i) {
					seqan::addSeed(seedSet2, IndSeed(vq[i].qStart, vq[i].tStart, vq[i].qEnd, vq[i].tEnd, i), seqan::Single());
				}
				
				seqan::Iterator<seqan::SeedSet<IndSeed, seqan::Unordered> >::Type it;
				for ( it = seqan::begin(seedSet2); it != seqan::end(seedSet2); ++it) {
					cerr << "beginPositionH(*it): " << beginPositionH(*it) << "  endPositionH(*it): " << endPositionH(*it)  << endl;
					cerr << "beginPositionV(*it): " << beginPositionV(*it) << "  endPositionV(*it): " << endPositionV(*it)  << endl;
					assert(beginPositionH(*it) <= endPositionH(*it) );
					assert(beginPositionV(*it) <= endPositionV(*it) );
				}	
				*/	


 				
 				v.clear();
				for (unsigned int i = 0; i < vq.size(); ++i) {
					std::set<unsigned int>::iterator ts, te; 
					ts = s2.lower_bound(vq[i].tStart);
					te = s2.lower_bound(vq[i].tEnd);
					std::set<unsigned int>::iterator l = ts;
					std::advance(l, 1);
					unsigned int ii = vq[i].start;
					unsigned int jj = ii;

					/*
					// debug code
					cout << "ts: " << *ts << endl;
					cout << "te: " << *te << endl;
					cout << "l: " << *l << endl;
					*/
					
					
					if (*l == *ts) {
						//debug code

						//cerr << "*l == *ts, push back" << "Cluster(" << vq[i].start << ", " << vq[i].end << ", " << vq[i].qStart << ", " << vq[i].qEnd << ", " << vq[i].tStart << ", " << vq[i].tEnd << ", 0)" << endl; 
						vt.push_back(Cluster(vq[i].start, vq[i].end));
					}
					else {
						while (ii != vq[i].end && l != te) { 
							while (refinedClusters[r].matches[ii].second.pos <= *l && ii != vq[i].end) {
								if (refinedClusters[r].matches[ii].second.pos + smallOpts.globalK - 1 >= *l) {
									// split 
									if (jj + 1 <= ii) {
										assert(jj < ii);
										vt.push_back(Cluster(jj, ii));
									}
									vt.push_back(Cluster(ii, ii + 1));
									jj = ii + 1;
								}
								++ii;
							}
							if (ii > jj && refinedClusters[r].matches[ii - 1].second.pos + smallOpts.globalK - 1 < *l) {
								assert(jj < ii);
								vt.push_back(Cluster(jj, ii));
								jj = ii;
							}
							std::advance(l, 1);
						}
						if (ii != vq[i].end && l == te) {
							assert(ii < vq[i].end);
							vt.push_back(Cluster(ii, vq[i].end));
						}
					}
				}
            

				vq.clear();
				for (unsigned int ij = 0; ij < vt.size(); ++ij) {
					int css = vt[ij].start; 
					int cee = vt[ij].start;
					cee = css + 1;
					GenomePos qStart = refinedClusters[r].matches[css].first.pos, 
					qEnd = refinedClusters[r].matches[css].first.pos + smallOpts.globalK, 
					tStart = refinedClusters[r].matches[css].second.pos, 
					tEnd = refinedClusters[r].matches[css].second.pos + smallOpts.globalK;	

					while(cee < vt[ij].end) {
						qStart = min(qStart, refinedClusters[r].matches[cee].first.pos);
						qEnd = max(qEnd, refinedClusters[r].matches[cee].first.pos + smallOpts.globalK);
						tStart = min(tStart, refinedClusters[r].matches[cee].second.pos);
						tEnd = max(tEnd, refinedClusters[r].matches[cee].second.pos + smallOpts.globalK);
						++cee;
					}
					vt[ij].qStart = qStart;
					vt[ij].qEnd = qEnd;
					vt[ij].tStart = tStart;
					vt[ij].tEnd = tEnd;
					assert(vt[ij].qStart < vt[ij].qEnd);
					assert(vt[ij].tStart < vt[ij].tEnd);
				}


				//-------------- debug code
				//cerr << "step 4 is finished" << endl;

				//----------------------------------------------------------------------------
				// Step 5:  Store the result in seedSet
				//----------------------------------------------------------------------------
				// (TODO: Jingwen)Not necessarily use seedSet as input of NaiveDp 
				seqan::clear(seedSet);
				for (unsigned int i = 0; i != vt.size(); ++i) {
					assert(vt[i].qStart < vt[i].qEnd);
					assert(vt[i].tStart < vt[i].tEnd);
					seqan::addSeed(seedSet, IndSeed(vt[i].qStart, vt[i].tStart, vt[i].qEnd, vt[i].tEnd, i), seqan::Single());
				}  
				/*
				// Debug code ----------- print out "seedSet"
				const string filename1("/home/cmb-16/mjc/jingwenr/lra/lra_test/TEST/MSseedSet" + std::to_string(r) + ".txt");  
				FILE *fd = fopen(filename1.c_str(), "w");
				SaveseedSet (seedSet, fd); 
				fclose(fd);	
				*/
					
		
			}

			else {
				continue;
			}	
		}
		else {
			for (int m=0; m< refinedClusters[r].matches.size(); m++) {
				seqan::addSeed(seedSet, seqan::Seed<IndexedSeed>(refinedClusters[r].matches[m].second.pos, refinedClusters[r].matches[m].first.pos, k), seqan::Single());
			}
		}
		
		/*
		//debug code
		cerr << "merge result" << endl;
		for (TIterator tt = begin(seedSet, seqan::Standard()); tt != end(seedSet, seqan::Standard()); ++tt) {
			cerr << *tt << endl;
		}
		*/

		/*
		seqan::Iterator<seqan::SeedSet<IndSeed, seqan::Unordered> >::Type it;
		for ( it = seqan::begin(seedSet); it != seqan::end(seedSet); ++it) {
			cerr << "beginPositionH(*it): " << beginPositionH(*it) << "  endPositionH(*it): " << endPositionH(*it)  << endl;
			cerr << "beginPositionV(*it): " << beginPositionV(*it) << "  endPositionV(*it): " << endPositionV(*it)  << endl;
			assert(beginPositionH(*it) <= endPositionH(*it) );
			assert(beginPositionV(*it) <= endPositionV(*it) );
		}		
		*/




		// Perform sparse chaining, uses time O(n log n).
		//
		// Merge anchors that are not overlapping

		//
		// The input is a set of seqan::Seed<seqan::Simple> > seeds
		// 
		seqan::String<IndSeed> chain;


		if (opts.NaiveDP) {
			// Perform sparse dp with convex gap cost
			seqan::clear(chain);
			if (seqan::length(seedSet) < 30000) {
				NaiveDP (seedSet, chain);
			}
			else {
				cerr << "Skipping naivedp on seed set of size " << seqan::length(seedSet) << endl;
			}

			/*
			//Debug code ------ print out "chain"
			/*
			if (r == 0) {
				//cout << "0" << endl;
				const string filename2("/home/cmb-16/mjc/jingwenr/lra/lra_test/TEST/NaiveDP." + std::to_string(r) + ".txt");  
				FILE *fz = fopen(filename2.c_str(), "w");
				SaveSparse (chain, fz);
				fclose(fz);		
			}		
			*/
		}
		else {
			seqan::clear(chain);
			if (seqan::length(seedSet) > 0) {
				seqan::chainSeedsGlobally(chain, seedSet, seqan::SparseChaining());
				cerr << "finished the chainSeedsGlobally" << endl;
			}

			/*
			//Debug code ------ print out "chain"
			const string filename4("/home/cmb-16/mjc/jingwenr/lra/lra_test/SeqanChaining_result.txt");  
			FILE *fi = fopen(filename4.c_str(), "w");
			SaveSparse (chain, fi);
			fclose(fi);	
			*/
		}
			
		if (opts.dotPlot) {
			stringstream outNameStrm;
			outNameStrm << baseName + "." << r << ".first-sdp.dots";
			ofstream baseDots(outNameStrm.str().c_str());
			for (int c=0; c < length(chain); c++) {
				int p =seqan::beginPositionH(chain[c]);
				baseDots << seqan::beginPositionH(chain[c]) << "\t" 
								 << seqan::beginPositionV(chain[c]) << "\t" 
								 << seqan::endPositionH(chain[c]) - seqan::beginPositionH(chain[c]) << "\t" << 7 << "\t0" << endl;				
			}
			baseDots.close();
		}

		vector<GenomePair> tupChain;
		int qPrev=0, tPrev=0;
		int csg = seqan::length(chain); 
		if (opts.mergeClusters) {
			//
			// Add small anchors to tupChain. (Use greedy algorithm to make small anchors overlap with each other)
			//
			/*
				tupChain.push_back(MatchingPos(GenomeTuple(0, beginPositionV(chain[ch])),
																			 GenomeTuple(0, beginPositionH(chain[ch])), length_of_split_kmer));
			*/

			for (int ch=0; ch < seqan::length(chain); ch++) { 

				//cerr << "ch: "<< ch<< endl;
				//cerr << "chain[" << ch <<"]   " << chain[ch] << endl; 

				unsigned int fprev = vt[chain[ch].index].start;
				unsigned int fcur = vt[chain[ch].index].start;
				tupChain.push_back(GenomePair(GenomeTuple(0, refinedClusters[r].matches[fcur].first.pos), 
																GenomeTuple(0, refinedClusters[r].matches[fcur].second.pos)));	
				
				//cerr << "tupChain[" << tupChain.size() << "]     push_back: " << "first: " << refinedClusters[r].matches[fcur].first.pos<<"  " << refinedClusters[r].matches[fcur].first.pos + smallOpts.globalK 
				//<< "   second: " << refinedClusters[r].matches[fcur].second.pos << " "<< refinedClusters[r].matches[fcur].second.pos + smallOpts.globalK - 1 <<endl;
				

				++fcur;					
				while (fcur < vt[chain[ch].index].end) {
					while (fcur < vt[chain[ch].index].end && (refinedClusters[r].matches[fcur].first.pos < refinedClusters[r].matches[fprev].first.pos + smallOpts.globalK - 1 ||
												refinedClusters[r].matches[fcur].second.pos < refinedClusters[r].matches[fprev].second.pos + smallOpts.globalK - 1)) {

						++fcur;

					}
					if (fcur != vt[chain[ch].index].end) {
						tupChain.push_back(GenomePair(GenomeTuple(0, refinedClusters[r].matches[fcur].first.pos), 
																GenomeTuple(0, refinedClusters[r].matches[fcur].second.pos)));	
				
						//cerr << "tupChain[" << tupChain.size() << "]     push_back: " << "first: " << refinedClusters[r].matches[fcur].first.pos<<"  " << refinedClusters[r].matches[fcur].first.pos + smallOpts.globalK 
						//<< "   second: " << refinedClusters[r].matches[fcur].second.pos << " "<< refinedClusters[r].matches[fcur].second.pos + smallOpts.globalK - 1 <<endl;
				
						assert(refinedClusters[r].matches[fcur].first.pos >= refinedClusters[r].matches[fprev].first.pos + smallOpts.globalK - 1);
						assert(refinedClusters[r].matches[fcur].second.pos >= refinedClusters[r].matches[fprev].second.pos + smallOpts.globalK - 1);								
					}


					fprev = fcur;
					++fcur;
				}
			}
			vt.clear();
			
			/*
			// debug code
			for (int ch = 0; ch < tupChain.size(); ++ch) {
				unsigned int fprev = 0;
				unsigned int fcur = 0;	
				++fcur;
				while(fcur < tupChain.size()) {
					assert(tupChain[fcur].first.pos >= tupChain[fprev].first.pos + smallOpts.globalK -1);
					assert(tupChain[fcur].second.pos >= tupChain[fprev].second.pos + smallOpts.globalK -1);
					fprev = fcur;
					++fcur;
				}		
			}
			*/





			/*
			// Debug code ----------- print out "seedSet"
			/*
			if (r == 0) {
				//cout << "1" << endl;
				const string filename4("/home/cmb-16/mjc/jingwenr/lra/lra_test/TEST/tupChain" + std::to_string(r) + ".txt");
				FILE *fir = fopen(filename4.c_str(), "w");
				SavetupChain (tupChain, fir, opts.globalK); 
				fclose(fir);
			}
			*/
			/*
			cerr << "put back seeds is finished " << endl;

			// debug code 
			for (int ch=0; ch < tupChain.size(); ++ ch) {
				cerr<<      "tuple chain["<< ch << "]: "
									  << tupChain[ch].first.pos << ", "
									  << tupChain[ch].first.pos + smallOpts.globalK - 1 << ", "
									  << tupChain[ch].second.pos << ", "
									  << tupChain[ch].second.pos + smallOpts.globalK - 1 << endl;
			
			}
			*/
				/*
			// Debug code ----------- print out "seedSet"
			//cout << "1" << endl;
			const string filename4("/home/cmb-16/mjc/jingwenr/lra/lra_test/TEST/tupChain" + std::to_string(r) + ".txt");
			FILE *fir = fopen(filename4.c_str(), "w");
			SavetupChain (tupChain, fir, opts.globalK); 
			fclose(fir);
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

		RemovePairedIndels(tupChain, chainClust, smallOpts);  
			
		seqan::clear(chain);
		for (int m=0; m< tupChain.size(); m++) {
			IndSeed s(tupChain[m].second.pos, tupChain[m].first.pos, glIndex.k - 1); // the reason of using  glIndex.k - 1 instead of glIndex.k here is that in later part we dont use "half open" seeds 
 			seqan::append(chain,s);

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
			assert(nextReadStart >= curReadEnd);
			GenomePos subreadLength = nextReadStart-curReadEnd;
			assert(nextGenomeStart >= curGenomeEnd);
			GenomePos subgenomeLength = nextGenomeStart-curGenomeEnd;

			if (nextReadStart > curReadEnd and nextGenomeStart > curGenomeEnd) {

				if (subreadLength > 50 and 
						subgenomeLength > 50 and
						opts.refineLevel & REF_DYN ) {

					GenomePos maxLen = max(subreadLength, subgenomeLength);
					if (maxLen < 500) {
						tinyOpts.globalK=5;
					}
					else if (maxLen < 2000) {
						tinyOpts.globalK=7;
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
						if (opts.NaiveDP) {
   							// Perform sparse dp with convex gap cost
							seqan::clear(gapChain);
							//NaiveDP< seqan::SeedSet<seqan::Seed<seqan::Simple>, seqan::Unordered>,seqan::Seed<seqan::Simple>(gapSeedSet, gapChain);
							NaiveDP(gapSeedSet, gapChain);
						}
						else {
							seqan::chainSeedsGlobally(gapChain, gapSeedSet, seqan::SparseChaining());
						}
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
			if (opts.dotPlot) {
				dotFile << seqan::beginPositionH(chain[c]) << "\t" 
								<< seqan::beginPositionV(chain[c]) << "\t" 
								<< seqan::endPositionH(chain[c]) - seqan::beginPositionH(chain[c]) << "\t1\t0" << endl;
			}

			alignment->blocks.push_back(Block(seqan::beginPositionV(chain[c]),
																				seqan::beginPositionH(chain[c]), glIndex.k - 1)); // use  glIndex.k - 1 instead of glIndex.k here, because we don't use half open seeds later
			if (alignment->blocks.size() > 1) {
				int last=alignment->blocks.size();
				assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].qPos);
				assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].tPos);
			}

			//			string curAnchor = string(genome.seqs[chromIndex], seqan::beginPositionV(chain[c]), glIndex.k );
			for (int cs = 0; cs < seqan::length(refinedChains[c]); cs++) {
				//
				// Refined anchors are with respect to the chained sequence
				nextRefinedReadStart   = seqan::beginPositionV(refinedChains[c][cs]);
				nextRefinedGenomeStart = seqan::beginPositionH(refinedChains[c][cs]);
				
				if (opts.dotPlot) {
					dotFile << seqan::beginPositionH(refinedChains[c][cs]) << "\t" 
									<< seqan::beginPositionV(refinedChains[c][cs]) << "\t" 
									<< seqan::endPositionH(refinedChains[c][cs]) - seqan::beginPositionH(refinedChains[c][cs]) << "\t2\t0" << endl;
				}

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
						int b;
						for (b=1; b < betweenAnchorAlignment.blocks.size(); b++) {
							assert(betweenAnchorAlignment.blocks[b-1].qPos + 
										 betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].qPos);
							assert(betweenAnchorAlignment.blocks[b-1].tPos + 
										 betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].tPos);						
						}

						betweenAnchorAlignment.blocks.clear();
					}
				}

				curRefinedReadEnd   = seqan::endPositionV(refinedChains[c][cs]);
				curRefinedGenomeEnd = seqan::endPositionH(refinedChains[c][cs]);
				alignment->blocks.push_back(Block(nextRefinedReadStart, nextRefinedGenomeStart, 
																					curRefinedReadEnd - nextRefinedReadStart));
				if (alignment->blocks.size() > 1) {
					int last=alignment->blocks.size();
					assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].qPos);
					assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].tPos);
				}

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
		if (opts.dotPlot) {
			dotFile.close();
		}

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
		pthread_mutex_unlock(semaphore);
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
