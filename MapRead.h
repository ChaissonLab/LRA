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
#include <iterator>
#include <ctime>
#include <cmath>	// std::log 

#include <sstream>
#include <thread>


#include "Clustering.h"
#include <thread>

#include "AffineOneGapAlign.h"
#include "GlobalChain.h"
#include "TupleOps.h"
#include "SparseDP.h"
#include "Merge.h"
#include "MergeAnchors.h"
#include "SparseDP_Forward.h"
#include "Chain.h"
#include "overload.h"


using namespace std;


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


int Matched(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te) {
	return min(qe-qs+1, te-ts+1); // TODO(Jingwen): check whether should add 1
}

void SetMatchAndGaps(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te, int &m, int &qg, int &tg) {
	m=Matched(qs, qe, ts, te);
	qg=qe-qs+1-m; //TODO(Jingwen): check whether should add 1
	tg=te-ts+1-m;
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

		while ( ovp < clusters.size() and clusters[a].Overlaps(clusters[ovp], 0.8 ) ) {
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

void SimpleMapQV(vector<SegAlignmentGroup> &alignments) {
	int a=0;
	int ovp=a;
	if (alignments.size() == 0) {
		return;
	}
	if (alignments.size() == 1) {
		alignments[a].mapqv=255;
		alignments[a].SetMapqv(); // give mapq 255 to every segment
		return;
	}
	while (a < alignments.size() - 1) {
		ovp=a+1;
		float num=1.0;
		float denom=1.0;

		while (ovp < alignments.size() and alignments[a].Overlaps(alignments[ovp], 0.9) ) {
			int nmDiff = alignments[a].nm - alignments[ovp].nm;
			if (nmDiff < 10) {
				denom += pow(0.5,nmDiff);
			}
			ovp++;
		}
		if (ovp == a+1){
			alignments[a].mapqv=255;
		}
		else {
			// comparing float ok because it is potentially not set.
			if (denom == 1.0) {
				alignments[a].mapqv=255;
			}
			else {
				alignments[a].mapqv=-10*log10(1-num/denom);
			}
		}
		alignments[a].SetMapqv();
		a=ovp;
	}			
}

int AlignSubstrings(char *qSeq, GenomePos &qStart, GenomePos &qEnd, char *tSeq, GenomePos &tStart, GenomePos &tEnd,
										vector<int> &scoreMat, vector<Arrow> &pathMat, Alignment &aln, Options &options) {
	
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

	int score = AffineOneGapAlign(readSeq, chromSeq, options.localMatch, options.localMismatch, options.localIndel, options.localBand, aln);
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
	bool operator()(const SegAlignmentGroup a, const SegAlignmentGroup b) const {
		return a.nm > b.nm;
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
	//cerr << "merging " << order.size() << " clusters" << endl;
	while(c < order.size()) {
		//cerr << "c: " << c << endl;
		//cerr << "Order[c]: " << order[c].qStart << "\t" << order[c].qEnd << "\t" << order[c].tStart << "\t" << order[c].tEnd << endl;
  		cn=c+1;
		int curEndChrom = genome.header.Find(order[c].tEnd);
		while (cn < order.size()) {
			//cerr <<"cn: " << cn<< endl;
			//cerr << "Order[cn]: " << order[cn].qStart << "\t" << order[cn].qEnd << "\t" << order[cn].tStart << "\t" << order[cn].tEnd << endl;
			int nextStartChrom = genome.header.Find(order[cn].tStart);
			int gap;
			//cerr << "(int)(order[cn].tStart - order[cn].qStart): " << (int)(order[cn].tStart - order[cn].qStart) << endl;
			//cerr << "(int)(order[c].tEnd-order[c].qEnd): " << (int)(order[c].tEnd-order[c].qEnd) << endl;
			gap = abs((int)((int)(order[cn].tStart - order[cn].qStart) - (int)(order[c].tEnd-order[c].qEnd)));
			//cerr << "gap: "<< gap << endl;
			/*
			for (int ci=0; ci < order[c].matches.size(); ci++) {
				cerr << c << "\t" 
						 << order[c].matches[ci].first.pos << "\t" 
						 << order[c].matches[ci].first.t << "\t" 
						 << order[c].matches[ci].second.pos << "\t"
						 << order[c].matches[ci].second.t << endl;
			}
			for (int ci=0; ci < order[cn].matches.size(); ci++) {
				cerr << cn << "\t" << order[cn].matches[ci].first.pos << "\t" 
						 << order[cn].matches[ci].first.t << "\t" 
						 << order[cn].matches[ci].second.pos << "\t"
						 << order[cn].matches[ci].second.t << endl;
			}*/
			// TODO(Jingwen): (gap < opts.maxDiag or order[c].Encompasses(order[cn],0.5)) or gap < opts.maxDiag???
			// (gap < opts.maxDiag or order[c].Encompasses(order[cn],0.5)) --> repetitive region will be merged into one cluster
			if (nextStartChrom == curEndChrom and order[c].strand == order[cn].strand and (gap < opts.maxDiag or order[c].Encompasses(order[cn],0.5))) {
				//cerr << "if happened " << endl;
				order[c].matches.insert(order[c].matches.end(), order[cn].matches.begin(), order[cn].matches.end());
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
				if (cn2 < order.size() and
						cn2 - cn < MAX_AHEAD and
						cn2 > cn and gap < opts.maxGap) {
					cn=cn2;
					order[c].matches.insert(order[c].matches.end(), order[cn].matches.begin(), order[cn].matches.end());
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

void RefineSubstrings(char *read, GenomePos readSubStart, GenomePos readSubEnd, char *genome, GenomePos genomeSubStart, GenomePos genomeSubEnd, 
											vector<int> &scoreMat, vector<Arrow> &pathMat, Alignment &aln, Options &opts) {

	aln.blocks.clear();
	AlignSubstrings(read, readSubStart, readSubEnd, genome, genomeSubStart, genomeSubEnd, scoreMat, pathMat, aln, opts);
	for (int b = 0; b < aln.blocks.size(); b++) {
		aln.blocks[b].qPos += readSubStart;
		aln.blocks[b].tPos += genomeSubStart;

	}
	
}

void SeparateMatchesByStrand(Read &read, Genome &genome, int k, vector<pair<GenomeTuple, GenomeTuple> > &allMatches,  vector<pair<GenomeTuple, GenomeTuple> > &forMatches,
								vector<pair<GenomeTuple, GenomeTuple> > &revMatches) {
	//
	// A value of 0 implies forward strand match.
	//
	vector<bool> strand(allMatches.size(), 0);
	int nForward=0;
	for (int i=0; i < allMatches.size(); i++) {
		int readPos = allMatches[i].first.pos;
		uint64_t refPos = allMatches[i].second.pos;
		char *genomePtr=genome.GlobalIndexToSeq(refPos);
		//
		// Read and genome are identical, the match is in the forward strand
		if (strncmp(&read.seq[readPos], genomePtr, k) == 0) {
			nForward++;
		}
		else {
			//
			// The k-mers are not identical, but a match was stored between
			// readPos and *genomePtr, therefore the match must be reverse.
			//
			strand[i] = true;
		}
	}
	//
	// Populate two lists, one for forward matches one for reverse.
	//
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


void traceback(vector<int> &clusterchain, int &i, vector<int> &clusters_predecessor, vector<bool> &used) {

	if (used[i] == 0) {
		clusterchain.push_back(i);	
		used[i] = 1;	
		if (clusters_predecessor[i] != -1) {
			if (used[clusters_predecessor[i]] == 0) {
				traceback(clusterchain, clusters_predecessor[i], clusters_predecessor, used);			
			}
			else {
				for (int lr = 0; lr < clusterchain.size(); ++lr) {
					used[clusterchain[lr]] = 0;
				}	
				clusterchain.clear();		
			}
		}
	}
}

// TODO(Jingwen): delete this function
void traceback (vector<int> &onechain, int &i, vector<int> &clusters_predecessor) {

	onechain.push_back(i);
	if (clusters_predecessor[i] != -1) {
		traceback(onechain, clusters_predecessor[i], clusters_predecessor);
	}
}


// This function removes spurious MERGED anchors in chain after 1st SDP
void RemoveSpuriousAnchors (vector<unsigned int> &chain, Options &opts, const vector<ClusterCoordinates> &Anchors, 
								const LogCluster &logcluster) {
	int cs = 0, ce = 1;
	vector<bool> remove(chain.size(), false);

	while (ce < chain.size()) {

		int dist = 0;
		int diffstrand = 0;
		if (Anchors[chain[ce-1]].strand == Anchors[chain[ce]].strand and Anchors[chain[ce-1]].strand == 0) { // forward stranded
			dist = max((int)Anchors[chain[ce]].qStart - Anchors[chain[ce-1]].qEnd, (int)Anchors[chain[ce]].tStart - Anchors[chain[ce-1]].tEnd);
		}
		else if (Anchors[chain[ce-1]].strand == Anchors[chain[ce]].strand and Anchors[chain[ce-1]].strand == 1) { // reverse stranded
			dist = max((int)Anchors[chain[ce]].qStart - Anchors[chain[ce-1]].qEnd, (int)Anchors[chain[ce-1]].tStart - Anchors[chain[ce]].tEnd);
		}
		else diffstrand = 1;

		//cerr << "ce-1: " << ce-1 << "  ce: " << ce << "  abs(dist): " << abs(dist) << "  diffstrand: " << diffstrand <<  endl;
		if (diffstrand == 0 and abs(dist) <= opts.maxRemoveSpuriousAnchorsDist) ce++;
		else {
			//cerr << "else" << endl;
			int anchorNum = 0, Length = 0;
			for (int i = cs; i < ce; i++) {
				anchorNum += Anchors[chain[i]].end - Anchors[chain[i]].start;
			}
			Length = min((int)Anchors[chain[ce - 1]].qEnd - Anchors[chain[cs]].qStart, (int)Anchors[chain[ce - 1]].tEnd - Anchors[chain[cs]].tStart);
			int coarse = Anchors[chain[ce-1]].coarseSubCluster; // Assume chain[cs]... chain[ce-1] are in the same SubCluster
			if ((anchorNum < opts.minRemoveSpuriousAnchorsNum or abs(Length) < opts.minRemoveSpuriousAnchorsLength) and 
				       (float)anchorNum/(logcluster.SubCluster[coarse].end - logcluster.SubCluster[coarse].start) < 0.05) { 
				//
				// The third condition enables small inversion or small stretches of matches 

				//cerr << "cs: " << cs << "  ce: " << ce << endl;
				//cerr << "anchorNum: " << anchorNum << " Length: " << Length << endl;
				//cerr << "ce - cs < min" << endl;
				for (int i = cs; i < ce; i++) {
					remove[i] = true;
				}
			}
			cs = ce;
			ce++;
		}
		//cerr << "cs: " << cs << "  ce: " << ce << endl;
	}
	if (ce == chain.size() and cs < chain.size()) {
		int anchorNum = 0, Length = 0;
		for (int i = cs; i < ce; i++) {
			anchorNum += Anchors[chain[i]].end - Anchors[chain[i]].start;
		}
		Length = min((int)Anchors[chain[ce - 1]].qEnd - Anchors[chain[cs]].qStart, (int)Anchors[chain[ce - 1]].tEnd - Anchors[chain[cs]].tStart);
		int coarse = Anchors[chain[ce-1]].coarseSubCluster;
		if ((anchorNum < opts.minRemoveSpuriousAnchorsNum or abs(Length) < opts.minRemoveSpuriousAnchorsLength) and 
			 (float)anchorNum/(logcluster.SubCluster[coarse].end - logcluster.SubCluster[coarse].start) < 0.05) {
			for (int i = cs; i < ce; i++) {
				remove[i] = true;
			}
		}		
	}

	int m = 0;
	for (int s = 0; s < chain.size(); s++) {
		if (remove[s] == false) {
			chain[m] = chain[s];
			m++;
		}
	}
	chain.resize(m);
}



// This function removes spurious UNMERGED anchors in chain after 1st SDP
void RemoveSpuriousAnchors (vector<unsigned int> &chain, Options &opts, const Cluster &Anchors, const LogCluster &logcluster) {
	int cs = 0, ce = 1;
	vector<bool> remove(chain.size(), false);

	while (ce < chain.size()) {

		int dist = 0;
		int diffstrand = 0;
		if (Anchors.strands[chain[ce-1]] == Anchors.strands[chain[ce]] and Anchors.strands[chain[ce-1]] == 0) { // forward stranded
			dist = max((int)Anchors.matches[chain[ce]].first.pos - Anchors.matches[chain[ce-1]].first.pos - opts.globalK, 
					   (int)Anchors.matches[chain[ce]].second.pos - Anchors.matches[chain[ce-1]].second.pos - opts.globalK);
		}
		else if (Anchors.strands[chain[ce-1]] == Anchors.strands[chain[ce]] and Anchors.strands[chain[ce-1]] == 1) { // reverse stranded
			dist = max((int)Anchors.matches[chain[ce]].first.pos - Anchors.matches[chain[ce-1]].first.pos - opts.globalK, 
				       (int)Anchors.matches[chain[ce-1]].second.pos - Anchors.matches[chain[ce]].second.pos - opts.globalK);
		}
		else diffstrand = 1;

		//cerr << "ce-1: " << ce-1 << "  ce: " << ce << "  dist: " << dist << endl;
		if (diffstrand == 0 and abs(dist) <= opts.maxRemoveSpuriousAnchorsDist) ce++;
		else {
			int coarse = Anchors.coarseSubCluster[chain[ce-1]];
			if (ce - cs < opts.minRemoveSpuriousAnchorsNum and 
					(float)(ce - cs)/(logcluster.SubCluster[coarse].end - logcluster.SubCluster[coarse].start) < 0.05) {
				//cerr << "cs: " << cs << "  ce: " << ce << endl;
				//cerr << "anchorNum: " << anchorNum << " Length: " << Length << endl;
				//cerr << "ce - cs < min" << endl;
				for (int i = cs; i < ce; i++) {
					remove[i] = true;
				}
			}
			cs = ce;
			ce++;
		}
	}
	if (ce == chain.size() and cs < chain.size()) {
		int coarse = Anchors.coarseSubCluster[chain[ce-1]];
		if (ce - cs < opts.minRemoveSpuriousAnchorsNum and 
				(float)(ce - cs)/(logcluster.SubCluster[coarse].end - logcluster.SubCluster[coarse].start) < 0.05) {
			for (int i = cs; i < ce; i++) {
				remove[i] = true;
			}
		}		
	}

	int m = 0;
	for (int s = 0; s < chain.size(); s++) {
		if (remove[s] == false) {
			chain[m] = chain[s];
			m++;
		}
	}
	chain.resize(m);
}


// This function removes spurious UNMERGED anchors in chain after 2nd SDP
void RemoveSpuriousAnchors(vector<unsigned int> &chain, Options &opts, const GenomePairs &Anchors) {
	int cs = 0, ce = 1;
	vector<bool> remove(chain.size(), false);

	while (ce < chain.size()) {

		int dist = 0;
		dist = max((int)Anchors[chain[ce]].first.pos - Anchors[chain[ce-1]].first.pos - opts.globalK, 
				   (int)Anchors[chain[ce]].second.pos - Anchors[chain[ce-1]].second.pos - opts.globalK);

		//cerr << "ce-1: " << ce-1 << "  ce: " << ce << "  dist: " << dist << endl;
		if (abs(dist) <= opts.maxRemoveSpuriousAnchorsDist) ce++;
		else {
			if (ce - cs < opts.minRemoveSpuriousAnchorsNum) {
				//cerr << "cs: " << cs << "  ce: " << ce << endl;
				//cerr << "anchorNum: " << anchorNum << " Length: " << Length << endl;
				//cerr << "ce - cs < min" << endl;
				for (int i = cs; i < ce; i++) {
					remove[i] = true;
				}
			}
			cs = ce;
			ce++;
		}
	}
	if (ce == chain.size() and cs < chain.size()) {
		if (ce - cs < opts.minRemoveSpuriousAnchorsNum) {
			for (int i = cs; i < ce; i++) {
				remove[i] = true;
			}
		}		
	}

	int m = 0;
	for (int s = 0; s < chain.size(); s++) {
		if (remove[s] == false) {
			chain[m] = chain[s];
			m++;
		}
	}
	chain.resize(m);	
}

int MapRead(const vector<float> & LookUpTable, Read &read, Genome &genome, vector<GenomeTuple> &genomemm, LocalIndex &glIndex, Options &opts, ostream *output, pthread_mutex_t *semaphore=NULL) {
	
	string baseName = read.name;

	for (int i=0; i < baseName.size(); i++) {	if (baseName[i] == '/') baseName[i] = '_';	}


	vector<GenomeTuple> readmm; // readmm stores minimizers
	vector<pair<GenomeTuple, GenomeTuple> > allMatches, forMatches, revMatches, matches;

	if (opts.storeAll) {
		Options allOpts = opts;
		allOpts.globalW=1;
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, allOpts.globalK, allOpts.globalW, readmm);			
	}
	else {
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, opts.globalK, opts.globalW, readmm);
	}
	sort(readmm.begin(), readmm.end()); //sort kmers in readmm(minimizers)
	//
	// Add matches between the read and the genome.
	//
	CompareLists(readmm, genomemm, allMatches, opts);
	DiagonalSort<GenomeTuple>(allMatches); // sort fragments in allMatches by forward diagonal, then by first.pos(read)

	// TODO(Jinwen): delete this after debug
	if (opts.dotPlot) {
		ofstream clust("all-matches.dots");
		for (int m=0; m < allMatches.size(); m++) {
			clust << allMatches[m].first.pos << "\t" << allMatches[m].second.pos << "\t" << allMatches[m].first.pos+ opts.globalK << "\t" 
				<< allMatches[m].second.pos + opts.globalK << "\t0\t0"<<endl;
		}
		clust.close();
	}


	SeparateMatchesByStrand(read, genome, opts.globalK, allMatches, forMatches, revMatches);


	// TODO(Jingwen): only for debug, delete later
	if (opts.dotPlot) {
		ofstream clust("for-matches0.dots");
		for (int m=0; m < forMatches.size(); m++) {
			clust << forMatches[m].first.pos << "\t" << forMatches[m].second.pos << "\t" << opts.globalK + forMatches[m].first.pos << "\t"
					<< forMatches[m].second.pos + opts.globalK << "\t0\t0"<<endl;
		}
		clust.close();
		ofstream rclust("rev-matches0.dots");
		for (int m=0; m < revMatches.size(); m++) {			
			rclust << revMatches[m].first.pos << "\t" << revMatches[m].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[m].first.pos << "\t"
					<< revMatches[m].second.pos << "\t0\t0"<<endl;
		}
		rclust.close();
	}


	int minDiagCluster = 0; // This parameter will be set inside function CleanOffDiagonal, according to anchors density
	CleanOffDiagonal(forMatches, opts, minDiagCluster);
	vector<Cluster> clusters;
	vector<Cluster> roughclusters;
	int forwardStrand=0;
	// maxDiag must be large enough for the following function "StoreDiagonalClusters". 
	// Otherwise fragments on the same line (line with little curve) might end up in several clusters[i], instead of one clusters[i]

	// The strategy we are taking: 
	// first let maxDiag be a small number but not too small(like 500), which also alleviate the cases where anchors are on a curvy line.
	// Then break the clusters[i] into two if any two anchors are father than maxGap. 
	StoreDiagonalClusters(forMatches, roughclusters, opts, 0, forMatches.size(), true, false, minDiagCluster, forwardStrand); // rough == true means only storing "start and end" in every clusters[i]
	for (int c = 0; c < roughclusters.size(); c++) {
		CartesianSort(forMatches, roughclusters[c].start, roughclusters[c].end);
		StoreDiagonalClusters(forMatches, clusters, opts, roughclusters[c].start, roughclusters[c].end, false, false, minDiagCluster, forwardStrand);
	}


	AntiDiagonalSort<GenomeTuple>(revMatches, genome.GetSize());
	CleanOffDiagonal(revMatches, opts, minDiagCluster, 1);

	// TODO(JIngwen): Only for debug
	if (opts.dotPlot) {
		ofstream rclust("rev-matches1.dots");
		for (int m=0; m < revMatches.size(); m++) {			
			rclust << revMatches[m].first.pos << "\t" << revMatches[m].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[m].first.pos << "\t"
					<< revMatches[m].second.pos << "\t0\t0"<<endl;
		}
		rclust.close();
	}

	vector<Cluster> revroughClusters;
	int reverseStrand=1;
	StoreDiagonalClusters(revMatches, revroughClusters, opts, 0, revMatches.size(), true, false, minDiagCluster, reverseStrand);

	for (int c = 0; c < revroughClusters.size(); c++) {
		CartesianSort(revMatches, revroughClusters[c].start, revroughClusters[c].end);
		StoreDiagonalClusters(revMatches, clusters, opts, revroughClusters[c].start, revroughClusters[c].end, false, false, minDiagCluster, reverseStrand);
	}

	if (opts.dotPlot) {
		ofstream clust("for-matches.dots");
		for (int m=0; m < forMatches.size(); m++) {
			clust << forMatches[m].first.pos << "\t" << forMatches[m].second.pos << "\t" << opts.globalK + forMatches[m].first.pos << "\t"
					<< forMatches[m].second.pos + opts.globalK << "\t0\t0"<<endl;
		}
		clust.close();
		ofstream rclust("rev-matches.dots");
		for (int m=0; m < revMatches.size(); m++) {			
			rclust << revMatches[m].first.pos << "\t" << revMatches[m].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[m].first.pos << "\t"
					<< revMatches[m].second.pos << "\t0\t0"<<endl;
		}
		rclust.close();

		ofstream wclust("roughclusters-matches.dots");
		for (int m=0; m < roughclusters.size(); m++) {
			for (int c = roughclusters[m].start; c < roughclusters[m].end; ++c) {
				wclust << forMatches[c].first.pos << "\t" << forMatches[c].second.pos << "\t" << opts.globalK + forMatches[c].first.pos << "\t"
					<< forMatches[c].second.pos + opts.globalK << "\t" << m << "\t0"<<endl;				
			}

		}

		ofstream revclust("revroughClusters-matches.dots");
		for (int m=0; m < revroughClusters.size(); m++) {
			for (int c = revroughClusters[m].start; c < revroughClusters[m].end; ++c) {
				revclust << revMatches[c].first.pos << "\t" << revMatches[c].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[c].first.pos << "\t"
					 << revMatches[c].second.pos<< "\t" << m << "\t0"<<endl;				
			}

		}
		revclust.close();
	}





	//
	// Add pointers to seq that make code more readable.
	//
	char *readRC;

	CreateRC(read.seq, read.length, readRC);
	char *strands[2] = { read.seq, readRC };


	if (clusters.size() != 0) {

		if (opts.dotPlot) {
			ofstream clust("clusters-pre-merge.tab");
			for (int c =0; c < clusters.size(); c++) {

				for (int m=0; m < clusters[c].matches.size(); m++) {
					clust << clusters[c].matches[m].first.pos << "\t" 
						  << clusters[c].matches[m].second.pos << "\t" 
						  << clusters[c].matches[m].first.pos + opts.globalK << "\t"
						  << clusters[c].matches[m].second.pos + opts.globalK << "\t" 
						  << c << "\t" 
						  << clusters[c].strand << endl;
				}

				//clust << clusters[c].qStart << "\t" << clusters[c].tStart << "\t" << clusters[c].qEnd << "\t" << clusters[c].tEnd << "\t" << c  << "\t" << clusters[c].strand << endl;
			}
			clust.close();
		}


		ClusterOrder clusterOrder(&clusters);  // clusterOrder is sorted first by tStart, then by qStart
		//MergeAdjacentClusters(clusterOrder, genome, opts);


		//
		// Merge clusters
		//
		vector<float> clusters_value(clusters.size(), 0);
		vector<float> clusters_value_ford(clusters.size(), 0);
		vector<float> clusters_value_rev(clusters.size(), 0);

		vector<int> clusters_predecessor(clusters.size(), -1);
		vector<int> clusters_predecessor_ford(clusters.size(), -1);
		vector<int> clusters_predecessor_rev(clusters.size(), -1);

		//clusters_value[0] = (float)(max(clusterOrder[0].qEnd - clusterOrder[0].qStart, clusterOrder[0].tEnd - clusterOrder[0].tStart));
		//vector<int> repetitivenum(clusters.size(), 0);

		// Get clusters_value 
		for (int cv = 0; cv < clusters_value.size(); cv++) {
			clusters_value[cv] = (float)(max(clusterOrder[cv].qEnd - clusterOrder[cv].qStart, clusterOrder[cv].tEnd - clusterOrder[cv].tStart));
			clusters_value_ford[cv] = clusters_value[cv];
			clusters_value_rev[cv] = clusters_value[cv];
		}

		if (clusters.size() > 1) {
			for (int c = 1; c < clusters.size(); ++c) {

				for (int s = 0; s <= c - 1; ++s) {

					int gap = 0;
					int ReverseStrand = -1; // ReverseStrand depicts the direction of clustersOrder[c] and clustersOrder[s]

					if (clusterOrder[c].strand == clusterOrder[s].strand and clusterOrder[c].strand == 1) {
						ReverseStrand = -1;
						gap = abs((int)(clusterOrder[c].qEnd + clusterOrder[c].tStart) - (int)(clusterOrder[s].qStart + clusterOrder[s].tEnd));
					}
					else {
						ReverseStrand = 1;
						gap = abs((int)(clusterOrder[c].qStart - clusterOrder[c].tStart) - (int)(clusterOrder[s].qEnd - clusterOrder[s].tEnd));
					}

					if ( ((clusterOrder[s].qEnd - ((clusterOrder[s].qEnd - clusterOrder[s].qStart)/3)*2) < clusterOrder[c].qStart or ReverseStrand == -1)
							and  ((clusterOrder[s].tEnd - ((clusterOrder[s].tEnd - clusterOrder[s].tStart)/3)*2) < clusterOrder[c].tStart or ReverseStrand == -1)
							and ((!clusterOrder[c].Encompasses(clusterOrder[s], 0.7) and !clusterOrder[s].Encompasses(clusterOrder[c], 0.4)) or 
								(!clusterOrder[c].Encompasses(clusterOrder[s], 0.4) and !clusterOrder[s].Encompasses(clusterOrder[c], 0.7)))
							and ((clusterOrder[c].qEnd - ((clusterOrder[c].qEnd - clusterOrder[c].qStart)/3)*2) < clusterOrder[s].qStart or ReverseStrand == 1)
							and ((clusterOrder[c].tEnd - ((clusterOrder[c].tEnd - clusterOrder[c].tStart)/3)*2) > clusterOrder[s].tStart or ReverseStrand == 1)
							//and clusterOrder[c].strand == clusterOrder[s].strand and gap < opts.maxGap
							) {

						float objective;
						float rate = max(clusterOrder[c].OverlapsRate(clusterOrder[s]), clusterOrder[s].OverlapsRate(clusterOrder[c]));
						int overlap = clusterOrder[c].Overlaps(clusterOrder[s]);
						//objective = (float)(max(clusterOrder[c].qEnd - clusterOrder[c].qStart, clusterOrder[c].tEnd - clusterOrder[c].tStart)) +
							//			clusters_value[s] - 0.5*((float)gap) - 2*rate*((float)overlap);

						//objective = clusters_value[c] + clusters_value[s] - log((float)gap) - 0.5*((float)gap);
						//objective = clusters_value[c] + clusters_value[s] - LookUpTable[(int)floor((gap-501)/5)] - 0.5*((float)gap);

						if (ReverseStrand == -1) {
							objective = (float)(max(clusterOrder[c].qEnd - clusterOrder[c].qStart, clusterOrder[c].tEnd - clusterOrder[c].tStart)) +
										clusters_value_rev[s] - 0.2*((float)gap) - 2*rate*((float)overlap);
							if (objective >= clusters_value[c]) {
								clusters_value_rev[c] = objective;
								clusters_predecessor_rev[c] = s;
							}
						}
						else if (ReverseStrand == 1) {
							objective = (float)(max(clusterOrder[c].qEnd - clusterOrder[c].qStart, clusterOrder[c].tEnd - clusterOrder[c].tStart)) +
										clusters_value_ford[s] - 0.2*((float)gap) - 2*rate*((float)overlap);
							if (objective >= clusters_value[c]) {
								clusters_value_ford[c] = objective;
								clusters_predecessor_ford[c] = s;
							}
						}

						/*if (objective >= clusters_value[c]) {
							clusters_value[c] = objective;
							clusters_predecessor[c] = s;
						}	*/	
					}
				} 
			}	
		}

		for (int r = 0; r < clusters_value.size(); ++r) {
			clusters_value[r] = max(clusters_value_ford[r], clusters_value_rev[r]);

			if (clusters_value_ford[r] >= clusters_value_rev[r]) {
				clusters_predecessor[r] = 0;
			}
			else clusters_predecessor[r] = 1;
		}


		Clusters_valueOrder clusters_valueOrder(&clusters_value); // sort clusters_value in descending order
		//
		// traceback chains until chains overlap less than 30% on the read
		// chains that are overlapped over 70% are considered as secondary alignments
		// chains that are overlapped less than 50% are considered as different primary alignments
		//
		vector<Primary_chain> Primary_chains;
		Primary_chains.clear();
		float maxValue = clusters_valueOrder[0];
		vector<bool> used(clusters.size(), 0); // used[i] == 1 means clusterOrder[i] has been used

		for (int r = 0; r < clusters_value.size() ; ++r) {

			vector<int> onechain;
			int idx = clusters_valueOrder.index[r];
			if (clusters_predecessor[idx] == 0) {
				traceback(onechain, idx, clusters_predecessor_ford, used);
			}
			else traceback(onechain, idx, clusters_predecessor_rev, used);

			if (onechain.size() != 0) {
				// Note: onechain store index from the last one to the first one
				int TrueIndex = clusterOrder.index[onechain[0]];
				GenomePos qEnd = clusters[TrueIndex].qEnd;
				GenomePos tEnd = clusters[TrueIndex].tEnd;

				TrueIndex = clusterOrder.index[onechain.back()];
				GenomePos qStart = clusters[TrueIndex].qStart;
				GenomePos tStart = clusters[TrueIndex].tStart;

				//
				// If this chain overlap with read greater than 30%, insert it to clusterchain
				if (((float)(qEnd - qStart)/read.length) > 0.3) {
					//
					// Compare onechain to all the Primary_chains we've found. 
					// If onechain overlaps with one primary chain over 70% ---> onechain is a secondary chain 
					// If onechain overlaps with all the primary chains less than 50% ---> onechain is another primary chain
					if (Primary_chains.size() == 0) {
						Primary_chains.push_back(Primary_chain(qStart, qEnd, tStart, tEnd, onechain));
					} 
					else {
						bool newpr = 1, inserted = 0;
						int p = 0;
						while (p < Primary_chains.size()) {
							if (Primary_chains[p].Overlaps(qStart, qEnd, 0.7)) {
								Primary_chains[p].chains.push_back(onechain);
								inserted = 1;
								break;
							}
							else if (Primary_chains[p].Overlaps(qStart, qEnd, 0.5)) {
								newpr = 0;
							}
							++p;
						}			
						if (p == Primary_chains.size() - 1 and inserted == 0 and newpr == 1) {
							Primary_chains.push_back(Primary_chain(qStart, qEnd, tStart, tEnd, onechain));
						}	
					}

				}
				else break;				
			}
		}

		clusters_value.clear();

		// change the index in Primary_chains[i].chains[j] to the index refering to clusters
		int chainNum = 0;
		for (int r = 0; r < Primary_chains.size(); r++) {
			for (int p = 0; p < Primary_chains[r].chains.size(); p++) {
				chainNum++;
				for (int l = 0; l < Primary_chains[r].chains[p].size(); l++) {
					int id = Primary_chains[r].chains[p][l]; 
					Primary_chains[r].chains[p][l] = clusterOrder.index[id];
				}
			}
		}
/*
		// Output the primary chains
		if (opts.dotPlot {
			for (int c = 0; c < Primary_chains.size(); c++) {
				stringstream outNameStrm;
				outNameStrm << "clusters-sdp." << c << ".dots";
				ofstream baseDots(outNameStrm.str().c_str());
				for (int s = 0; s < Primary_chains[c].chains[0].size(); ++s) {
					for (int m = 0; m < clusters[Primary_chains[c].chains[0][s]].matches.size(); ++m) {
						baseDots << clusters[Primary_chains[c].chains[0][s]].matches[m].first.pos << "\t" 
								 << clusters[Primary_chains[c].chains[0][s]].matches[m].second.pos << "\t" 
								 << clusters[Primary_chains[c].chains[0][s]].matches[m].first.pos + opts.globalK << "\t"
								 << clusters[Primary_chains[c].chains[0][s]].matches[m].second.pos + opts.globalK << "\t"
								 << s << "\t" 
								 << clusters[Primary_chains[c].chains[0][s]].strand << endl;					
					}				
				}
				baseDots.close();
			}
		}
*/

		// TODO(Jingwen): only for debug and delete this later. Output the 1st every chain subject to the 1st primary chain
		if (opts.dotPlot) {
			for (int c = 0; c < Primary_chains.size(); c++) {

				for (int ss = 0; ss < Primary_chains[c].chains.size(); ++ss) {
					stringstream outNameStrm;
					outNameStrm << "clusters-dp." << ss << ".dots";
					ofstream baseDots(outNameStrm.str().c_str());
					for (int s = 0; s < Primary_chains[c].chains[ss].size(); ++s) {
						for (int m = 0; m < clusters[Primary_chains[c].chains[ss][s]].matches.size(); ++m) {
							baseDots << clusters[Primary_chains[c].chains[ss][s]].matches[m].first.pos << "\t" 
									 << clusters[Primary_chains[c].chains[ss][s]].matches[m].second.pos << "\t" 
									 << clusters[Primary_chains[c].chains[ss][s]].matches[m].first.pos + opts.globalK << "\t"
									 << clusters[Primary_chains[c].chains[ss][s]].matches[m].second.pos + opts.globalK << "\t"
									 << s << "\t" 
									 << clusters[Primary_chains[c].chains[ss][s]].strand << endl;					
						}				
					}					
					baseDots.close();
				}
			}
		}


		//
		// Record every chain's information in logClusters
		//
		vector<LogCluster> logClusters(chainNum);
		vector<int> ind(clusters.size(), 0); // ind[i] == 1 means clusters[ind[i]] should be kept
		int num = 0;
		for (int c = 0; c < Primary_chains.size(); ++c) {

			int pr = num;
			for (int  p = 0; p < Primary_chains[c].chains.size(); ++p) {
				
				int id = Primary_chains[c].chains[p][0];
				ind[id] = 1;

				if (p != 0) {
					logClusters[num].ISsecondary = 1;
					logClusters[num].primary = pr;
					logClusters[pr].secondary.push_back(num);
				}
				if (clusters[id].strand == 1) {
					logClusters[num].SubCluster.push_back(Cluster(0, clusters[id].matches.size(), read.length - clusters[id].qEnd, 
																	read.length - clusters[id].qStart, clusters[id].tStart, 
																		clusters[id].tEnd, clusters[id].strand, id));	
					// Swap rev anchors	
					for (int m = 0; m < clusters[id].matches.size(); m++) {
						clusters[id].matches[m].first.pos = read.length - (clusters[id].matches[m].first.pos + opts.globalK);
					}						
				}
				else {
					logClusters[num].SubCluster.push_back(Cluster(0, clusters[id].matches.size(), clusters[id].qStart, clusters[id].qEnd,
																clusters[id].tStart, clusters[id].tEnd, clusters[id].strand, id));
				}

				if (Primary_chains[c].chains[p].size() > 1) {

					for (int s = 1; s < Primary_chains[c].chains[p].size(); ++s) {

						int lc = clusters[id].matches.size();
						int ids = Primary_chains[c].chains[p][s];
						clusters[id].matches.insert(clusters[id].matches.end(), clusters[ids].matches.begin(), clusters[ids].matches.end());
						clusters[id].qStart = clusters[ids].qStart;
						clusters[id].tStart = clusters[ids].tStart;

						// If this is a reverse strand, then swap it, which will make the later step "refine clusters" easier
						// and also swap the qStart, qEnd
						if (clusters[ids].strand == 1) {
							logClusters[num].SubCluster.push_back(Cluster(lc, clusters[id].matches.size(), read.length - clusters[ids].qEnd,
																		 read.length - clusters[ids].qStart, clusters[ids].tStart, 
																		 clusters[ids].tEnd, clusters[ids].strand));

							for (int m = lc; m < clusters[id].matches.size(); m++) {
								clusters[id].matches[m].first.pos = read.length - (clusters[id].matches[m].first.pos + opts.globalK);
							}				
						}
						else {
							logClusters[num].SubCluster.push_back(Cluster(lc, clusters[id].matches.size(), clusters[ids].qStart, clusters[ids].qEnd,
												clusters[ids].tStart, clusters[ids].tEnd, clusters[ids].strand));
						}
					}			
				}
				logClusters[num].SetCoarse();	
				++num;			
			}
		}		

		vector<int> clustersempty = ind;
		if (clusters.size() > 0) {
			for (int c = 1; c < clusters.size(); ++c) {
				clustersempty[c] += clustersempty[c-1];
			}				
		}

		for (int c = 0; c < logClusters.size(); ++c) {
			logClusters[c].coarse = clustersempty[logClusters[c].coarse] - 1;
		}

		/*
		// delete clusters with small size
		for (int c= 0; c < clusterchain.size(); ++c) {
			if (((float)(clusters[clusterchain[c][0]].qEnd - clusters[clusterchain[c][0]].qStart))/read.length < 0.3) ind[clusterchain[c][0]] = 0;
		}
		*/
		

		int cl=0;
		int cn;
		for(cl=0, cn=0; cn < clusters.size(); cn++) {
			if (ind[cn] == 1) {
				clusters[cl] = clusters[cn];
				cl++;
			}
		}


		clusters.resize(cl);
		clusters_predecessor.clear();
		Primary_chains.clear();

		/* TODO(Jingwen): I think we don't need Cartesian sort clusters before RemoveOverlappingClusters. 
		Just sorting them first by forward diag and then by read is okay
		/////////////TODO(Jingwen): is the following redundant?	Can't just Cartesian sort clusters?
		ClusterOrder reducedOrder(&clusters);
		vector<Cluster> reducedClusters;
		for (int i = 0; i < reducedOrder.size(); i++) {
			reducedClusters.push_back(reducedOrder[i]);
		}
		clusters= reducedClusters; // Now clusters is Cartesian sorted(first by strand, then by tStart, then by qStart)
		//////////////////
		*/
		
		//TODO(Jingwen): check whether we need the following sort
		//sort(clusters.begin(), clusters.end(), OrderClusterBySize()); // clusters are sorted in descending order
		

		if (opts.dotPlot) {
			ofstream matchfile("long_matches.tab");
			for (int m =0; m < matches.size(); m++) {
				matchfile << matches[m].first.pos << "\t" << matches[m].second.pos << "\t" << opts.globalK << "\t0\t0" << endl;
			}
			matchfile.close();
			ofstream clust("clusters.tab");
			for (int c =0; c < logClusters.size(); c++) {
		
				for (int n=0; n<logClusters[c].SubCluster.size();n++) {

					for (int m=logClusters[c].SubCluster[n].start; m<logClusters[c].SubCluster[n].end; m++) {

						if (logClusters[c].SubCluster[n].strand == 0) {
							clust << clusters[logClusters[c].coarse].matches[m].first.pos << "\t"
								  << clusters[logClusters[c].coarse].matches[m].second.pos << "\t" 
								  << clusters[logClusters[c].coarse].matches[m].first.pos + opts.globalK << "\t" 
								  << clusters[logClusters[c].coarse].matches[m].second.pos + opts.globalK  << "\t" 
								 << c << "\t" 
								 << logClusters[c].SubCluster[n].strand << endl;
						}
						else {
							clust 	<< read.length - clusters[logClusters[c].coarse].matches[m].first.pos - 1 << "\t" 
									<< clusters[logClusters[c].coarse].matches[m].second.pos << "\t" 
									<< read.length - clusters[logClusters[c].coarse].matches[m].first.pos - 1 - opts.globalK << "\t"
									<< clusters[logClusters[c].coarse].matches[m].second.pos + opts.globalK << "\t"
									<<  c << "\t" 
									<< logClusters[c].SubCluster[n].strand << endl;
						}
					}
				}
			}
			clust.close();	
		}

		//RemoveOverlappingClusters(clusters, opts); //TODO(Jingwen): check whether need to keep this

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

		// Set the parameters for merging anchors and 1st SDP
		Options smallOpts = opts;
		smallOpts.globalK=glIndex.k;
		smallOpts.globalW=glIndex.w;
		smallOpts.globalMaxFreq=6;
		smallOpts.cleanMaxDiag=10;// used to be 25
		smallOpts.maxDiag=50;
		smallOpts.maxGapBtwnAnchors=100; // used to be 200 // 200 seems a little bit large
		smallOpts.minDiagCluster=3;

		Options tinyOpts = smallOpts;
		tinyOpts.globalMaxFreq=3;
		tinyOpts.maxDiag=5;
		tinyOpts.minDiagCluster=2;
		tinyOpts.globalK=smallOpts.globalK-3;
		tinyOpts.minRemoveSpuriousAnchorsNum=5;
		tinyOpts.maxRemoveSpuriousAnchorsDist=50;


		vector<Cluster> refinedClusters(clusters.size());
		vector<LogCluster> refinedLogClusters(clusters.size());

		//logClusters[c] stores all anchors from one chain
		for (int c = 0; c < logClusters.size(); c++) {
			
			if (clusters[logClusters[c].coarse].matches.size() == 0) {
				continue;
			}			

			if (opts.dotPlot) {
				stringstream outNameStrm;
				outNameStrm << baseName + "." << c << ".clust.dots";
				ofstream baseDots(outNameStrm.str().c_str());
				for (int n=0; n<logClusters[c].SubCluster.size();n++) {

					for (int m=logClusters[c].SubCluster[n].start; m<logClusters[c].SubCluster[n].end; m++) {

						if (logClusters[c].SubCluster[n].strand == 0) {
							baseDots << clusters[logClusters[c].coarse].matches[m].first.pos << "\t" << clusters[logClusters[c].coarse].matches[m].second.pos << "\t" 
									<< clusters[logClusters[c].coarse].matches[m].first.pos + opts.globalK << "\t" << clusters[logClusters[c].coarse].matches[m].second.pos + opts.globalK  << "\t" 
									<<  c << "\t" << logClusters[c].SubCluster[n].strand << endl;
						}
						else {
							baseDots << read.length - clusters[logClusters[c].coarse].matches[m].first.pos - 1 << "\t" << clusters[logClusters[c].coarse].matches[m].second.pos << "\t" 
									<< read.length - clusters[logClusters[c].coarse].matches[m].first.pos - 1 - opts.globalK << "\t"
									<< clusters[logClusters[c].coarse].matches[m].second.pos + opts.globalK << "\t"
									<<  c << "\t" << logClusters[c].SubCluster[n].strand << endl;
						}
					}
				}
				baseDots.close();
			}

			//
			// Get the boundaries of the cluster in genome sequence.
			//
			int nMatch = clusters[logClusters[c].coarse].matches.size();
			GenomePos tPos=clusters[logClusters[c].coarse].tStart;
			int firstChromIndex = genome.header.Find(tPos);
			int lastChromIndex;
			if (nMatch > 1 ) {
				tPos = clusters[logClusters[c].coarse].tEnd;
				lastChromIndex = genome.header.Find(tPos);
			} else { 
				lastChromIndex = firstChromIndex; 
			}
			clusters[logClusters[c].coarse].chromIndex = firstChromIndex;  
			if (firstChromIndex != lastChromIndex ) {
				clusters[logClusters[c].coarse].matches.clear();
				continue;
			}
	
			//
			// Make the anchors reference this chromosome for easier bookkeeping 
			//
			// TODO(Jingwen): remember to add chromOffset back in refinedClusters
			// Note: the qStart, qEnd, tStart, tEnd in refinedClusters[c] are referring to forward strand direction,
			// while the qStart, qEnd, tStart, tEnd in refinedLogClusters[c].SubCluster are referring to forward strand direction;
			// and the coordinates of reversed anchors in refinedClusters[c].matches are in forward strand direction
			GenomePos chromOffset = genome.header.pos[firstChromIndex];
			for (int m=0; m < clusters[logClusters[c].coarse].matches.size(); m++) {
				clusters[logClusters[c].coarse].matches[m].second.pos-=chromOffset;
			}
			GenomePos GenomeClusterEnd = clusters[logClusters[c].coarse].tEnd;
			GenomePos chromEndOffset = genome.header.GetNextOffset(GenomeClusterEnd);
			refinedLogClusters[c].setHp(refinedClusters[c]);


			for (int sc = 0; sc < logClusters[c].SubCluster.size(); sc++) {

				//
				// Get the boundaries of cluster fragment in both sequences.
				//
				// sorted by second.pos and then first.pos
				//
				CartesianTargetSort<GenomeTuple>(clusters[logClusters[c].coarse].matches, logClusters[c].SubCluster[sc].start, logClusters[c].SubCluster[sc].end); 

				//
				// Get shorthand access to alignment boundaries.
				//
				GenomePos genomeClusterSegStart, genomeClusterSegEnd;
				genomeClusterSegStart = logClusters[c].SubCluster[sc].tStart;
				genomeClusterSegEnd = logClusters[c].SubCluster[sc].tEnd;

				int cl = logClusters[c].SubCluster[sc].end - logClusters[c].SubCluster[sc].start;

				int ls, le;
				// Search region starts in window, or beginning of chromosome
				GenomePos wts, wte;
				if ( chromOffset + smallOpts.window > genomeClusterSegStart ) {
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
				readIndex = localIndexes[logClusters[c].SubCluster[sc].strand];

				int lmIndex=0;
				GenomePos alnStart=-1;
				GenomePos alnEnd=0;
				int rfCsize = refinedClusters[c].size();



				for (int lsi=ls; lsi <= le; lsi++) {
					//
					// Find the coordinates in the cluster fragment that start in this local index.
					//
					GenomePos genomeLocalIndexStart = glIndex.seqOffsets[lsi]  - chromOffset;
					GenomePos genomeLocalIndexEnd   = glIndex.seqOffsets[lsi+1] - 1 - chromOffset;

					int matchStart = CartesianTargetLowerBound<GenomeTuple>(clusters[logClusters[c].coarse].matches.begin() + logClusters[c].SubCluster[sc].start,
																clusters[logClusters[c].coarse].matches.begin() + logClusters[c].SubCluster[sc].end, genomeLocalIndexStart);

					int matchEnd = CartesianTargetUpperBound<GenomeTuple>(clusters[logClusters[c].coarse].matches.begin() + logClusters[c].SubCluster[sc].start, 
																clusters[logClusters[c].coarse].matches.begin() + logClusters[c].SubCluster[sc].end, genomeLocalIndexEnd);

					matchStart += logClusters[c].SubCluster[sc].start;
					matchEnd += logClusters[c].SubCluster[sc].start;

					//
					// If there is no overlap with this cluster
					if (matchStart >= logClusters[c].SubCluster[sc].end) {
						continue;
					}
					GenomePos readStart = clusters[logClusters[c].coarse].matches[matchStart].first.pos;
					if (lsi == ls) {
						if (readStart < smallOpts.window) {
							readStart = 0;
						}
						else {
							readStart -= smallOpts.window;
						}
					}
					GenomePos readEnd;
					if (matchEnd > matchStart) {
						readEnd = clusters[logClusters[c].coarse].matches[matchEnd - 1].first.pos;
					}
					else {
						readEnd = clusters[logClusters[c].coarse].matches[matchStart].first.pos + smallOpts.globalK;
					}
					//
					// Expand boundaries of read to match.
					if (lsi == le) {
						if (readEnd + smallOpts.window > read.length) {
							readEnd = read.length; 
						}
						else { 
							readEnd += smallOpts.window;	
						}
					}			
						
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

						CompareLists<LocalTuple>(readIndex->minimizers.begin() + qStartBoundary, readIndex->minimizers.begin() + qEndBoundary, 
														glIndex.minimizers.begin()+ glIndex.tupleBoundaries[lsi], glIndex.minimizers.begin()+ glIndex.tupleBoundaries[lsi+1], 
																smallMatches, smallOpts);

						lmIndex+=smallMatches.size();


						//
						// Do local processing of matches to ensure the region that is searched returns reasonable anchors.
						//
						//DiagonalSort<LocalTuple>(smallMatches); // sort by forward diagonal
						//int smallminDiagCluster = 0;
						//CleanOffDiagonal(smallMatches, smallOpts, smallminDiagCluster);					
						AppendValues<LocalPairs>(refinedClusters[c].matches, smallMatches.begin(), smallMatches.end(), readSegmentStart, genomeLocalIndexStart);
					}
				}
				if (refinedClusters[c].size() == 0) break;
				DiagonalSort<GenomeTuple>(refinedClusters[c].matches.begin() + rfCsize, refinedClusters[c].matches.begin() + refinedClusters[c].matches.size());
				int median = floor((refinedClusters[c].matches.size() - rfCsize)/2) + rfCsize;
				int diagOrigin = diagOrigin = refinedClusters[c].matches[median].first.pos - refinedClusters[c].matches[median].second.pos;
				CleanOffDiagonal(refinedClusters[c].matches, rfCsize, refinedClusters[c].matches.size(), smallOpts, 0, diagOrigin, 600);

				refinedLogClusters[c].SubCluster.push_back(Cluster(rfCsize, refinedClusters[c].matches.size(), logClusters[c].SubCluster[sc].strand));
				if (logClusters[c].SubCluster[sc].strand == 1) SwapStrand(read, smallOpts, refinedClusters[c].matches, rfCsize, refinedClusters[c].matches.size());	
				refinedLogClusters[c].SetSubClusterBoundariesFromMatches(smallOpts, sc);
				if (logClusters[c].ISsecondary == 1) {
					refinedLogClusters[c].ISsecondary = 1;
					refinedLogClusters[c].primary = logClusters[c].primary;
				}
				else {
					refinedLogClusters[c].secondary =  logClusters[c].secondary;
				}
			}

			SetClusterBoundariesFromSubCluster(refinedClusters[c], opts, refinedLogClusters[c]);
			refinedClusters[c].chromIndex = clusters[logClusters[c].coarse].chromIndex;
			refinedClusters[c].coarse = c;
			refinedClusters[c].strands.resize(refinedClusters[c].matches.size()); 
			refinedClusters[c].coarseSubCluster.resize(refinedClusters[c].matches.size());
			SetCoarseFromSubClusters(refinedClusters[c], refinedLogClusters[c]);
			// refinedClusters[c].strands keeps track of the strand direction of every anchors in refinedClusters[c].matches

			if (opts.dotPlot) {
				stringstream outNameStrm;
				outNameStrm << baseName + "." << c << ".orig.dots";
				ofstream baseDots(outNameStrm.str().c_str());
				for (int m=0; m<refinedLogClusters[c].SubCluster.size(); ++m) {

					for (int n=refinedLogClusters[c].SubCluster[m].start; n<refinedLogClusters[c].SubCluster[m].end; n++) {
						if (refinedLogClusters[c].SubCluster[m].strand == 0) {
							baseDots << refinedClusters[c].matches[n].first.pos << "\t" << refinedClusters[c].matches[n].second.pos << "\t" 
										<< refinedClusters[c].matches[n].first.pos + smallOpts.globalK << "\t" 
										<< refinedClusters[c].matches[n].second.pos + smallOpts.globalK << "\t"  
										<< m << "\t"
										<< refinedLogClusters[c].SubCluster[m].strand << endl;
						}
						else {
							baseDots << refinedClusters[c].matches[n].first.pos << "\t" << refinedClusters[c].matches[n].second.pos + smallOpts.globalK << "\t"
										<< refinedClusters[c].matches[n].first.pos + smallOpts.globalK << "\t" 
										<< refinedClusters[c].matches[n].second.pos  << "\t"  
										<< m << "\t"
										<< refinedLogClusters[c].SubCluster[m].strand << endl;							
						}
					}
				}
				baseDots.close();
			}
		}


		//
		// Remove clusters under some minimal number of anchors. By default this is one. 
		//
		//RemoveEmptyClusters(refinedClusters);
		
		vector<SegAlignmentGroup> alignments(refinedClusters.size()); // alignments[i] stores the alignment for one chain
		for (int r = 0; r < refinedClusters.size(); ++r) {   

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
			// TODO(Jingwen): already did the above DiagonalSort and CleanoffDiagonal above
			/*
			DiagonalSort<GenomeTuple>(refinedClusters[r].matches); // sort first by forward diagonal and then by first.pos
			CleanOffDiagonal(refinedClusters[r].matches, smallOpts);
			*/

			if (opts.dotPlot) {
				stringstream outNameStrm;
				outNameStrm << baseName + "." << r << ".clean.dots";
				ofstream baseDots(outNameStrm.str().c_str());
				for (int m=0; m < refinedClusters[r].matches.size(); m++) {
					baseDots  << refinedClusters[r].matches[m].first.pos << "\t" 
							  << refinedClusters[r].matches[m].second.pos << "\t" 
							  << smallOpts.globalK  << "\t" 
							  << r << endl;
				}
				baseDots.close();
			}


			if (refinedClusters[r].matches.size() == 0) {
				continue;
			}
			if (refinedClusters[r].matches.size() > read.length or opts.refineLevel & REF_LOC == 0) {
				refinedClusters[r].matches = clusters[refinedClusters[r].coarse].matches;
				continue;
			}

			bool ReverseOnly = 1; // ReverseOnly == 1 means there are only reverse matches
								  // If ReverseOnly == 1, then only inserting two points s2, e2 for every reversed match in 1-st SDP.
			for (int m = 0; m < refinedLogClusters[r].SubCluster.size(); ++m) {
				if (refinedLogClusters[r].SubCluster[m].strand != 1) ReverseOnly = 0;
			} 

			//
			// At this point in the code, refinedClusters[r].matches is a list
			// of k-mer matches between the read and the genome. 
			
			//
			// This is where the code should go for merging colinear matches

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


			vector<ClusterCoordinates> mergedAnchors;
			if (opts.mergeClusters) {

				mergedAnchors.clear();
				bool WholeReverseDirection = 1; // WholeReverseDirection == 1 means the read is reversely mapped to reference;
												// This is important when deciding the boundary in MergeAnchor function
				GenomePos prev_qEnd = 0, cur_qEnd = 0;
				for (int m = 0; m < refinedLogClusters[r].SubCluster.size(); ++m) {

					if (m == 0) prev_qEnd = refinedLogClusters[r].SubCluster[m].qEnd;
					else {
						cur_qEnd = refinedLogClusters[r].SubCluster[m].qEnd;
						if (prev_qEnd >= cur_qEnd) WholeReverseDirection = 0;
					}
				}

				//MergeClusters(smallOpts, refinedClusters, vt, r);
				//mergeClusters (smallOpts, refinedClusters[r].matches, vt, r, baseName);
				MergeAnchors (smallOpts, refinedClusters, refinedLogClusters, r, mergedAnchors, WholeReverseDirection);
				int MergeAnchorsNum = 0;
				for (int m = 0; m < mergedAnchors.size(); m++) {
					MergeAnchorsNum += mergedAnchors[m].end - mergedAnchors[m].start;
				}


				//
				// Output the anchor efficiency file
				/*
				if (refinedClusters[r].matches.size() != 0) {
					
					stringstream outNameStrm;
					outNameStrm << "AnchorEfficiency.tab";
					ofstream baseDots;
					baseDots.open(outNameStrm.str().c_str(), std::ios::app);
					baseDots << mergedAnchors.size() << "\t" << refinedClusters[r].matches.size() << endl;
					baseDots.close();						
				}
				*/

				// TODO(Jingwen): only for debug the new MergeAnchors function from MergeAnchors.h
				if (opts.dotPlot) {
					stringstream outNameStrm;
					outNameStrm << baseName + "." << r << ".merged.dots";
					ofstream baseDots(outNameStrm.str().c_str());
					for (int m=0; m < mergedAnchors.size(); m++) {

						if (mergedAnchors[m].strand == 0) {
							// chain stores indices which refer to elments in vt
							baseDots << mergedAnchors[m].qStart << "\t" 
									 << mergedAnchors[m].tStart << "\t" 
									 << mergedAnchors[m].qEnd << "\t" 
									 << mergedAnchors[m].tEnd << "\t"
									 << r << "\t"
									 << mergedAnchors[m].strand << endl;									
						}		
						else {
							baseDots << mergedAnchors[m].qStart << "\t" 
									 << mergedAnchors[m].tEnd << "\t" 
									 << mergedAnchors[m].qEnd << "\t" 
									 << mergedAnchors[m].tStart << "\t"
									 << r << "\t"
									 << mergedAnchors[m].strand << endl;								
						}	

					}
					baseDots.close();
				}


				/*
				//TODO(Jingwen): Only for debug
				if (opts.dotPlot) {
					stringstream outNameStrm;
					outNameStrm << baseName + "." << r << ".merged.dots";
					ofstream baseDots(outNameStrm.str().c_str());
					for (int m=0; m < vt.size(); m++) {
						baseDots << vt[m].qStart << "\t" 
								 << vt[m].tStart << "\t" 
								 << vt[m].qEnd << "\t" 
								 << vt[m].tEnd << "\t" 
								 << m << endl;
					}
					baseDots.close();
				}
				*/
				/*
				if (opts.dotPlot) {
					stringstream outNameStrm;
					outNameStrm << baseName + "." << r << ".merged.real.dots";
					ofstream baseDots(outNameStrm.str().c_str());
					for (int m=0; m < vt.size(); m++) {
						for (int n = vt[m].start; n < vt[m].end; ++n) {
							baseDots << refinedClusters[r].matches[n].first.pos << "\t" << refinedClusters[r].matches[n].second.pos << "\t" << smallOpts.globalK << "\t" << m << endl;

						}
					}
					baseDots.close();
				}
				*/


			}


			// Perform sparse chaining, uses time O(n (log n)^2).
			//
			// Merge anchors that are not overlapping

			//
			// The input are a buntch of fragments which are either stored in vector<Cluster> vt or refinedClusters[r].matches
			// 
			vector<unsigned int> chain; // chain contains the index of the fragments involved in the result of SparseDP.
			if (opts.SparseDP) {
				chain.clear();
				if (opts.mergeClusters and mergedAnchors.size() < 30000 and mergedAnchors.size() > 0) {
					if (refinedClusters[r].matches.size()/((float)(min(refinedClusters[r].qEnd - refinedClusters[r].qStart, refinedClusters[r].tEnd - refinedClusters[r].tStart))) < 0.1) {
						SparseDP(mergedAnchors, chain, smallOpts, LookUpTable, ReverseOnly, 5); 
					}
					else {
						SparseDP(mergedAnchors, chain, smallOpts, LookUpTable, ReverseOnly);
					}
				}
				else if (refinedClusters[r].matches.size() < 30000) {
					if (refinedClusters[r].matches.size()/((float)(min(refinedClusters[r].qEnd - refinedClusters[r].qStart, refinedClusters[r].tEnd - refinedClusters[r].tStart))) < 0.1) {
						SparseDP(refinedClusters[r].matches, chain, smallOpts, LookUpTable, refinedLogClusters[r], refinedClusters[r].strands, ReverseOnly, 10);
					}
					else { 
						// If anchors are unmerged, then we need to give a higher anchor value to every anchor
						// Since gap cost of chaining is higher.
						SparseDP(refinedClusters[r].matches, chain, smallOpts, LookUpTable, refinedLogClusters[r], refinedClusters[r].strands, ReverseOnly, 5);
					}
				}
			}
			else { 
				chain.clear();
				// TODO(Jingwen): add linear sdp here?
			}


			// TODO(Jingwen): Only for debug
			if (opts.dotPlot) {
				stringstream outNameStrm;
				outNameStrm << baseName + "." << r << ".first-sdp.dots";
				ofstream baseDots(outNameStrm.str().c_str());
				for (int c = 0; c < chain.size(); c++) {
					if (opts.mergeClusters) {

						if (mergedAnchors[chain[c]].strand == 0) {
							// chain stores indices which refer to elments in vt
							baseDots << mergedAnchors[chain[c]].qStart << "\t" 
									 << mergedAnchors[chain[c]].tStart << "\t" 
									 << mergedAnchors[chain[c]].qEnd << "\t" 
									 << mergedAnchors[chain[c]].tEnd << "\t"
									 << r << "\t"
									 << mergedAnchors[chain[c]].strand << endl;									
						}		
						else {
							baseDots << mergedAnchors[chain[c]].qStart << "\t" 
									 << mergedAnchors[chain[c]].tEnd << "\t" 
									 << mergedAnchors[chain[c]].qEnd << "\t" 
									 << mergedAnchors[chain[c]].tStart << "\t"
									 << r << "\t"
									 << mergedAnchors[chain[c]].strand << endl;								
						}	
						//TODO(Jingwen): delete the following code later
						if (c != chain.size() - 1) {
							assert(mergedAnchors[chain[c]].qStart < mergedAnchors[chain[c+1]].qStart);
							//assert(refinedClusters[r].matches[chain[c]].second.pos >= refinedClusters[r].matches[chain[c+1]].second.pos + smallOpts.globalK);
						}			
					}
					else {
						// chain stores indices which refer to elements in refinedClusters[r].matches
						if (refinedClusters[r].strands[chain[c]] == 0) {
							baseDots << refinedClusters[r].matches[chain[c]].first.pos << "\t" 
									 << refinedClusters[r].matches[chain[c]].second.pos << "\t" 
									 << refinedClusters[r].matches[chain[c]].first.pos + smallOpts.globalK << "\t"
									 << refinedClusters[r].matches[chain[c]].second.pos + smallOpts.globalK << "\t"
									 << r << "\t"
									 << refinedClusters[r].strands[chain[c]] << endl;								
						}
						else {
							baseDots << refinedClusters[r].matches[chain[c]].first.pos << "\t" 
									 << refinedClusters[r].matches[chain[c]].second.pos + smallOpts.globalK << "\t" 
									 << refinedClusters[r].matches[chain[c]].first.pos + smallOpts.globalK << "\t"
									 << refinedClusters[r].matches[chain[c]].second.pos << "\t"
									 << r << "\t"
									 << refinedClusters[r].strands[chain[c]] << endl;	
						}		
						//TODO(Jingwen): delete the following code later
						if (c != chain.size() - 1) {
							assert(refinedClusters[r].matches[chain[c]].first.pos < refinedClusters[r].matches[chain[c+1]].first.pos);
							//assert(refinedClusters[r].matches[chain[c]].second.pos >= refinedClusters[r].matches[chain[c+1]].second.pos + smallOpts.globalK);
						}				
					}	
				}
				baseDots.close();
			}


			GenomePairs tupChain;
			vector<Cluster> tupChainClusters;
			tupChain.clear();
			tupChainClusters.clear();

			if (opts.mergeClusters and mergedAnchors.size() > 0) {
				//
				RemovePairedIndels(chain, mergedAnchors, opts);	

				//
				RemoveSpuriousAnchors(chain, opts, mergedAnchors, refinedLogClusters[r]);

				//
				// Add small anchors to tupChain. (Use greedy algorithm to make small anchors not overlap with each other)
				// 
				vector<int> tupChainStrand;
				tupChainStrand.clear();
				tupChain.clear();
				for (int ch=0; ch < chain.size(); ch++) { 

					//cerr << "ch: "<< ch<< endl<<endl;
					//cerr << "chain[" << ch <<"]   " << chain[ch] << endl; 
					//
					// Use greedy algorithm to make small anchors not overlap with each other
					int id = mergedAnchors[chain[ch]].start;
					int c = 0;
					//cerr << "ch: " << ch << endl;
					//cerr << "id: " << id << endl;
					//cerr << "mergedAnchors[chain[ch]].strand: " << mergedAnchors[chain[ch]].strand << endl;
					if (mergedAnchors[chain[ch]].strand == 1) { // rev strand
						//cerr << "id: " << id << endl;
						tupChain.push_back(GenomePair(GenomeTuple(0, read.length - (refinedClusters[r].matches[id].first.pos + smallOpts.globalK)), 
												GenomeTuple(0, refinedClusters[r].matches[id].second.pos)));	
						tupChainStrand.push_back(mergedAnchors[chain[ch]].strand);
						c++;

						GenomePos qEnd = refinedClusters[r].matches[id].first.pos + smallOpts.globalK;
						GenomePos tEnd = refinedClusters[r].matches[id].second.pos; // should not use -1 here, because refinedClusters[r].matches[id].second.pos might be 0
						//cerr << "qEnd: " << qEnd << endl;
						//cerr << "tEnd: " << tEnd << endl;
						int ce = id + 1;
						while (ce < mergedAnchors[chain[ch]].end) {

							while (ce < mergedAnchors[chain[ch]].end and 
									(refinedClusters[r].matches[ce].first.pos < qEnd or refinedClusters[r].matches[ce].second.pos + smallOpts.globalK >= tEnd)) { 
								++ce;
							}
							if (ce < mergedAnchors[chain[ch]].end) {

								// The following if condition aims to remove paired indel inside each merged fragment
							/*	if (c >= 3) {
									int si = tupChain.size() - 1;
									int prevGap = (int)(tupChain[si].second.pos - tupChain[si].first.pos) - 
											  (int)(tupChain[si-1].second.pos - tupChain[si-1].first.pos);
									int nextGap = (int)(refinedClusters[r].matches[ce].second.pos - read.length 
														+ (refinedClusters[r].matches[ce].first.pos + smallOpts.globalK)) - 
									              (int)(tupChain[si].second.pos - tupChain[si].first.pos);
									if (sign(prevGap)!= sign(nextGap) and abs(prevGap) + abs(nextGap) > abs(prevGap + nextGap)) {
										tupChain.pop_back();
										tupChainStrand.pop_back();
										c--;
									}
								}
							*/
								tupChain.push_back(GenomePair(GenomeTuple(0, read.length - (refinedClusters[r].matches[ce].first.pos + smallOpts.globalK)), 
														GenomeTuple(0, refinedClusters[r].matches[ce].second.pos)));
								tupChainStrand.push_back(mergedAnchors[chain[ch]].strand);	
								c++;

								//
								// TODO(Jingwen): only for deubg and delete later
								int s = tupChain.size() - 1;
								assert(tupChain[s].first.pos + smallOpts.globalK <= tupChain[s-1].first.pos);
								assert(tupChain[s].second.pos + smallOpts.globalK <= tupChain[s-1].second.pos);

								id = ce;
								qEnd = refinedClusters[r].matches[id].first.pos + smallOpts.globalK;
								tEnd = refinedClusters[r].matches[id].second.pos;
								++ce;
							}
						}

					}
					else { // forward strand
						tupChain.push_back(GenomePair(GenomeTuple(0, refinedClusters[r].matches[id].first.pos), 
													GenomeTuple(0, refinedClusters[r].matches[id].second.pos)));	
						tupChainStrand.push_back(mergedAnchors[chain[ch]].strand);
						c++;
						//cerr << "push back: (" << refinedClusters[r].matches[id].first.pos << ", " <<  refinedClusters[r].matches[id].second.pos << ")" << endl;
						//cerr << "id: " << id << endl;

						GenomePos qEnd = refinedClusters[r].matches[id].first.pos + smallOpts.globalK;
						GenomePos tEnd = refinedClusters[r].matches[id].second.pos + smallOpts.globalK;	

						int ce = id + 1;
						while (ce < mergedAnchors[chain[ch]].end) {
							
							while (ce < mergedAnchors[chain[ch]].end and 
									(refinedClusters[r].matches[ce].first.pos < qEnd or refinedClusters[r].matches[ce].second.pos < tEnd)) { 
								++ce;
							}
							if (ce < mergedAnchors[chain[ch]].end) {

								// The following if condition aims to remove paired indel inside each merged fragment
							/*	if (c >= 3) {
									int si = tupChain.size() - 1;
									int prevGap = (int)(tupChain[si].second.pos - tupChain[si].first.pos) - 
											  (int)(tupChain[si-1].second.pos - tupChain[si-1].first.pos);
									int nextGap = (int)(refinedClusters[r].matches[ce].second.pos - 
														refinedClusters[r].matches[ce].first.pos) - 
									              (int)(tupChain[si].second.pos - tupChain[si].first.pos);
									if (sign(prevGap)!= sign(nextGap) and abs(prevGap) + abs(nextGap) > abs(prevGap + nextGap)) {
										tupChain.pop_back();
										tupChainStrand.pop_back();
										c--;
									}
								}
							*/
								tupChain.push_back(GenomePair(GenomeTuple(0, refinedClusters[r].matches[ce].first.pos), 
														GenomeTuple(0, refinedClusters[r].matches[ce].second.pos)));
								tupChainStrand.push_back(mergedAnchors[chain[ch]].strand);
								c++;
								//cerr << "push back: (" << refinedClusters[r].matches[ce].first.pos << ", " <<  refinedClusters[r].matches[ce].second.pos << ")" << endl;
								//cerr << "ce: " << ce << endl;

								//
								// TODO(Jingwen): Only for debug and delete this later
								int s = tupChain.size() - 1;
								assert(tupChain[s].first.pos >= tupChain[s-1].first.pos + smallOpts.globalK);
								assert(tupChain[s].second.pos >= tupChain[s-1].second.pos + smallOpts.globalK);								
	
								id = ce;
								qEnd = refinedClusters[r].matches[id].first.pos + smallOpts.globalK;
								tEnd = refinedClusters[r].matches[id].second.pos + smallOpts.globalK;	
								++ce;
							}
						}
					}
				}

				assert(tupChain.size() == tupChainStrand.size());// TODO(Jingwen): only for debug and delete this later
				// Remove paired indels in merged fragments
				RemovePairedIndels(tupChain, tupChainStrand, smallOpts);

				int cs = 0, ce = 0;
				while (ce < tupChainStrand.size()) {
					if (tupChainStrand[cs] == tupChainStrand[ce]) ce++;
					else {
						tupChainClusters.push_back(Cluster(cs, ce, tupChainStrand[cs]));
						//if (ce - cs >= smallOpts.mintupChainClustersize) tupChainClusters.push_back(Cluster(cs, ce, tupChainStrand[cs]));
						cs = ce;
					}
				}
				if (ce == tupChainStrand.size() and tupChain.size()) {
						tupChainClusters.push_back(Cluster(cs, ce, tupChainStrand[cs]));
				}
			}
			else {
				//
				// Remove paired indels
				RemovePairedIndels(chain, refinedClusters[r].matches, smallOpts);

				//
				RemoveSpuriousAnchors(chain, smallOpts, refinedClusters[r], refinedLogClusters[r]);
				//
				// chain stores indices which refer to elements in refinedClusters[r].matches
				for (int ch=0; ch < chain.size(); ch++) {

					assert(refinedClusters[r].matches[chain[ch]].first.pos <= read.length); // TODO(Jingwen): delete this after debug
					// Swap anchors of reverse strand back to reverse strand coordinates, making it easier to fill up the gap
					if (refinedClusters[r].strands[chain[ch]] == 1) {
						tupChain.push_back(GenomePair(GenomeTuple(0, read.length - (refinedClusters[r].matches[chain[ch]].first.pos + smallOpts.globalK)), 
													GenomeTuple(0, refinedClusters[r].matches[chain[ch]].second.pos)));
					}
					else {
						tupChain.push_back(GenomePair(GenomeTuple(0, refinedClusters[r].matches[chain[ch]].first.pos), 
															GenomeTuple(0, refinedClusters[r].matches[chain[ch]].second.pos)));						
					}
				}

				int cs = 0, ce = 0;
				while (ce < chain.size()) {
					if (refinedClusters[r].strands[chain[cs]] == refinedClusters[r].strands[chain[ce]]) ce++;
					else {
						tupChainClusters.push_back(Cluster(cs, ce, refinedClusters[r].strands[chain[cs]]));
						//if (ce - cs >= smallOpts.mintupChainClustersize) tupChainClusters.push_back(Cluster(cs, ce, refinedClusters[r].strands[chain[cs]]));
						cs = ce;
					}
				}
				if (ce == chain.size() and cs < chain.size()) {
						tupChainClusters.push_back(Cluster(cs, ce, refinedClusters[r].strands[chain[cs]]));
				}
			}


			//(TODO)Jingwen: For Debug(remove this later)
			if (opts.dotPlot) {
				stringstream outNameStrm;
				outNameStrm << baseName + "." << r << ".first-sdp-clean.dots";
				ofstream baseDots(outNameStrm.str().c_str());
				for (int m = 0; m < tupChainClusters.size(); m++) {
					for (int c = tupChainClusters[m].start; c < tupChainClusters[m].end; c++) {
	
						if (tupChainClusters[m].strand == 0) {
							baseDots << tupChain[c].first.pos << "\t" 
									 << tupChain[c].second.pos << "\t" 
									 << tupChain[c].first.pos + smallOpts.globalK << "\t" 
									 << tupChain[c].second.pos + smallOpts.globalK << "\t"
									 << m << "\t"
									 << tupChainClusters[m].strand << endl;								
						}
						else {
							baseDots << read.length - tupChain[c].first.pos - smallOpts.globalK << "\t" 
									 << tupChain[c].second.pos + smallOpts.globalK<< "\t" 
									 << read.length - tupChain[c].first.pos << "\t" 
									 << tupChain[c].second.pos << "\t"
									 << m << "\t"
									 << tupChainClusters[m].strand << endl;	

						}	
						//TODO(Jingwe): the following code is only for debug
						if (tupChainClusters[m].strand == 1 and c != tupChainClusters[m].end - 1) {
							assert(tupChain[c].first.pos > tupChain[c+1].first.pos);
							assert(tupChain[c].second.pos >= tupChain[c+1].second.pos + smallOpts.globalK);
						}
					}
				}
			}



			if (tupChain.size() == 0) {
				refinedClusters[r].matches.clear();
				continue;
			}

			// TODO(Jingwen): the following code is for removepairedIndels, check later how to modify it
			vector<Cluster> chainClust;
			Options diagOpts;
			diagOpts = smallOpts;
			diagOpts.maxDiag=15;
			diagOpts.minClusterSize=1;

			// remove fragments which are in the middle of an insertion and a deletion OR a deletion and an insertion. 
			//StoreDiagonalClusters(tupChain, chainClust, diagOpts, true); ///Jingwen adds this here, otherwise chainClust is empty
			//RemovePairedIndels(tupChain, chainClust, smallOpts);  
			// TODO(Jingwen): make this work
			//RemovePairedIndels(curReadEnd, curGenomeEnd, nextReadStart, nextGenomeStart, tupChain, refinedClusters[r].matches, smallOpts);

			//(TODO)Jingwen: For Debug(remove this later)
			if (opts.dotPlot) {
				stringstream outNameStrm;
				outNameStrm << baseName + "." << r << ".first-sdp-clean-removeIndel.dots";
				ofstream baseDots(outNameStrm.str().c_str());
				for (int c = 0; c < tupChain.size(); c++) {
					// chain stores indices which refer to elments in vt
					baseDots << tupChain[c].first.pos << "\t" 
							 << tupChain[c].second.pos << "\t" 
							 << smallOpts.globalK << "\t"
							 << r << endl;							
				}
				baseDots.close();
			}	
			/*
			int prevq=0;
			int prevt=0;
			
			//		cout << "Chain is on " << chainClust.size() << " diagonals " << endl;

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


			// tupChain stores the fragments from first sdp. the length of fragments in tupChain is smallOpts.globalK
			//GenomePos chainGenomeStart = tupChain[0].second.pos;
			//GenomePos chainGenomeEnd   = tupChain[tupChain.size()-1].second.pos + smallOpts.globalK;
			//GenomePos chainReadStart = tupChain[0].first.pos;
			//GenomePos chainReadEnd   = tupChain[tupChain.size()-1].first.pos + smallOpts.globalK;

			int chromIndex = refinedClusters[r].chromIndex;
				
			//
			// Create subsequences that will be used to generate the alignment.  Gaps should be inserted 
			// with respect to an offset from chainGenomeStart and chainReadStart
			//

			// TODO(Jingwen): gapOpts should be replaced by tinyOpts
			Options gapOpts=opts;
			gapOpts.globalMaxFreq=5;
			gapOpts.globalK=7;


			for (int s = 0; s < tupChainClusters.size(); s++) {
				if (tupChainClusters[s].end == tupChainClusters[s].start) {
					continue;
				}
				if (tupChainClusters[s].strand == 1) reverse(tupChain.begin() + tupChainClusters[s].start, tupChain.begin() + tupChainClusters[s].end);
				vector<GenomeTuple> gapReadTup, gapGenomeTup;
				GenomePairs gapPairs;
				vector<GenomePairs> refinedChains(tupChainClusters[s].end - tupChainClusters[s].start - 1); // Note: refinedChains[i] stores the gap-fragments which locate after chain[i]
				vector<int> refinedChainsLength(tupChainClusters[s].end - tupChainClusters[s].start - 1, -1); // refinedChainsLength[i] stores the tinyOpts.globalK of the gap-fragments which locate after chain[i]
				
				vector<int> scoreMat;
				vector<Arrow> pathMat;

				int chainLength = tupChainClusters[s].end - tupChainClusters[s].start;
				for (int c = tupChainClusters[s].start; chainLength > 0 and c < tupChainClusters[s].end - 1; c++) {
					//cerr << "tupChainClusters[s].start: " << tupChainClusters[s].start << "  tupChainClusters[s].end: " << tupChainClusters[s].end << endl;
					//cerr << "c: " << c << endl;

					GenomePos curGenomeEnd = tupChain[c].second.pos + smallOpts.globalK;
					GenomePos curReadEnd = tupChain[c].first.pos + smallOpts.globalK;

					GenomePos nextGenomeStart = tupChain[c+1].second.pos;
					GenomePos nextReadStart = tupChain[c+1].first.pos;

					assert(nextReadStart >= curReadEnd);
					GenomePos subreadLength = nextReadStart - curReadEnd; 
					assert(nextGenomeStart >= curGenomeEnd);
					GenomePos subgenomeLength = nextGenomeStart - curGenomeEnd;


					if (nextReadStart > curReadEnd and nextGenomeStart > curGenomeEnd) {

						if (subreadLength > 50 and subgenomeLength > 50 and opts.refineLevel & REF_DYN ) {

							// TODO(Jingwen): should we only consider minLen? because min(subreadLength, subgenomeLength) is the length of possible matches
							GenomePos maxLen = max(subreadLength, subgenomeLength);
							//GenomePos maxLen = min(subreadLength, subgenomeLength);
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
							StoreMinimizers<GenomeTuple, Tuple>(genome.seqs[chromIndex] + curGenomeEnd, subgenomeLength, tinyOpts.globalK, 1, gapGenomeTup, false);

							sort(gapGenomeTup.begin(), gapGenomeTup.end());
							StoreMinimizers<GenomeTuple, Tuple>(strands[tupChainClusters[s].strand] + curReadEnd, subreadLength, tinyOpts.globalK, 1, gapReadTup, false);
							sort(gapReadTup.begin(), gapReadTup.end());
							CompareLists(gapReadTup.begin(), gapReadTup.end(), gapGenomeTup.begin(), gapGenomeTup.end(), gapPairs, tinyOpts);
							//CompareLists(gapReadTup.begin(), gapReadTup.end(), gapGenomeTup.begin(), gapGenomeTup.end(), gapPairs, gapOpts);
							//cerr << "gapGenomeTup.size(): " << gapGenomeTup.size() << endl;
							//cerr << "gapReadTup.size(): " << gapReadTup.size() << endl;
							//
							// Remove egregious off diagonal seeds
							//
							DiagonalSort<GenomeTuple>(gapPairs); // sort gapPairs by forward diagonals and then by first.pos
							
							int tinyDiagStart = curReadEnd - curGenomeEnd;
							int tinyDiagEnd =  nextReadStart - nextGenomeStart;
							int diagDiff = abs((int) tinyDiagStart - (int) tinyDiagEnd);

							// TODO(Jingwen): after applying this CleanOffDiagonal, no gapPairs left
							// CleanOffDiagonal(gapPairs, tinyOpts, 0, diagDiff);
							//CartesianTargetSort<GenomeTuple>(gapPairs);

							for(int rm=0; rm < gapPairs.size(); rm++) {
								gapPairs[rm].first.pos  += curReadEnd;
								gapPairs[rm].second.pos += curGenomeEnd;
								assert(gapPairs[rm].first.pos < read.length); // TODO(Jingwen): delete this after debug
							}

							// gapPairs stores all the fragments in the gap. the length of the fragment here is tinyOpts.globalK
							// gapChain stores the index of fragments that are involved in the result of sdp
							vector<unsigned int> gapChain; 
							if (gapPairs.size() > 0) {
								if (opts.SparseDP) {
									gapChain.clear();
									if (gapPairs.size() < 20000) {
										if (gapPairs.size()/((float)min(subreadLength, subgenomeLength)) < 0.1) {
											SparseDP_ForwardOnly(gapPairs, gapChain, tinyOpts, LookUpTable, 10); // TODO(Jingwen): change the rate to customized number
										}
										else { 
											SparseDP_ForwardOnly(gapPairs, gapChain, tinyOpts, LookUpTable); 
										}
									}
								}
							}

							RemovePairedIndels(curReadEnd, curGenomeEnd, nextReadStart, nextGenomeStart, gapChain, gapPairs, tinyOpts);
							RemoveSpuriousAnchors(gapChain, tinyOpts, gapPairs);
							//cerr << "gapChain.size()/((float)min(subreadLength, subgenomeLength)): " << gapChain.size()/((float)min(subreadLength, subgenomeLength)) << endl;
							if (gapChain.size()/((float)max(subreadLength, subgenomeLength)) < 0.02) gapChain.clear();

							if (opts.dotPlot) {
								stringstream outNameStrm;
								outNameStrm << baseName + "." << r << ".second-sdp.dots";
								ofstream baseDots;
								baseDots.open(outNameStrm.str().c_str(), std::ios::app);
								for (int c = 0; c < gapChain.size(); c++) {
									// chain stores indices which refer to elements in refinedClusters[r].matches
									if (tupChainClusters[s].strand == 0) {
										baseDots << gapPairs[gapChain[c]].first.pos << "\t" 
												 << gapPairs[gapChain[c]].second.pos << "\t" 
												 << tinyOpts.globalK + gapPairs[gapChain[c]].first.pos << "\t"
												 << tinyOpts.globalK + gapPairs[gapChain[c]].second.pos << "\t"
												 << r << "\t"
												 << tupChainClusters[s].strand << endl;														
									}
									else {	
										baseDots << read.length - gapPairs[gapChain[c]].first.pos - tinyOpts.globalK << "\t" 
												 << gapPairs[gapChain[c]].second.pos + tinyOpts.globalK << "\t" 
												 << read.length - gapPairs[gapChain[c]].first.pos << "\t"
												 << gapPairs[gapChain[c]].second.pos << "\t"
												 << r << "\t"
												 << tupChainClusters[s].strand << endl;												
									}
								}
								baseDots.close();
							}

							for (unsigned int u = 0; u < gapChain.size(); ++u) {
								refinedChains[c - tupChainClusters[s].start].push_back(gapPairs[gapChain[u]]);
								assert(refinedChains[c - tupChainClusters[s].start].back().first.pos <= read.length);

								// TODO(Jingwen): delete the debug code
								if (refinedChains[c - tupChainClusters[s].start].size() > 1) { 
									int last = refinedChains[c - tupChainClusters[s].start].size();
									assert(refinedChains[c - tupChainClusters[s].start][last -2].first.pos + tinyOpts.globalK <= refinedChains[c - tupChainClusters[s].start][last -1].first.pos);
									assert(refinedChains[c - tupChainClusters[s].start][last -2].second.pos + tinyOpts.globalK <= refinedChains[c - tupChainClusters[s].start][last -1].second.pos);						
								}
							}
							refinedChainsLength[c - tupChainClusters[s].start] = tinyOpts.globalK;	
							//cerr<< "s: " << s << "  c: " << c << "  gapPairs.size(): " << gapPairs.size() << "  refinedChains[c].size(): "<< refinedChains[c].size() << endl;
						}
					}	
				}



				//
				// Refine and store the alignment
				//
				//
				// The alignment is on a substring that starts at the beginning of the first chain.
				//
				Alignment *alignment = new Alignment(strands[tupChainClusters[s].strand], read.seq, read.length, read.name, tupChainClusters[s].strand, genome.seqs[chromIndex],  
													genome.lengths[chromIndex], genome.header.names[chromIndex], chromIndex);

				alignments[r].SegAlignment.push_back(alignment);
				if (refinedLogClusters[r].ISsecondary == 0) alignments[r].secondary = refinedLogClusters[r].secondary;
				else {
					alignments[r].ISsecondary = refinedLogClusters[r].ISsecondary;
					alignments[r].primary = refinedLogClusters[r].primary;					
				}

				for (int c = tupChainClusters[s].start; chainLength > 0 and c < tupChainClusters[s].end - 1; c++) {
					//
					// Chain is with respect to full sequence
					//
					GenomePos curGenomeEnd     = tupChain[c].second.pos + smallOpts.globalK;
					GenomePos curReadEnd       = tupChain[c].first.pos + smallOpts.globalK;
					GenomePos nextGenomeStart  = tupChain[c+1].second.pos;
					GenomePos nextReadStart    = tupChain[c+1].first.pos;
					int curRefinedReadEnd      = curReadEnd;
					int curRefinedGenomeEnd    = curGenomeEnd;
					int nextRefinedReadStart   = nextReadStart;
					int nextRefinedGenomeStart = nextGenomeStart;

					if (opts.dotPlot) {
						if (tupChainClusters[s].strand == 0) {
							dotFile << tupChain[c].first.pos << "\t" 
									<< tupChain[c].second.pos << "\t" 
									<< tupChain[c].first.pos + smallOpts.globalK << "\t"
									<< tupChain[c].second.pos + smallOpts.globalK << "\t"
									<< r << "\t"
									<< tupChainClusters[s].strand << endl;					
						}
						else {
							dotFile << read.length - tupChain[c].first.pos - smallOpts.globalK << "\t" 
									<< tupChain[c].second.pos  + smallOpts.globalK << "\t" 
									<< read.length - tupChain[c].first.pos << "\t"
									<< tupChain[c].second.pos << "\t" 
									<< r << "\t"
									<< tupChainClusters[s].strand << endl;							
						}
					}

					alignment->blocks.push_back(Block(tupChain[c].first.pos, tupChain[c].second.pos, smallOpts.globalK)); 
					//cerr << "alignment->blocks.back().size(): " << alignment->blocks.back().size() << endl;
					//cerr << "alignments[r].SegAlignment.back().size(): " << alignments[r].SegAlignment.back().size() << endl;

					if (alignment->blocks.size() > 1) {
						int last=alignment->blocks.size();
						assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].qPos);
						assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].tPos);
					}

					for (int cs = 0; cs < refinedChains[c - tupChainClusters[s].start].size(); cs++) {
						//
						// Refined anchors are with respect to the chained sequence

						nextRefinedReadStart = refinedChains[c - tupChainClusters[s].start][cs].first.pos;
						nextRefinedGenomeStart = refinedChains[c - tupChainClusters[s].start][cs].second.pos;
						
						if (opts.dotPlot) {
							if (tupChainClusters[s].strand == 0) {
								dotFile << refinedChains[c - tupChainClusters[s].start][cs].first.pos << "\t" 
										<< refinedChains[c - tupChainClusters[s].start][cs].second.pos << "\t" 
										<< refinedChains[c - tupChainClusters[s].start][cs].first.pos + refinedChainsLength[c - tupChainClusters[s].start] << "\t"
										<< refinedChains[c - tupChainClusters[s].start][cs].second.pos + refinedChainsLength[c - tupChainClusters[s].start] << "\t" 
										<< r << "\t"
										<< tupChainClusters[s].strand << endl;								
							}
							else {
								dotFile << read.length - refinedChains[c - tupChainClusters[s].start][cs].first.pos - refinedChainsLength[c - tupChainClusters[s].start] << "\t" 
										<< refinedChains[c - tupChainClusters[s].start][cs].second.pos + refinedChainsLength[c - tupChainClusters[s].start]<< "\t" 
										<< read.length - refinedChains[c - tupChainClusters[s].start][cs].first.pos << "\t"
										<< refinedChains[c - tupChainClusters[s].start][cs].second.pos << "\t" 
										<< r << "\t"
										<< tupChainClusters[s].strand << endl;								
							}

						}

						// find small matches between fragments in gapChain
						int m, rg, gg;
						SetMatchAndGaps(curRefinedReadEnd, nextRefinedReadStart, curRefinedGenomeEnd, nextRefinedGenomeStart, m, rg, gg);
						if (m > 0) {
							Alignment betweenAnchorAlignment;
							if (opts.refineLevel & REF_DP) {						
								RefineSubstrings(strands[tupChainClusters[s].strand], curRefinedReadEnd, nextRefinedReadStart, genome.seqs[chromIndex], 
																 curRefinedGenomeEnd, nextRefinedGenomeStart, scoreMat, pathMat, betweenAnchorAlignment, opts);
								alignment->blocks.insert(alignment->blocks.end(), betweenAnchorAlignment.blocks.begin(), betweenAnchorAlignment.blocks.end());
								int b;
								for (b = 1; b < betweenAnchorAlignment.blocks.size(); b++) {
									assert(betweenAnchorAlignment.blocks[b-1].qPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].qPos);
									assert(betweenAnchorAlignment.blocks[b-1].tPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].tPos);						
								}
								betweenAnchorAlignment.blocks.clear();
							}
						}

						curRefinedReadEnd = refinedChains[c - tupChainClusters[s].start][cs].first.pos + refinedChainsLength[c - tupChainClusters[s].start];
						curRefinedGenomeEnd = refinedChains[c - tupChainClusters[s].start][cs].second.pos + refinedChainsLength[c - tupChainClusters[s].start];
						alignment->blocks.push_back(Block(nextRefinedReadStart, nextRefinedGenomeStart, curRefinedReadEnd - nextRefinedReadStart));

						if (alignment->blocks.size() > 1) {
							int last=alignment->blocks.size();
							assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].qPos);
							assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].tPos);
						}

					}
					// Add the last gap, or the only one if no refinements happened here.				 
					int match, readGap, genomeGap;
					SetMatchAndGaps(curRefinedReadEnd, nextReadStart, curRefinedGenomeEnd, nextGenomeStart, match, readGap, genomeGap);
					if (match > 0) {
						if (opts.refineLevel & REF_DP) {
							Alignment aln;
							assert(curRefinedReadEnd < read.length);
							assert(nextReadStart <read.length);
							assert(curRefinedGenomeEnd < genome.lengths[chromIndex]);
							assert(nextGenomeStart < genome.lengths[chromIndex]);
							RefineSubstrings(strands[tupChainClusters[s].strand], curRefinedReadEnd, nextReadStart, genome.seqs[chromIndex], 
															 curRefinedGenomeEnd, nextGenomeStart, scoreMat, pathMat, aln, opts);
							alignment->blocks.insert(alignment->blocks.end(), aln.blocks.begin(), aln.blocks.end());
							aln.blocks.clear();			
						}		
					}
				}
				alignment->blocks.push_back(Block(tupChain[tupChainClusters[s].end - 1].first.pos, tupChain[tupChainClusters[s].end - 1].second.pos, smallOpts.globalK));

				int nm=0;
				for(int b=0; b < alignment->blocks.size(); b++) {
					nm+= alignment->blocks[b].length;
				}
				alignment->nblocks = tupChainClusters[s].end - tupChainClusters[s].start;
				if (opts.dotPlot) {
					dotFile.close();
				}
			}
		}

		int cm = 0;
		for (int a = 0; a < alignments.size(); a++) {
			if (alignments[a].SegAlignment.size() != 0) {
				alignments[cm] = alignments[a];
				for (int b = 0; b < alignments[cm].SegAlignment.size(); b++) {
					alignments[cm].SegAlignment[b]->CalculateStatistics(alignments[cm].SegAlignment.size(), b);	
				}
				alignments[cm].SetBoundariesFromSegAlignmentAndnm(read);
				cm++;
			} 
		}
		alignments.resize(cm);
		//alignments.SetBoundariesFromSegAlignmentAndnm(read);

		sort(alignments.begin(), alignments.end(), SortAlignmentsByMatches());

		SimpleMapQV(alignments);



		if (opts.dotPlot) {
			stringstream outNameStrm;
			//outNameStrm << baseName + "." << a << ".alignment.dots";
			ofstream baseDots;
			//baseDots.open(outNameStrm.str().c_str());
		

			for (int a=0; a < min(opts.bestn, (int) alignments.size()); a++){

				outNameStrm << baseName + "." << a << ".alignment.dots";
				baseDots.open(outNameStrm.str().c_str());

				for (int s = 0; s < alignments[a].SegAlignment.size(); s++) {

					// for debug. TODO(Jingwen): delete this later fix the following code!!!
					for (int c = 0; c < alignments[a].SegAlignment[s]->blocks.size(); c++) {
						if (alignments[a].SegAlignment[s]->strand == 0) {
							baseDots << alignments[a].SegAlignment[s]->blocks[c].qPos << "\t" 
									 << alignments[a].SegAlignment[s]->blocks[c].tPos << "\t" 
									 << alignments[a].SegAlignment[s]->blocks[c].qPos + alignments[a].SegAlignment[s]->blocks[c].length << "\t" 
									 << alignments[a].SegAlignment[s]->blocks[c].tPos + alignments[a].SegAlignment[s]->blocks[c].length << "\t"
									 << a << "\t"
									 << s << "\t"
									 << alignments[a].SegAlignment[s]->strand << endl;							
						} 
						else {
							baseDots << read.length - alignments[a].SegAlignment[s]->blocks[c].qPos - alignments[a].SegAlignment[s]->blocks[c].length << "\t" 
									 << alignments[a].SegAlignment[s]->blocks[c].tPos + alignments[a].SegAlignment[s]->blocks[c].length << "\t" 
									 << read.length - alignments[a].SegAlignment[s]->blocks[c].qPos << "\t" 
									 << alignments[a].SegAlignment[s]->blocks[c].tPos << "\t"
									 << a << "\t"
									 << s << "\t"
									 << alignments[a].SegAlignment[s]->strand << endl;
						}
					}		
				}
			}
			
			baseDots.close();
		}


		//
		// Lock for output.
		//
		/*
		if (semaphore != NULL) {
			pthread_mutex_lock(semaphore);
		}
		*/

		for (int a=0; a < min(opts.bestn, (int) alignments.size()); a++){

			for (int s = 0; s < alignments[a].SegAlignment.size(); s++) {

				if (opts.printFormat == 'b') {
					alignments[a].SegAlignment[s]->PrintBed(*output);
				}
				else if (opts.printFormat == 's') {
					alignments[a].SegAlignment[s]->PrintSAM(*output, opts);
				}
				else if (opts.printFormat == 'p') {
					alignments[a].SegAlignment[s]->PrintPairwise(*output);
				}
			}
		}

		/*
		if (semaphore != NULL ) {
			pthread_mutex_unlock(semaphore);
		}
		*/

		//
		// Done with one read. Clean memory.
		//
		delete[] readRC;
		for (int a = 0; a < alignments.size(); a++) {
			for (int s = 0; s < alignments[a].SegAlignment.size(); s++) {
				delete alignments[a].SegAlignment[s];
			}
		}
		
		//read.Clear();

		/*
		// get the time for the program
		clock_t end = std::clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
		cerr << "Time: " << elapsed_secs << endl;
		*/
		if (alignments.size() > 0 ) {
			return 1;
		}
		else {
			return 0;
		}
	}
	return 0;
}

#endif
