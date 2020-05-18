#ifndef MAP_READ_H_
#define MAP_READ_H_
#include "MMIndex.h"
#include "Genome.h"
#include "Read.h"
#include "Options.h"
#include "CompareLists.h"
#include "Sorting.h"
#include "TupleOps.h"
#include "Clustering.h"
#include "AffineOneGapAlign.h"
#include "GlobalChain.h"
#include "TupleOps.h"
#include "SparseDP.h"
#include "Merge.h"
#include "MergeAnchors.h"
#include "Chain.h"
#include "overload.h"
#include "LinearExtend.h"
#include "SplitClusters.h"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cmath>	// std::log 
#include <sstream>
#include <thread>
#include <thread>
#include <climits>

using namespace std;

class Timing {
public:
	vector<long> ticks;
	vector<string> labels;
	int index;
	clock_t curTime;
	void Start() {
		index=0;
		curTime=clock()/1000;
	}
	void Tick(string label) {
		if (index +1 > ticks.size()) {
			ticks.push_back(0);
		}
		long prev=ticks[index];
		ticks[index] += clock()/1000 - curTime;
		if (prev > ticks[index]) {
			cerr << "Warning, tick has wrapped over long" << endl;
		}
		curTime=clock()/1000;
		if (index + 1 > labels.size()) {
			labels.push_back(label);
		}
		index+=1;
	}
	void Add(Timing &t) {
		//
		// Only add ticks the same size -- sometimes a thread 
		// will be created without adding ticks.
		//
		if (t.ticks.size() == ticks.size()) {
			for (int i =0; i < t.ticks.size(); i++) {
				ticks[i] += t.ticks[i];
			}
		}
	}
	void Summarize(const string &outFileName) {
		long total=0;
		ofstream outFile(outFileName.c_str());
		for (int i=0;i<ticks.size(); i++) {
			total+= ticks[i];
		}
		for (int i=0; i < ticks.size(); i++) {
			outFile << labels[i] << "\t" << ticks[i] << "\t" << ((double)ticks[i])/total << endl;
		}
		outFile.close();
	}
		
};


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

/*
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
*/

// void SimpleMapQV(AlignmentsOrder &alignmentsOrder, Read &read) {
// 	//
// 	// Store the index of primary alignments in vector<int> pry_index;
// 	//
// 	vector<int> pry_index;
// 	for (int sl = 0; sl < alignmentsOrder.size(); sl++) {
// 		if (!alignmentsOrder[sl].ISsecondary) { // This is primary alignment;
// 			pry_index.push_back(sl);
// 		}
// 	}
// 	//
// 	// compute mapqv for each alignment;
// 	//
// 	for (int pi = 0; pi < pry_index.size(); pi++) {
		
// 		int first = pry_index[pi];
// 		int last;
// 		if (pi == pry_index.size() - 1) last = alignmentsOrder.size();
// 		else last = pry_index[pi+1];

// 		if (last - first == 1) { // No secondary alignments
// 			alignmentsOrder[first].mapqv = 60;
// 			alignmentsOrder[first].SetMapqv();
// 		}
// 		else {
// 			// assign mapqv to primary alignment;
// 			int nmmdiff = alignmentsOrder[first].nmm - alignmentsOrder[first+1].nmm;
// 			if (nmmdiff < -20) { //-10
// 				alignmentsOrder[first].mapqv = 60;				
// 			}
// 			else if (nmmdiff >= -20 and nmmdiff <= 0) {
// 				if (alignmentsOrder[first].value < alignmentsOrder[first+1].value + 200) { 
// 					float denom = pow(0.4,20-abs(nmmdiff));  
// 					alignmentsOrder[first].mapqv = min(60, (int) (60 + 5*log10(denom)));// when nmmdiff=0, mapqv=20
// 				}
// 				else { alignmentsOrder[first].mapqv = 60;}				
// 			}
// 			else {
// 				//alignmentsOrder[first].mapqv = (int) 35*(1-alignmentsOrder[first+1].value/alignmentsOrder[first].value)*
// 				// (1 + alignmentsOrder[first].nm/(float)read.length);
// 				if (alignmentsOrder[first].value < alignmentsOrder[first+1].value + 200) { 
// 					float denom = pow(0.3,abs(nmmdiff));  //0.3 is the probability of nmmdiff are real snps. 
// 					int cpr=(int) (60 + 5*log10(denom)); //60 + 10*log10(denom) // when nmmdiff=20, mapqv=0
// 					if (cpr>=0) alignmentsOrder[first].mapqv = min(60, (int) (60 + 5*log10(denom)));
// 					else alignmentsOrder[first].mapqv = 0;
// 				}
// 				else { alignmentsOrder[first].mapqv = 60;}
// 			}
// 			alignmentsOrder[first].SetMapqv();

// 			assert(alignmentsOrder[first].mapqv <= 60);
// 			// assign mapqv to secondary alignments;
// 			for (int sd = first + 1; sd < last; sd++) {
// 				alignmentsOrder[sd].mapqv = (int) alignmentsOrder[first].mapqv/(2*(last-first-1));
// 				alignmentsOrder[sd].SetMapqv();
// 			}
// 		}
// 	}
// }


void SimpleMapQV(AlignmentsOrder &alignmentsOrder, Read &read) {
	//
	// Store the index of primary alignments in vector<int> pry_index;
	//
	vector<int> pry_index;
	for (int sl = 0; sl < alignmentsOrder.size(); sl++) {
		if (!alignmentsOrder[sl].ISsecondary) { // This is primary alignment;
			pry_index.push_back(sl);
		}
	}
	//
	// compute mapqv for each alignment;
	//
	for (int pi = 0; pi < pry_index.size(); pi++) {
		
		int first = pry_index[pi];
		int last;
		if (pi == pry_index.size() - 1) last = alignmentsOrder.size();
		else last = pry_index[pi+1];

		if (last - first == 1) { // No secondary alignments
			alignmentsOrder[first].mapqv = 60;
			alignmentsOrder[first].SetMapqv();
		}
		else {
			// assign mapqv to primary alignment;
			int nmmdiff = alignmentsOrder[first].nmm - alignmentsOrder[first+1].nmm;
			int nsmallgap = (alignmentsOrder[first].ndel + alignmentsOrder[first].nins) - 
								(alignmentsOrder[first+1].ndel + alignmentsOrder[first+1].nins);
			//cerr << "nmmdiff: " << nmmdiff << " nsmallgap: " << nsmallgap << endl;
			float denom_1 = 1, denom_2 = 1;
			if (nmmdiff <= 0 and nsmallgap <= 0) { // nmmdiff=0 and nsmallgap=0 ==> mapqv=52
				if (alignmentsOrder[first].value < alignmentsOrder[first+1].value + 300) { 
					if (nmmdiff > -20) denom_1 = 1;
					else denom_1 = pow(0.9,20-abs(nmmdiff));  
					if (nsmallgap > -20) denom_2 = 1;
					else  denom_2 = pow(0.9,20-abs(nsmallgap));  
					alignmentsOrder[first].mapqv = min(60, (int) (60 + 5*log10(denom_1) + 5*log10(denom_2)));
				}
				else { alignmentsOrder[first].mapqv = 60;}			
			}
			else if (nmmdiff <= 0 and nsmallgap >= 0) { // when nsmallgap=30 and nmmdiff=0; mapqv=0
				if (alignmentsOrder[first].value < alignmentsOrder[first+1].value + 300) { 
					if (nmmdiff < -20) denom_1 = 1;
					else denom_1 = pow(0.9,20-abs(nmmdiff));  
					denom_2 = pow(0.5, abs(nsmallgap));  
					int cpr=(int) (60 + 5*log10(denom_1) + 5*log10(denom_2)); //60 + 10*log10(denom) 
					if (cpr>=0) alignmentsOrder[first].mapqv = min(60, cpr);
					else alignmentsOrder[first].mapqv = 0;					
				}
				else { alignmentsOrder[first].mapqv = 60;}					

			}
			else if (nmmdiff >= 0 and nsmallgap <= 0) { // when nsmallgap=0 and nmmdiff=30, mapqv=0
				if (alignmentsOrder[first].value < alignmentsOrder[first+1].value + 300) { 
					denom_1 = pow(0.5, abs(nmmdiff));  
					if (nsmallgap < -20) denom_2 = 1;
					else  denom_2 = pow(0.9,20-abs(nsmallgap)); 
					int cpr=(int) (60 + 5*log10(denom_1) + 5*log10(denom_2)); //60 + 10*log10(denom) 
					if (cpr>=0) alignmentsOrder[first].mapqv = min(60, cpr);
					else alignmentsOrder[first].mapqv = 0;						
				}
				else { alignmentsOrder[first].mapqv = 60;}	
			}
			else { 
				if (alignmentsOrder[first].value < alignmentsOrder[first+1].value + 300) { 
					denom_1 = pow(0.5, abs(nmmdiff));  
					denom_2 = pow(0.5, abs(nsmallgap));  
					int cpr=(int) (60 + 5*log10(denom_1) + 5*log10(denom_2));
					if (cpr>=0) alignmentsOrder[first].mapqv = min(60, cpr);
					else alignmentsOrder[first].mapqv = 0;					
				}
				else { alignmentsOrder[first].mapqv = 60;}	
			}
			alignmentsOrder[first].SetMapqv();
			//
			// assign mapqv to secondary alignments;
			//
			assert(alignmentsOrder[first].mapqv <= 60);
			for (int sd = first + 1; sd < last; sd++) {
				alignmentsOrder[sd].mapqv = (int) alignmentsOrder[first].mapqv/(2*(last-first-1));
				alignmentsOrder[sd].SetMapqv();
			}
		}
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
void UpdateBoundaries(T &matches, GenomePos &qStart, GenomePos &qEnd, GenomePos &tStart, GenomePos &tEnd) {
	for (int i =0; i< matches.size(); i++) {
		qStart=min(qStart, matches[i].first.pos);
		qEnd=max(qEnd, matches[i].first.pos);
		tStart=min(tStart, matches[i].second.pos);
		tEnd=max(tEnd, matches[i].second.pos);
	}
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

// Revised: output strand for each matches
void SeparateMatchesByStrand(Read &read, Genome &genome, int k, vector<pair<GenomeTuple, GenomeTuple> > &allMatches, vector<bool> & strand) {
	//
	// A value of 0 implies forward strand match.
	//
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
}

//
// This function removes spurious anchors after SDP;
//
void RemoveSpuriousAnchors(FinalChain & chain, Options & opts) {
	int cs = 0, ce = 1;
	vector<bool> remove(chain.size(), false);

	while (ce < chain.size()) {

		int dist = (int)(chain[ce].first.pos + chain.length(ce) - chain[ce-1].first.pos); 
		if (chain[ce].second.pos + chain.length(ce) <= chain[ce - 1].second.pos) {
			dist = max(dist, (int)(chain[ce].second.pos - chain[ce].second.pos - chain.length(ce)));			
		}
		else if (chain[ce].second.pos >= chain[ce-1].second.pos + chain.length(ce-1)) {
			dist = max(dist, (int)(chain[ce].second.pos - chain[ce-1].second.pos - chain.length(ce-1)));			
		}
		else dist = opts.maxRemoveSpuriousAnchorsDist + 1;

		if (dist <= opts.maxRemoveSpuriousAnchorsDist) ce++;
		else {
			if (ce - cs < opts.minRemoveSpuriousAnchorsNum) {
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
			chain.chain[m] = chain.chain[s];
			m++;
		}
	}
	chain.resize(m);
}

//
// This function removes paired indels in the finalchain after SDP;
//
void RemovePairedIndels (FinalChain & chain, Options & opts) {

	if (chain.size() < 3) return;
	vector<bool> remove(chain.size(), false); // If remove[i] == true, then remove chain[i]

	for (int c = 1; c < chain.size() - 1; c++) {

		if (chain.strand(c-1) == chain.strand(c) and chain.strand(c) == chain.strand(c+1)) {
			int prevGap, nextGap;
			if (chain.strand(c) == 0) {
				prevGap = (int)(((long int)chain[c-1].second.pos - (long int)chain[c-1].first.pos) 
								- ((long int)chain[c].second.pos - (long int)chain[c].first.pos));
				nextGap = (int)(((long int)chain[c].second.pos - (long int)chain[c].first.pos) 
								- ((long int)chain[c+1].second.pos - (long int)chain[c+1].first.pos));
			}
			else {
				prevGap = (int)((long int)(chain[c-1].first.pos + chain[c-1].second.pos + chain.length(c-1) - 1) 
								- (long int)(chain[c].first.pos + chain[c].second.pos));
				nextGap = (int)((long int)(chain[c].first.pos + chain[c].second.pos + chain.length(c) - 1) 
								- (long int)(chain[c+1].first.pos + chain[c+1].second.pos));
			}
			if (sign(prevGap)!= sign(nextGap) and abs(prevGap) + abs(nextGap) > abs(prevGap + nextGap)
					and chain.length(c) < opts.minRemovePairedIndelsLength) { // the second condition filter out when prevGap == 0 or nextGap == 0
				remove[c] = true;
			}	
		}
	}

	int m = 0;
	for (int i = 0; i < chain.size(); i++) {
		if (remove[i] == false) {
			chain.chain[m] = chain.chain[i];
			m++;
		}
	}
	chain.resize(m);
}

//
// This function switches index in splitclusters back 
//
int
switchindex (vector<Cluster> & splitclusters, vector<Primary_chain> & Primary_chains, vector<Cluster> & clusters) {
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
				Primary_chains[p].chains[h].ch[c] = splitclusters[Primary_chains[p].chains[h].ch[c]].coarse;
			}
		}
	}

	// Change "vector<bool> link" accordingly
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			if (Primary_chains[p].chains[h].link.size() > 0) {
				vector<bool> rm(Primary_chains[p].chains[h].link.size(), 0);
				for (int c = 1; c < Primary_chains[p].chains[h].ch.size(); c++) {
					if (Primary_chains[p].chains[h].ch[c] == Primary_chains[p].chains[h].ch[c-1]) rm[c-1] = 1;
				}	

				int sm = 0; 
				for (int c = 0; c < Primary_chains[p].chains[h].link.size(); c++) {
					if (rm[c] == 0) {
						Primary_chains[p].chains[h].link[sm] = Primary_chains[p].chains[h].link[c];	
						sm++;					
					}
				}	
				Primary_chains[p].chains[h].link.resize(sm);						
			}
		}
	}

	//
	// Remove the dupplicates in consecutive group of elements
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
	  		vector<unsigned int>::iterator itp;
	  		itp = unique(Primary_chains[p].chains[h].ch.begin(), Primary_chains[p].chains[h].ch.end()); 
	  		int pdist = distance(Primary_chains[p].chains[h].ch.begin(), itp);
	  		Primary_chains[p].chains[h].ch.resize(pdist);			
		}
	}

	/*
	//Remove nonconsecutive dupplicates
	vector<int> count(clusters.size(), 0);
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			vector<bool> rm(Primary_chains[p].chains[h].ch.size(), 0);
			for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
				int ic = Primary_chains[p].chains[h].ch[c];
				if (count[ic] == 0) {
					count[ic] += 1;
				}
				else {
					rm[c] = 1;
				}
			}

			int sm = 0;
			for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
				if (rm[c] == 0){
					Primary_chains[p].chains[h].ch[sm] = Primary_chains[p].chains[h].ch[c];
					if (sm >= 1) Primary_chains[p].chains[h].link[sm-1] = Primary_chains[p].chains[h].link[c-1];
					sm++;
				}
				else {
					//Primary_chains[p].chains[h].link.erase(Primary_chains[p].chains[h].link.begin() + c - 1);
				}
			}
			Primary_chains[p].chains[h].ch.resize(sm);
			Primary_chains[p].chains[h].link.resize(sm-1);
		}
	}
	*/
	return 0;
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
		long long int maxDN, minDN;
		maxDN = (long long int) clusters[ph].matches[0].first.pos - (long long int) clusters[ph].matches[0].second.pos;
		minDN = maxDN;
		for (int db = 0; db < clusters[ph].matches.size(); db++) {
			maxDN = max(maxDN, (long long int)clusters[ph].matches[db].first.pos - (long long int)clusters[ph].matches[db].second.pos);
			minDN = min(minDN, (long long int)clusters[ph].matches[db].first.pos - (long long int)clusters[ph].matches[db].second.pos);
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

			int matchEnd = CartesianTargetUpperBound<GenomeTuple>(clusters[ph].matches.begin(), clusters[ph].matches.end(), 
																		genomeLocalIndexEnd);
			//
			// If there is no overlap with this cluster
			//
			if (matchStart >= clusters[ph].matches.size()) {
				continue;
			}
			GenomePos readStart = clusters[ph].matches[matchStart].first.pos;
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
				readEnd = clusters[ph].matches[matchEnd - 1].first.pos;
			}
			else {
				readEnd = clusters[ph].matches[matchStart].first.pos + opts.globalK;
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

				CompareLists<LocalTuple>(readIndex->minimizers.begin()+qStartBoundary, readIndex->minimizers.begin()+qEndBoundary, 
										glIndex.minimizers.begin()+ glIndex.tupleBoundaries[lsi], 
										glIndex.minimizers.begin()+ glIndex.tupleBoundaries[lsi+1], smallMatches, smallOpts);
				//
				// Add refined anchors if they fall into the diagonal band and cluster box
				//
				AppendValues<LocalPairs>(refinedclusters[ph].matches, smallMatches.begin(), smallMatches.end(), readSegmentStart, 
										genomeLocalIndexStart,clusters[ph].maxDiagNum, clusters[ph].minDiagNum, clusters[ph].qStart, 
										clusters[ph].qEnd, clusters[ph].tStart-chromOffset, clusters[ph].tEnd-chromOffset);
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


//
// This function find anchors btwn two adjacent Clusters;
//
int 			
RefineBtwnSpace(Cluster * cluster, Options & opts, Genome & genome, Read & read, char *strands[2], GenomePos qe, GenomePos qs, 
							GenomePos te, GenomePos ts, bool st, int &cur, GenomePos lrts=0, GenomePos lrlength=0) {

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
	// Decide the diagonal band for this space
	//
	long long int minDiagNum, maxDiagNum; 
	long long int diag1, diag2;
	diag1 = 0;
	diag2 = (long long int) (qe - qs) - (long long int) (te - (ts - lrts)); // scale diag1 and diag2 to the local coordinates
	minDiagNum = min(diag1, diag2) - 300; // used 50 before
	maxDiagNum = max(diag1, diag2) + 300; 
 
	//
	// Find matches in read and reference 
	//
	vector<GenomeTuple> EndReadTup, EndGenomeTup;
	GenomePairs EndPairs;
	StoreMinimizers<GenomeTuple, Tuple>(genome.seqs[ChromIndex]+(ts-lrts), te-ts+lrlength, opts.globalK, opts.globalW, EndGenomeTup, false);
	sort(EndGenomeTup.begin(), EndGenomeTup.end());
	StoreMinimizers<GenomeTuple, Tuple>(strands[st] + qs, qe - qs, opts.globalK, opts.globalW, EndReadTup, false);
	sort(EndReadTup.begin(), EndReadTup.end());
	CompareLists(EndReadTup.begin(), EndReadTup.end(), EndGenomeTup.begin(), EndGenomeTup.end(), EndPairs, opts, maxDiagNum, minDiagNum); // By passing maxDiagNum and minDiagNum, this function 																															// filters out anchors that are outside the diagonal band;

	for (int rm = 0; rm < EndPairs.size(); rm++) {
		EndPairs[rm].first.pos  += qs;
		EndPairs[rm].second.pos += ts-lrts;
		assert(EndPairs[rm].first.pos < read.length);
		if (st == 1) EndPairs[rm].first.pos  = read.length - EndPairs[rm].first.pos - opts.globalK;
		
	}	

	if (EndPairs.size() > 0) {
		cluster->matches.insert(cluster->matches.end(), EndPairs.begin(), EndPairs.end());  // TODO(Jingwen): Time consuming???????
		cluster->SetClusterBoundariesFromMatches(opts);
		cluster->refinespace = 1;
	}

	return 0;
}


//
// This function splits the chain if Clusters on the chain are mapped to different chromosomes or different locations (quite far, default: 100000) on the same chromosome;
//
void
SPLITChain(vector<Cluster> ExtendClusters, vector<SplitChain> & splitchains, vector<bool> & link, Options & opts) {

	//
	// Split the chain if Clusters are mapped to different chromosomes;
	//
	vector<int> Index; // Index stores the chromIndex of Clusters on chain
	vector<int> posIndex(ExtendClusters.size(), 0); // posIndex[i] means ExtendClusters[i] has chromIndex Index[posIndex[i]];

	//
	// Get Index and posIndex;
	//
	for (int im = 0; im < ExtendClusters.size(); im++) {

		int curChromIndex = ExtendClusters[im].chromIndex;
		if (Index.empty()) {
			Index.push_back(curChromIndex);
			posIndex[im] = Index.size() - 1;
		}
		else {
			int ex = 0;
			while (ex < Index.size()) {
				if (Index[ex] == curChromIndex) {
					posIndex[im] = ex;
					break;
				}
				ex++;
			}
			if (ex == Index.size()) {
				Index.push_back(curChromIndex);
				posIndex[im] = Index.size() - 1;
			}
		}
	}

	//
	// Get sptchain based on Index and posIndex;
	//
	vector<vector<unsigned int>> sptchain(Index.size()); //// TODO(Jingwen): make sure this initialization is right!
	for (int im = 0; im < posIndex.size(); im++) {
		sptchain[posIndex[im]].push_back(im);
	}

	//
	// split the chain if Clusters are mapped to different locations of the same chromosome;
	//
	int bf = 0;
	for (int im = 0; im < sptchain.size(); im++) {
		int in = 1;
		vector<unsigned int> onec; 
		vector<bool> lk;
		onec.push_back(sptchain[im][0]);

		while (in < sptchain[im].size()) {

			int prev = sptchain[im][in-1];
			int cur = sptchain[im][in];

			if (ExtendClusters[cur].tStart > ExtendClusters[prev].tEnd + opts.splitdist
				or ExtendClusters[cur].tEnd + opts.splitdist < ExtendClusters[prev].tStart) {
				
				splitchains.push_back(SplitChain(onec, lk));
				onec.clear();
				lk.clear();
				onec.push_back(cur);
			}
			else {
				onec.push_back(cur);
				lk.push_back(link[bf+in-1]);
			}
			in++;
		}
		splitchains.push_back(SplitChain(onec, lk));
		bf += sptchain[im].size() - 1;
	}
}

//
// This function seperates finalchain by different strands; 
//
void
SeperateChainByStrand(FinalChain & finalchain, vector<vector<int>> & finalSeperateChain, const vector<Cluster> & ExtendClusters) {

	vector<int> sep;
	int fl = 0;
	sep.push_back(fl);
	while (fl < finalchain.chain.size() - 1) {

		if (finalchain.strand(fl) == finalchain.strand(fl+1)) {
			fl++;
		}
		else {
			sep.push_back(fl+1);
			finalSeperateChain.push_back(sep);
			sep.clear();
			sep.push_back(fl+1);
			fl++;
		}
	}
	if (fl == finalchain.chain.size() - 1) {
		sep.push_back(fl+1);
		finalSeperateChain.push_back(sep);
	}	
}

/*
//
// This function sorts the SeperateChain by the number of anchors in the descending order;
//
template <typename T> 
class SortChainByAnchorNumOp {
public: 
	int operator() (const T & a, const T & b) {
		return a.back() - a[0] > b.back() - b[0];
	}
};
*/

int LargestSplitChain(vector<SplitChain> &splitchains) {
	int maxi = 0;
	for (int mi = 1; mi < splitchains.size(); mi++) {
		if (splitchains[mi].sptc.size() > splitchains[maxi].sptc.size()) {
			maxi = mi;
		}
	}
	return maxi;
}


int LargestFinalSeperateChain(vector<vector<int>> &finalSeperateChain) {
	int maxi = 0;
	for (int mi = 1; mi < finalSeperateChain.size(); mi++) {
		if (finalSeperateChain[mi].back() - finalSeperateChain[mi][0] > finalSeperateChain[maxi].back() - finalSeperateChain[maxi][0]){
			maxi = mi;
		}
	}
	return maxi;	
}

//
// This function creates the alignment for a "segment"(one big chunk) on the chain;
//
void
RefinedAlignmentbtwnAnchors(int & cur, int & next, int & str, int & chromIndex, FinalChain & finalchain, Alignment * alignment, Read & read, Genome & genome, char *strands[2], 
				vector<int> & scoreMat, vector<Arrow> & pathMat, Options & tinyOpts, GenomePos & genomeThres) {

	if (str == 0) alignment->blocks.push_back(Block(finalchain[cur].first.pos, finalchain[cur].second.pos, finalchain.length(cur)));
	else alignment->blocks.push_back(Block(read.length - finalchain[cur].first.pos - finalchain.length(cur), finalchain[cur].second.pos, finalchain.length(cur)));

	// This is for debugging
	if (alignment->blocks.size() > 1) {
		int last=alignment->blocks.size();
		assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].qPos);
		assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].tPos);
	}

	//
	// Refine the alignment in the space between two adjacent anchors;
	GenomePos curGenomeEnd, curReadEnd, nextGenomeStart, nextReadStart;
	if (str == 0) {
		curReadEnd = finalchain[cur].first.pos + finalchain.length(cur);
		nextReadStart = finalchain[next].first.pos;
		curGenomeEnd = finalchain[cur].second.pos + finalchain.length(cur);
		nextGenomeStart = finalchain[next].second.pos; 
	}
	else { // flip the space if it is reversed stranded
		curReadEnd = read.length - finalchain[cur].first.pos;
		nextReadStart = read.length - finalchain[next].first.pos - finalchain.length(next);
		curGenomeEnd = finalchain[cur].second.pos + finalchain.length(cur);
		nextGenomeStart = finalchain[next].second.pos;
	}
	//assert(curGenomeEnd <= nextGenomeStart);
	if (curGenomeEnd <= nextGenomeStart) {
		//
		// find small matches between fragments in gapChain
		int m, rg, gg;
		SetMatchAndGaps(curReadEnd, nextReadStart, curGenomeEnd, nextGenomeStart, m, rg, gg);
		if (m > 0) {
			Alignment betweenAnchorAlignment;
			if (tinyOpts.refineLevel & REF_DP) {						
				RefineSubstrings(strands[str], curReadEnd, nextReadStart, genome.seqs[chromIndex], 
												 curGenomeEnd, nextGenomeStart, scoreMat, pathMat, betweenAnchorAlignment, tinyOpts);
				int b;
				for (b = 1; b < betweenAnchorAlignment.blocks.size(); b++) {
					assert(betweenAnchorAlignment.blocks[b-1].qPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].qPos);
					assert(betweenAnchorAlignment.blocks[b-1].tPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].tPos);						
				}

				alignment->blocks.insert(alignment->blocks.end(), betweenAnchorAlignment.blocks.begin(), betweenAnchorAlignment.blocks.end());
				betweenAnchorAlignment.blocks.clear();
			}
		}	
		genomeThres = 0;		
	}
	else genomeThres = curGenomeEnd;
}

int 
PassgenomeThres(int cur, GenomePos & genomeThres, FinalChain & finalchain) {
	if (genomeThres != 0 and finalchain[cur].second.pos < genomeThres) return 1;
	else return 0;
}

int MapRead(const vector<float> & LookUpTable, Read &read, Genome &genome, vector<GenomeTuple> &genomemm, LocalIndex &glIndex, Options &opts, 
						ostream *output, ostream *svsigstrm, Timing &timing, pthread_mutex_t *semaphore=NULL) {
	string baseName = read.name;

	for (int i=0; i < baseName.size(); i++) {	
		if (baseName[i] == '/') baseName[i] = '_';	
		if (baseName[i] == '|') baseName[i] = '_';
	}


	vector<GenomeTuple> readmm; // readmm stores minimizers
	vector<pair<GenomeTuple, GenomeTuple> > allMatches, forMatches, revMatches, matches;
	timing.Start();

	if (opts.storeAll) {
		Options allOpts = opts;
		allOpts.globalW=1;
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, allOpts.globalK, allOpts.globalW, readmm);			
	}
	else {
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, opts.globalK, opts.globalW, readmm);
	}
	sort(readmm.begin(), readmm.end()); //sort kmers in readmm(minimizers)

	timing.Tick("Store minimizers");
	//
	// Add matches between the read and the genome.
	//
	CompareLists(readmm, genomemm, allMatches, opts);
	DiagonalSort<GenomeTuple>(allMatches); // sort fragments in allMatches by forward diagonal, then by first.pos(read)
	timing.Tick("Sort minimizers");

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
	int minDiagCluster = 0; // This parameter will be set inside function CleanOffDiagonal, according to anchors density
	CleanOffDiagonal(forMatches, opts, minDiagCluster);


	vector<Cluster> clusters;
	vector<Cluster> roughClusters;
	int forwardStrand=0;
	// maxDiag must be large enough for the following function "StoreDiagonalClusters". 
	// Otherwise fragments on the same line (line with little curve) might end up in several clusters[i], instead of one clusters[i]

	// The strategy we are taking: 
	// first let maxDiag be a small number but not too small(like 500), which also alleviate the cases where anchors are on a curvy line.
	// Then break the clusters[i] into two if any two anchors are father than maxGap. 
	StoreDiagonalClusters(forMatches, roughClusters, opts, 0, forMatches.size(), true, false, forwardStrand); // rough == true means only storing "start and end" in every clusters[i]

	if (opts.dotPlot) {
		ofstream clust("for-matches.dots");
		for (int m=0; m < forMatches.size(); m++) {
			clust << forMatches[m].first.pos << "\t" << forMatches[m].second.pos << "\t" << opts.globalK + forMatches[m].first.pos << "\t"
					<< forMatches[m].second.pos + opts.globalK << "\t" << m << "\t0"<<endl;
		}
		clust.close();
		ofstream wclust("roughClusters-matches.dots");
		for (int m=0; m < roughClusters.size(); m++) {
			int rci = genome.header.Find(roughClusters[m].tStart);
			for (int c = roughClusters[m].start; c < roughClusters[m].end; ++c) {
				wclust << forMatches[c].first.pos << "\t" << forMatches[c].second.pos << "\t" << opts.globalK + forMatches[c].first.pos << "\t"
					<< forMatches[c].second.pos + opts.globalK << "\t" << m << "\t" << genome.header.names[rci]<< "\t0"<<endl;				
			}
		}
		wclust.close();

	}

	interval_map<GenomePos, int> xIntv;
	interval_map<GenomePos, int> yIntv;

	for (int c = 0; c < roughClusters.size(); c++) {
		CartesianSort(forMatches, roughClusters[c].start, roughClusters[c].end);
		//		StoreDiagonalClusters(forMatches, clusters, opts, roughclusters[c].start, roughclusters[c].end, false, false, forwardStrand);
		int rci = genome.header.Find(roughClusters[c].tStart);
		//		cout << "roughClusters: " << c << "\t" << roughClusters[c].start << "\t" << roughClusters[c].end << endl;
		StoreFineClusters(forMatches, clusters, opts, roughClusters[c].start, roughClusters[c].end, read.length, xIntv, yIntv, forwardStrand, c);
	}


	AntiDiagonalSort<GenomeTuple>(revMatches, genome.GetSize());
	minDiagCluster = 0; // This parameter will be set inside function CleanOffDiagonal, according to anchors density
	CleanOffDiagonal(revMatches, opts, minDiagCluster, 1);


	vector<Cluster> revroughClusters;
	int reverseStrand=1;
	StoreDiagonalClusters(revMatches, revroughClusters, opts, 0, revMatches.size(), true, false, reverseStrand);
	vector<int> ci(revroughClusters.size());
	for (int c = 0 ; c < revroughClusters.size(); c++ ) {
		
		for (int cii=revroughClusters[c].start; cii< revroughClusters[c].end; ++cii) {
			ci.push_back(c);
		}
	}

	if (opts.dotPlot) {
		ofstream rclust("rev-matches.dots");
		for (int m=0; m < revMatches.size(); m++) {			
			rclust << revMatches[m].first.pos << "\t" << revMatches[m].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[m].first.pos << "\t"
					<< revMatches[m].second.pos << "\t" << ci[m] << "\t0"<<endl;
		}
		rclust.close();


		ofstream revclust("revroughClusters-matches.dots");
		for (int m=0; m < revroughClusters.size(); m++) {
			int rci = genome.header.Find(revroughClusters[m].tStart);
			for (int c = revroughClusters[m].start; c < revroughClusters[m].end; ++c) {
				revclust << revMatches[c].first.pos << "\t" << revMatches[c].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[c].first.pos << "\t"
					 << revMatches[c].second.pos<< "\t" << m << "\t" << genome.header.names[rci]<<"\t0"<<endl;				
			}
		}
		revclust.close();
	}

	
	for (int c = 0; c < revroughClusters.size(); c++) {
		CartesianSort(revMatches, revroughClusters[c].start, revroughClusters[c].end);
		//StoreDiagonalClusters(revMatches, clusters, opts, revroughClusters[c].start, revroughClusters[c].end, false, false, reverseStrand);
		int rci = genome.header.Find(revroughClusters[c].tStart);
		//		cout << "revroughClusters: " << c << "\t" << revroughClusters[c].end - revroughClusters[c].start << "\tchr: " << genome.header.names[rci] << " " << revroughClusters[c].tStart << "\t" << revroughClusters[c].tEnd << "\t" << revroughClusters[c].qStart << "\t" << revroughClusters[c].qEnd << endl;
		StoreFineClusters(revMatches, clusters, opts, revroughClusters[c].start, revroughClusters[c].end, read.length, xIntv, yIntv, reverseStrand, c);
	}
	/*
	for (int c=0; c < clusters.size(); c++) {
		cout << "cluster\t" << c << "\t" << clusters[c].qStart << "\t" << clusters[c].qEnd << "\t" << clusters[c].tStart << "\t" << clusters[c].tEnd << "\t" << clusters[c].matches.size() << endl;
	}
	*/
	/*	
	for (int c=0; c < clusters.size(); c++) { 
		CartesianSort(clusters[c].matches, 0, clusters[c].matches.size());
	}
	*/

	//
	// Split clusters on x and y coordinates, vector<Cluster> splitclusters, add a member for each splitcluster to specify the original cluster it comes from
	//
	// INPUT: vector<Cluster> clusters   OUTPUT: vector<Cluster> splitclusters with member--"coarse" specify the index of the original cluster splitcluster comes from
	if (clusters.size() == 0) {	
		return 0; // This read cannot be mapped to the genome; 
	}


	if (opts.dotPlot) {
		ofstream clust("clusters-coarse.tab");
		for (int m = 0; m < clusters.size(); m++) {
			if (clusters[m].strand == 0) {
				clust << clusters[m].qStart << "\t"
							<< clusters[m].tStart << "\t"
							<< clusters[m].qEnd   << "\t"
							<< clusters[m].tEnd   << "\t"
							<< m << "\t"
							<< clusters[m].strand << "\t"
							<< clusters[m].outerCluster << "\t"
							<< clusters[m].matches.size() << endl;
			}
			else {
				clust << clusters[m].qStart << "\t"
							<< clusters[m].tEnd   << "\t"
							<< clusters[m].qEnd   << "\t"
							<< clusters[m].tStart << "\t"
							<< m << "\t"
							<< clusters[m].strand << "\t"
							<< clusters[m].outerCluster << "\t"
							<< clusters[m].matches.size() << endl;
			}
		}
		clust.close();
	}
	//
	// remove Cluster that firstChromIndex != lastChromIndex; or remove Cluster of 0 matches;
	//
	vector<bool> RV(clusters.size(), 0);
	for (int m = 0; m < clusters.size(); m++) {
		if (CHROMIndex(clusters[m], genome) == 1) {
			RV[m] = 1;
		}
	}
	int vm = 0;
	for (int m = 0; m < RV.size(); m++) {
		if (RV[m] == 0) {
			clusters[vm] = clusters[m];
			vm++;
		}
	}
	clusters.resize(vm);
	if (clusters.size() == 0) return 0; // This read cannot be mapped to the genome; 

	vector<Cluster> splitclusters;
	SplitClusters(clusters, splitclusters);
	DecideSplitClustersValue(clusters, splitclusters, opts);

	timing.Tick("Preprocess");

	if (opts.dotPlot) {
		ofstream clust("clusters.tab");
		for (int m = 0; m < clusters.size(); m++) {
			for (int n = 0; n < clusters[m].matches.size(); n++) {
				if (clusters[m].strand == 0) {
					clust << clusters[m].matches[n].first.pos << "\t"
						  << clusters[m].matches[n].second.pos << "\t"
						  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].second.pos + opts.globalK << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].strand << endl;
				}
				else {
					clust << clusters[m].matches[n].first.pos << "\t"
						  << clusters[m].matches[n].second.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].second.pos << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].strand << endl;
				}				
			}
		}
		clust.close();
	}


	if (opts.dotPlot) {
		ofstream clust("splitclusters-coarse.tab");
		for (int m = 0; m < splitclusters.size(); m++) {
			if (splitclusters[m].strand == 0) {
				clust << splitclusters[m].qStart << "\t" 
					  << splitclusters[m].tStart << "\t"
					  << splitclusters[m].qEnd   << "\t"
					  << splitclusters[m].tEnd   << "\t"
					  << m << "\t"
					  << splitclusters[m].coarse << "\t"
					  << splitclusters[m].strand << endl;
			}
			else {
				clust << splitclusters[m].qStart << "\t" 
					  << splitclusters[m].tEnd << "\t"
					  << splitclusters[m].qEnd   << "\t"
					  << splitclusters[m].tStart   << "\t"
					  << m << "\t"
					  << splitclusters[m].coarse << "\t"
					  << splitclusters[m].strand << endl;
			}
		}
		clust.close();
	}


	if (opts.dotPlot) {
		ofstream clust("splitclusters-decideval.tab");
		for (int m = 0; m < splitclusters.size(); m++) {
			if (splitclusters[m].Val !=0) {
				if (splitclusters[m].strand == 0) {
					clust << splitclusters[m].qStart << "\t" 
						  << splitclusters[m].tStart << "\t"
						  << splitclusters[m].qEnd   << "\t"
						  << splitclusters[m].tEnd   << "\t"
						  << m << "\t"
						  << splitclusters[m].coarse << "\t"
						  << splitclusters[m].strand << endl;
				}
				else {
					clust << splitclusters[m].qStart << "\t" 
						  << splitclusters[m].tEnd << "\t"
						  << splitclusters[m].qEnd   << "\t"
						  << splitclusters[m].tStart   << "\t"
						  << m << "\t"
						  << splitclusters[m].coarse << "\t"
						  << splitclusters[m].strand << endl;
				}				
			}
		}
		clust.close();
	}
	
	//
	// Apply SDP on splitclusters. Based on the chain, clean clusters to make it only contain clusters that are on the chain.   --- vector<Cluster> clusters
	// class: chains: vector<chain> chain: vector<vector<int>>     Need parameters: PrimaryAlgnNum, SecondaryAlnNum
	// NOTICE: chains in Primary_chains do not overlap on Cluster
	//

	////// TODO(Jingwen): customize a rate fro SparseDP
	vector<Primary_chain> Primary_chains, Primary_chains2;
	float rate = 4;
	if (splitclusters.size() > 10) { // mapping to repetitive region
		rate = 1; 
	}
	SparseDP (splitclusters, Primary_chains, opts, LookUpTable, read, rate);
	SparseDP (splitclusters, Primary_chains2, opts, LookUpTable, read, rate);

	if (opts.dotPlot) {
		ofstream clust("Chains.tab");
		for (int p = 0; p < Primary_chains.size(); p++) {
			for (int h = 0; h < Primary_chains[p].chains.size(); h++){
				for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
					int ph = Primary_chains[p].chains[h].ch[c];
					if (splitclusters[ph].strand == 0) {
						clust << splitclusters[ph].qStart << "\t" 
							  << splitclusters[ph].tStart << "\t"
							  << splitclusters[ph].qEnd   << "\t"
							  << splitclusters[ph].tEnd   << "\t"
							  << p << "\t"
							  << h << "\t"
							  << Primary_chains[p].chains[h].ch[c] << "\t"
							  << splitclusters[ph].strand << endl;
					} 
					else {
						clust << splitclusters[ph].qStart << "\t" 
							  << splitclusters[ph].tEnd << "\t"
							  << splitclusters[ph].qEnd   << "\t"
							  << splitclusters[ph].tStart  << "\t"
							  << p << "\t"
							  << h << "\t"
							  << Primary_chains[p].chains[h].ch[c] << "\t"
							  << splitclusters[ph].strand << endl;
					}
				}
			}
		}
		clust.close();
	}	

	switchindex(splitclusters, Primary_chains, clusters);
	timing.Tick("Sparse DP - clusters");
	//
	// Remove Clusters in "clusters" that are not on the chains;
	//
	int ChainNum = 0;
	for (int p = 0; p < Primary_chains.size(); p++) {
		ChainNum += Primary_chains[p].chains.size();
	}

	vector<bool> Remove(clusters.size(), 1);
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++){
			for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
				Remove[Primary_chains[p].chains[h].ch[c]] = 0;
			}
		}
	}

	int lm = 0;
	for (int s = 0; s < clusters.size(); s++) {
		if (Remove[s] == 0) {
			clusters[lm] = clusters[s];
			lm++;
		}
	}
	clusters.resize(lm);	

	//
	// Change the index stored in Primary_chains, since we remove some Clusters in "clusters";
	//
	vector<int> NumOfZeros(Remove.size(), 0);
	int num = 0;
	for (int s = 0; s < Remove.size(); s++) {
		if (Remove[s] == 0) {
			num++;
			NumOfZeros[s] = num;
		}
	}

	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++){
			for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
				Primary_chains[p].chains[h].ch[c] = NumOfZeros[Primary_chains[p].chains[h].ch[c]] - 1;
			}
		}
	}


	if (opts.dotPlot) {
		ofstream clust("CoarseChains.tab");

		for (int p = 0; p < Primary_chains.size(); p++) {
			for (int h = 0; h < Primary_chains[p].chains.size(); h++){
				//cerr << "p: " << p << " h: " << h << " chr: " << genome.header.names[genome.header.Find(Primary_chains[p].chains[h].tStart)] << 
				 // " value: " << Primary_chains[p].chains[h].value << " tStart: " <<  Primary_chains[p].chains[h].tStart << endl;
				for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
					int ph = Primary_chains[p].chains[h].ch[c];
					if (clusters[ph].strand == 0) {
						clust << clusters[ph].qStart << "\t" 
							  << clusters[ph].tStart << "\t"
							  << clusters[ph].qEnd   << "\t"
							  << clusters[ph].tEnd   << "\t"
							  << p << "\t"
							  << h << "\t"
							  << Primary_chains[p].chains[h].ch[c] << "\t"
							  << clusters[ph].strand << endl;
					} 
					else {
						clust << clusters[ph].qStart << "\t" 
							  << clusters[ph].tEnd << "\t"
							  << clusters[ph].qEnd   << "\t"
							  << clusters[ph].tStart  << "\t"
							  << p << "\t"
							  << h << "\t"
							  << Primary_chains[p].chains[h].ch[c] << "\t"
							  << clusters[ph].strand << endl;
					}
				}
			}
		}
		clust.close();
	}	

	if (opts.dotPlot) {
		ofstream clust("clusters_1stSDP.tab");
		for (int m = 0; m < clusters.size(); m++) {
			for (int n = 0; n < clusters[m].matches.size(); n++) {
				if (clusters[m].strand == 0) {
					clust << clusters[m].matches[n].first.pos << "\t"
						  << clusters[m].matches[n].second.pos << "\t"
						  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].second.pos + opts.globalK << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].strand << endl;
				}
				else {
					clust << clusters[m].matches[n].first.pos << "\t"
						  << clusters[m].matches[n].second.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].second.pos << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].strand << endl;
				}				
			}
		}
		clust.close();
	}

	//
	// Add pointers to seq that make code more readable.
	//
	char *readRC;
	CreateRC(read.seq, read.length, readRC);
	char *strands[2] = { read.seq, readRC };


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
	smallOpts.coefficient = 9; // used to be 9
	Options tinyOpts = smallOpts;
	tinyOpts.globalMaxFreq=3;
	tinyOpts.maxDiag=5;
	tinyOpts.minDiagCluster=2;
	tinyOpts.globalK=smallOpts.globalK-3;


	// Decide whether the number of anchors in each Cluster is enough to skip refining step;
	bool sparse = 0;
	for (int p = 0; p < clusters.size(); p++) {
		if (((float)(clusters[p].matches.size())/(clusters[p].qEnd - clusters[p].qStart)) < 0.005) sparse = 1;
		else sparse = 0;
	}

	//
	// Refining each cluster in "clusters" needed if CLR reads are aligned OR CCS reads with very few anchors. Otherwise, skip this step
	// After this step, the t coordinates in clusters and refinedclusters have been substract chromOffSet. 
	//
	vector<Cluster> refinedclusters(clusters.size());
 	vector<Cluster*> RefinedClusters(clusters.size());

	if (opts.HighlyAccurate == false or (opts.HighlyAccurate == true and sparse == 1) ) {
			
		smallOpts.globalK=glIndex.k;
		smallOpts.globalW=glIndex.w;
		smallOpts.coefficient = 12; // used to be 9
		smallOpts.globalMaxFreq=6;
		smallOpts.cleanMaxDiag=10;// used to be 25
		smallOpts.maxDiag=50;
		smallOpts.maxGapBtwnAnchors=100; // used to be 200 // 200 seems a little bit large
		smallOpts.minDiagCluster=3; // used to be 3

		REFINEclusters(clusters, refinedclusters, genome, read,  glIndex, localIndexes, smallOpts, opts);
		// refinedclusters have GenomePos, chromIndex, coarse, matches, strand, refinespace;
		for (int s = 0; s < clusters.size(); s++) {
			RefinedClusters[s] = &refinedclusters[s];
		}
		clusters.clear();
	}
	else {
		//// TODO(Jingwen): write a function to determine the chromIndex and subtract chromOffSet from t coord.
		for (int s = 0; s < clusters.size(); s++) {	

			// Determine the chromIndex of each Cluster;
			int pass = CHROMIndex(clusters[s], genome);
			if (pass == 1) {
				clusters[s].matches.clear();
				continue;				
			}

			// Subtract chromOffSet from t coord.
			GenomePos chromOffset = genome.header.pos[clusters[s].chromIndex];
			for (int m = 0; m < clusters[s].matches.size(); m++) {
				clusters[s].matches[m].second.pos -= chromOffset;
			}
			clusters[s].tStart -= chromOffset;
			clusters[s].tEnd -= chromOffset;
			RefinedClusters[s] = &clusters[s];
		}
	}

	if (opts.dotPlot) {
		ofstream clust("RefinedClusters.tab");

		for (int p = 0; p < RefinedClusters.size(); p++) {

			for (int h = 0; h < RefinedClusters[p]->matches.size(); h++) {

				if (RefinedClusters[p]->strand == 0) {
					clust << RefinedClusters[p]->matches[h].first.pos << "\t"
						  << RefinedClusters[p]->matches[h].second.pos << "\t"
						  << RefinedClusters[p]->matches[h].first.pos + smallOpts.globalK << "\t"
						  << RefinedClusters[p]->matches[h].second.pos + smallOpts.globalK << "\t"
						  << p << "\t"
						  << genome.header.names[RefinedClusters[p]->chromIndex] <<"\t"
						  << RefinedClusters[p]->strand << endl;
				}
				else {
					clust << RefinedClusters[p]->matches[h].first.pos << "\t"
						  << RefinedClusters[p]->matches[h].second.pos + smallOpts.globalK << "\t"
						  << RefinedClusters[p]->matches[h].first.pos + smallOpts.globalK << "\t"
						  << RefinedClusters[p]->matches[h].second.pos<< "\t"
						  << p << "\t"
						  << genome.header.names[RefinedClusters[p]->chromIndex] <<"\t"
						  << RefinedClusters[p]->strand << endl;					
				}
			}
		}
		clust.close();
	}	

	timing.Tick("Refine");
	int SizeRefinedClusters = 0;
	for (int p = 0; p < RefinedClusters.size(); p++) {
		SizeRefinedClusters += RefinedClusters[p]->matches.size();
	}
	if (SizeRefinedClusters == 0) return 0; // This read is not mapped!

	
	//
	// Remove RefinedCluster without matches;
	//
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			int cp = 0;
			for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
				if (RefinedClusters[Primary_chains[p].chains[h].ch[c]]->matches.size() != 0) {
					Primary_chains[p].chains[h].ch[cp] = Primary_chains[p].chains[h].ch[c];
					cp++;
				}
			}
			Primary_chains[p].chains[h].ch.resize(cp);
		}	
	}

	//
	// For each chain, check the two ends and spaces between adjacent clusters. If the spaces are too wide, go to find anchors in the banded region.
	// For each chain, we have vector<Cluster> btwnClusters to store anchors;
	//
	vector<SegAlignmentGroup> alignments;
	AlignmentsOrder alignmentsOrder(&alignments);

	for (int p = 0; p < Primary_chains.size(); p++) {
		
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			
			if (Primary_chains[p].chains[h].ch.size() == 0) continue;
			alignments.resize(alignments.size() + 1);
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
				qs = RefinedClusters[cur]->qEnd; 
				qe = RefinedClusters[prev]->qStart;

				if (RefinedClusters[cur]->tEnd <= RefinedClusters[prev]->tStart) {
					ts = RefinedClusters[cur]->tEnd;
					te = RefinedClusters[prev]->tStart;
					//st = 0;
					st = 1;
				}
				else if (RefinedClusters[cur]->tStart > RefinedClusters[prev]->tEnd) {
					ts = RefinedClusters[prev]->tEnd;
					te = RefinedClusters[cur]->tStart;
					//st = 1;
					st = 0;
				}
				else {
					c++;
					continue; // No need to refine the space!					
				}

				if (qe > qs and te > ts) {
					SpaceLength = max(qe - qs, te - ts);
					if (SpaceLength > 1000 and SpaceLength < 6000 and RefinedClusters[cur]->chromIndex == RefinedClusters[prev]->chromIndex) {
						// btwnClusters have GenomePos, st, matches, coarse
						// This function also set the "coarse" flag for RefinedClusters[cur]
						RefineBtwnSpace(RefinedClusters[cur], smallOpts, genome, read, strands, qe, qs, te, ts, st, cur);
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
				if (SpaceLength < 10000 and te+500 < genome.lengths[RefinedClusters[rh]->chromIndex]) { // used (1000, 6000)
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
					RefineBtwnSpace(RefinedClusters[rh], smallOpts, genome, read, strands, qe, qs, te, ts, st, rh, lrts, lrlength);
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
			//cerr << "left  p: " << p << " h: " << h << " qs: " << qs << " qe: " << qe << " ts: " << ts << " te: " << te << endl;
			if (qe > qs and te > ts) {
				SpaceLength = max(qe - qs, te - ts);
				if (SpaceLength < 10000 and te+500 < genome.lengths[RefinedClusters[lh]->chromIndex]) { // used (1000, 6000)
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
					RefineBtwnSpace(RefinedClusters[lh], smallOpts, genome, read, strands, qe, qs, te, ts, st, lh, lrts, lrlength);
				}				
			}
			//
			// Do linear extension for each anchors and avoid overlapping locations;
			// INPUT: RefinedClusters; OUTPUT: ExtendClusters;
			// NOTICE: ExtendClusters have members: strand, matches, matchesLengths, GenomePos, chromIndex;
			//
			vector<Cluster> ExtendClusters(Primary_chains[p].chains[h].ch.size());
			LinearExtend(RefinedClusters, ExtendClusters, Primary_chains[p].chains[h].ch, smallOpts, genome, read);

			int SizeExtendClusters = 0;
			for (int ep = 0; ep < ExtendClusters.size(); ep++) {
				SizeExtendClusters += ExtendClusters[ep].matches.size();
			}	

			if (opts.dotPlot) {
				ofstream clust("ExtendClusters.tab", std::ofstream::app);
				for (int ep = 0; ep < ExtendClusters.size(); ep++) {
					for (int eh = 0; eh < ExtendClusters[ep].matches.size(); eh++) {	
						if (ExtendClusters[ep].strand == 0) {
							clust << ExtendClusters[ep].matches[eh].first.pos << "\t"
								  << ExtendClusters[ep].matches[eh].second.pos << "\t"
								  << ExtendClusters[ep].matches[eh].first.pos + ExtendClusters[ep].matchesLengths[eh] << "\t"
								  << ExtendClusters[ep].matches[eh].second.pos + ExtendClusters[ep].matchesLengths[eh] << "\t"
								  << h << "\t"
								  << genome.header.names[ExtendClusters[ep].chromIndex]<< "\t"
								  << ExtendClusters[ep].strand << endl;
						}
						else {
							clust << ExtendClusters[ep].matches[eh].first.pos << "\t"
								  << ExtendClusters[ep].matches[eh].second.pos + ExtendClusters[ep].matchesLengths[eh] << "\t"
								  << ExtendClusters[ep].matches[eh].first.pos + ExtendClusters[ep].matchesLengths[eh] << "\t"
								  << ExtendClusters[ep].matches[eh].second.pos<< "\t"
								  << h << "\t"
								  << genome.header.names[ExtendClusters[ep].chromIndex]<< "\t"
								  << ExtendClusters[ep].strand << endl;					
						}
					}
				}
				clust.close();
			}	

			assert(SizeRefinedClusters != 0);
			/*
			cerr << "SizeRefinedClusters: " << SizeRefinedClusters << "   SizeExtendClusters: " << SizeExtendClusters
				   << "  read.name:"<< read.name <<  endl;
			cerr << "LinearExtend efficiency: " << (float)SizeExtendClusters/(float)SizeRefinedClusters << endl;
			*/

			//
			// Split the chain Primary_chains[p].chains[h] if clusters are aligned to different chromosomes; 
			// SplitAlignment is class that vector<* vector<Cluster>>
			// INPUT: vector<Cluster> ExtendClusters; OUTPUT:  vector<vector<unsigned int>> splitchain;
			//
			vector<SplitChain> splitchains;
			SPLITChain(ExtendClusters, splitchains, Primary_chains[p].chains[h].link, smallOpts);
			int LSC = LargestSplitChain(splitchains);
			//
			// Apply SDP on all splitchains to get the final rough alignment path;
			// store the result in GenomePairs tupChain; 
			// We need vector<Cluster> tupClusters for tackling anchors of different strands
			// NOTICE: Insert 4 points for anchors in the overlapping regions between Clusters;
			//
			for (int st = 0; st < splitchains.size(); st++) {
				//
				// Apply SparseDP on extended anchors on every splitchain;
				// INPUT: vector<unsigned int> splitchain, vector<Cluster> ExtendClusters; OUTPUT: FinalChain finalchain;
				//
				FinalChain finalchain(&ExtendClusters);
				SparseDP(splitchains[st], ExtendClusters, finalchain, smallOpts, LookUpTable, read);
				//
				// RemoveSpuriousAnchors and RemovePairedIndels; 
				//
				int orig = finalchain.size();
				//RemoveSpuriousAnchors(finalchain, smallOpts);
				RemovePairedIndels(finalchain, smallOpts);
				if (finalchain.size() == 0) continue; // cannot be mapped to the genome!
				if (opts.dotPlot) {
					ofstream clust("SparseDP.tab", std::ofstream::app);
					for (int ep = 0; ep < finalchain.chain.size(); ep++) {
						if (finalchain.strand(ep) == 0) {
							clust << finalchain[ep].first.pos << "\t"
								  << finalchain[ep].second.pos << "\t"
								  << finalchain[ep].first.pos + finalchain.length(ep) << "\t"
								  << finalchain[ep].second.pos + finalchain.length(ep) << "\t"
								  << h << "\t"
								  << finalchain.ClusterNum(ep) << "\t"
								  << finalchain.strand(ep) << endl;
						}
						else {
							clust << finalchain[ep].first.pos << "\t"
								  << finalchain[ep].second.pos + finalchain.length(ep) << "\t"
								  << finalchain[ep].first.pos + finalchain.length(ep) << "\t"
								  << finalchain[ep].second.pos << "\t"
								  << h << "\t"
								  << finalchain.ClusterNum(ep) << "\t"
								  << finalchain.strand(ep) << endl;					
						}
					}
					clust.close();
				}	
				//
				// If there are inversions in the path, then seperate finalchain into several parts. 
				// In each part anchors are in the same direction, which is for better manipulation;
				// INPUT: FinalChain finalchain, OUTPUT: vector<int> finalchainSeperateChain;
				//
				vector<vector<int>> finalSeperateChain;
				SeperateChainByStrand(finalchain, finalSeperateChain, ExtendClusters); 
				int LFC = LargestFinalSeperateChain(finalSeperateChain);
				//
				// Refine and store the alignment; NOTICE: when filling in the space between two adjacent anchors, 
				// the process works in forward direction, so we need to flip the small matches
				// found in the spaces before insert them into the alignment if necessary;
				// INPUT: vector<int> finalchainSeperateStrand; OUTPUT: Alignment *alignment;
				//
				for (int fsc = 0; fsc < finalSeperateChain.size(); fsc++) {

					assert(finalSeperateChain[fsc].size() == 2);

					int start = finalSeperateChain[fsc][0];
					int end = finalSeperateChain[fsc][1];
					assert(start < end); 
					int str = finalchain.strand(start);
					int cln = finalchain.ClusterNum(start);
					int chromIndex = ExtendClusters[cln].chromIndex;	
					Alignment *alignment = new Alignment(strands[str], read.seq, read.length, read.name, str, genome.seqs[chromIndex],  
																							 genome.lengths[chromIndex], genome.header.names[chromIndex], chromIndex, 0); 
					alignments.back().SegAlignment.push_back(alignment);
					vector<int> scoreMat;
					vector<Arrow> pathMat;
					//
					// Set the secondary or supplymentary flag for the alignment; 
					// h == 0 ---> primary alignment, otherwise secondary alignment;
					// st > 0 ----> supplymentary alignment; fsc > 0 ----> supplymeantary alignment;
					//
					if (h > 0) alignment->ISsecondary = 1;
					if (splitchains.size() > 1 or finalSeperateChain.size() > 1) alignment->split = 1;					
					if (st != LSC or fsc != LFC) alignment->Supplymentary = 1; //if (st > 0 or fsc > 0) alignment->Supplymentary = 1;	

					GenomePos genomeThres = 0;
					if (str == 0) {
						int fl = end - 1;
						while (fl > start) {
							int cur = fl;
							int next = fl - 1;
							if (!PassgenomeThres(cur, genomeThres, finalchain)) {
								RefinedAlignmentbtwnAnchors(cur, next, str, chromIndex, finalchain, alignment, read, genome,
															 strands, scoreMat, pathMat, tinyOpts, genomeThres);
							}
							fl--;
						}
						if (!PassgenomeThres(start, genomeThres, finalchain)) {
							alignment->blocks.push_back(Block(finalchain[start].first.pos, finalchain[start].second.pos, 
																	finalchain.length(start)));
						}
					}
					else {
						int fl = start; 
						while (fl < end - 1) {
							int cur = fl;
							int next = fl + 1;
							if (!PassgenomeThres(cur, genomeThres, finalchain)) {
								RefinedAlignmentbtwnAnchors(cur, next, str, chromIndex, finalchain, alignment, read, genome,
															strands, scoreMat, pathMat, tinyOpts, genomeThres);
							}
							fl++;
						}	
						if (!PassgenomeThres(end - 1, genomeThres, finalchain)) {
							alignment->blocks.push_back(Block(read.length-finalchain[end - 1].first.pos-finalchain.length(end - 1), 
																finalchain[end - 1].second.pos, finalchain.length(end - 1)));
						}
					}
				
					for (int bb = 1; bb < alignment->blocks.size(); bb++) {
						assert(alignment->blocks[bb-1].qPos + alignment->blocks[bb-1].length <= alignment->blocks[bb].qPos);
						assert(alignment->blocks[bb-1].tPos + alignment->blocks[bb-1].length <= alignment->blocks[bb].tPos);	
					}
					//
					// Set some parameters in Alignment alignment
					//
					int nm=0;
					for(int b=0; b < alignment->blocks.size(); b++) {
						nm+= alignment->blocks[b].length;
					}
					alignment->read = strands[str];
					alignment->strand = str;
					alignment->nblocks = end - 1 - start;
					alignment->CalculateStatistics(smallOpts,svsigstrm, LookUpTable);
				}
			}
			alignments.back().SetFromSegAlignment();
		}
		alignmentsOrder.Update(&alignments);
	}	
	//AlignmentsOrder alignmentsOrder(&alignments);
	SimpleMapQV(alignmentsOrder, read);

	timing.Tick("Local-SDP");
	if (opts.dotPlot) {
			ofstream baseDots("alignment.dots");
			for (int a=0; a < (int) alignmentsOrder.size(); a++){
				for (int s = 0; s < alignmentsOrder[a].SegAlignment.size(); s++) {

					for (int c = 0; c < alignmentsOrder[a].SegAlignment[s]->blocks.size(); c++) {
						if (alignmentsOrder[a].SegAlignment[s]->strand == 0) {
							baseDots << alignmentsOrder[a].SegAlignment[s]->blocks[c].qPos << "\t" 
									 << alignmentsOrder[a].SegAlignment[s]->blocks[c].tPos << "\t" 
									 << alignmentsOrder[a].SegAlignment[s]->blocks[c].qPos + alignmentsOrder[a].SegAlignment[s]->blocks[c].length << "\t" 
									 << alignmentsOrder[a].SegAlignment[s]->blocks[c].tPos + alignmentsOrder[a].SegAlignment[s]->blocks[c].length << "\t"
									 << a << "\t"
									 << s << "\t"
									 << alignmentsOrder[a].SegAlignment[s]->strand << endl;							
						} 
						else {
							baseDots << read.length - alignmentsOrder[a].SegAlignment[s]->blocks[c].qPos - alignmentsOrder[a].SegAlignment[s]->blocks[c].length << "\t" 
									 << alignmentsOrder[a].SegAlignment[s]->blocks[c].tPos + alignmentsOrder[a].SegAlignment[s]->blocks[c].length << "\t" 
									 << read.length - alignmentsOrder[a].SegAlignment[s]->blocks[c].qPos << "\t" 
									 << alignmentsOrder[a].SegAlignment[s]->blocks[c].tPos << "\t"
									 << a << "\t"
									 << s << "\t"
									 << alignmentsOrder[a].SegAlignment[s]->strand << endl;
						}
					}		
				}
			}
			
			baseDots.close();
	}
	if (alignmentsOrder.size() > 0 and alignmentsOrder[0].SegAlignment.size() > 0) {
		for (int a=0; a < (int) alignmentsOrder.size(); a++){
			for (int s = 0; s < alignmentsOrder[a].SegAlignment.size(); s++) {
				alignmentsOrder[a].SegAlignment[s]->order=s;
				
				if (opts.printFormat == "b") {
					alignmentsOrder[a].SegAlignment[s]->PrintBed(*output);
				}
				else if (opts.printFormat == "s") {
					alignmentsOrder[a].SegAlignment[s]->PrintSAM(*output, opts);
				}
				else if (opts.printFormat == "a") {
					alignmentsOrder[a].SegAlignment[s]->PrintPairwise(*output);
				}
				else if (opts.printFormat == "p" or opts.printFormat == "pc") {
					alignmentsOrder[a].SegAlignment[s]->PrintPAF(*output, opts.printFormat=="pc");
				}
			}
		}
	}
	else {
		if (opts.printFormat == "s") {
			Alignment unaligned;
			unaligned.read=read.seq;
			unaligned.readLen=read.length;
			unaligned.PrintSAM(*output,opts);
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
	if (alignments.size() > 0) {
		return 1;
	}
	else {
		return 0;
	}
	return 0;
}

#endif
