#ifndef MAP_READ_H_
#define MAP_READ_H_
#include <math.h>
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
#include "SparseDP_Forward.h"
#include "Merge.h"
#include "MergeAnchors.h"
#include "Chain.h"
#include "overload.h"
#include "LinearExtend.h"
#include "SplitClusters.h"
#include "Timing.h"

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

// Contain changes from master branch 
// void RemoveOverlappingClusters(vector<Cluster> &clusters, vector<int> &clusterOrder, Options &opts) {
// 	int a=0;
// 	int ovp=a;

// 	if (clusters.size() == 0) {
// 		return;
// 	}
// 	vector<long> forDiagonals, revDiagonals;
// 	vector<long> *diagPtr;
// 	int nForCandidates=0, nRevCandidates=0;
// 	int maxCand=opts.maxCandidates;
// 	vector<bool> keep(clusters.size(), true);
// 	std::map<long, vector<int> > diagToCluster;
// 	long targetDiag=0;
// 	for (a=0; a < clusters.size(); a++) { 
		
// 		int orderIndex=clusterOrder[a];
// 		clusters[orderIndex].rank=a;
// 		float num=1.0;
// 		float denom=1.0;
// 		long diag=(long)clusters[orderIndex].tStart - (long)clusters[orderIndex].qStart;
// 		bool foundDiag=false;
// 		long clusterDiag, clusterEndDiag;
// 		//		cerr << "processing cluster on query " << clusters[orderIndex].qStart << "\t" << clusters[orderIndex].qEnd << "\t" << diag << "\t" << orderIndex << "\t" << clusters[orderIndex].tStart << "\t" << clusters[orderIndex].tEnd << endl;
// 		if (clusters[orderIndex].strand == 0) {
// 			diagPtr = &forDiagonals;
// 		}
// 		else {
// 			diagPtr = &revDiagonals;
// 		}
// 		bool encompassed=false;		
// 		bool onDiag=false;
// 		bool nearPoint=false;
// 		long curClusterDiag=0;
// 		long diagDist=0;
// 		long targetClusterDist=0;
// 		long targetDiagDist=1000;
// 		for (int d=0; d < diagPtr->size() and encompassed == false; d++) {							
// 			curClusterDiag=(*diagPtr)[d];
// 			assert (diagToCluster.find(curClusterDiag) != diagToCluster.end());
// 			assert (diagToCluster[curClusterDiag].size() > 0);
// 			for (int di = 0; di < diagToCluster[curClusterDiag].size(); di++) {
// 				int c=diagToCluster[curClusterDiag][di];
// 				clusterDiag=(long)clusters[c].tStart - (long) clusters[c].qStart;
// 				clusterEndDiag=(long)clusters[c].tEnd - (long) clusters[c].qEnd;
// 				long fey=(long)clusters[c].tStart - (long)clusters[orderIndex].tEnd;
// 				long fex=(long)clusters[c].qStart - (long)clusters[orderIndex].qEnd;
// 				long efy=(long)clusters[orderIndex].tStart - (long)clusters[c].tEnd;
// 				long efx =(long)clusters[orderIndex].tStart - (long)clusters[c].tEnd;
// 				long fe=(long) sqrt(fex*fex+fey*fey);
// 				long ef=(long) sqrt(efx*efx+efy*efy);
// 				diagDist=min(fe,ef);
					
// 				if (clusters[c].EncompassesInRectangle(clusters[orderIndex],0.5)) {
// 					//					cout << "cluster " << c << " encompasses " << orderIndex << endl;
// 					encompassed=true;
// 					break;
// 				}
// 				else {
// 					//					cout << "cluster " << c << " does not encompass " << orderIndex << "\t" << clusters[c].tEnd-clusters[c].tStart << "\t" << clusters[orderIndex].tEnd - clusters[orderIndex].tStart << endl;
// 				}					
// 				if ((abs(clusterDiag - diag) < 1000) or (abs(clusterEndDiag - diag) < 1000)) {
// 					foundDiag=true;					
// 					onDiag = true;
// 					targetDiag=curClusterDiag;
// 					targetDiagDist = min(abs(clusterDiag - diag), abs(clusterEndDiag - diag));
// 					//break;
// 				}
// 				if (abs(diagDist) < 1000) {
// 					nearPoint=true;
// 					targetClusterDist=diagDist;
// 					targetDiag=curClusterDiag;
// 					//break;
// 				}
// 			}
// 		}
// 		/*

// 			if (encompassed == false and (onDiag==true or nearPoint==true)) {
// 				foundDiag=true;
// 				break;
// 			}
// 			else {
// 				if (encompassed) {
// 					break;
// 				}
// 			}
// 			}*/
// 		//
// 		// Add hard-coded check for long matches for whole-contig alignments.
// 		if (foundDiag == false and 
// 				( ( diagPtr->size() < maxCand  and encompassed == false) or (abs(clusters[orderIndex].tEnd - clusters[orderIndex].tStart) > 1000) ) ) {
// 			(*diagPtr).push_back(diag);
// 		 //cerr << "Creating diagonal " << diag << "\t" << clusters[orderIndex].matches.size() << "\t" << clusters[orderIndex].tEnd - clusters[orderIndex].tStart << endl;
// 			diagToCluster[diag].push_back(orderIndex);
// 			foundDiag=true;
// 		}
// 		else if (foundDiag == true and (encompassed == false or (abs(clusters[orderIndex].tEnd - clusters[orderIndex].tStart) > 2000) )) { // change from master
// 		//else if (foundDiag == true and encompassed == false) {
// 			/*			cerr << "Keeping match " << clusters[orderIndex].matches.size() << "\t" << orderIndex 
// 							<< "\ton diag " << diag << "\t" << diagDist << "\t" << (int) nearPoint << "\t" << (int) encompassed << "\t" << targetDiagDist << "\t" << targetClusterDist << "\t" << clusters[orderIndex].qStart << "\t" << clusters[orderIndex].tStart << endl;*/
// 			assert(targetDiag != 0);
// 			diagToCluster[targetDiag].push_back(orderIndex);
// 		}			
// 		else {
// 			//			cerr << "Discarding cluster of size " << clusters[orderIndex].matches.size() << " on diag " << diag << endl;
// 			clusters[orderIndex].matches.resize(0);
// 		}
// 	}
// 	int c=0;
// 	for (int i=0; i < clusters.size(); i++) {
// 		if (clusters[i].matches.size() > 0) {
// 			clusters[c] = clusters[i];
// 			c++;
// 		}
// 	}
// 	clusters.resize(c);
// }

void RemoveOverlappingClusters(vector<Cluster> &clusters, vector<int> &clusterOrder, Options &opts) {
	int a=0;
	int ovp=a;

	if (clusters.size() == 0) {
		return;
	}
	vector<long> forDiagonals, revDiagonals;
	vector<long> *diagPtr;
	int nForCandidates=0, nRevCandidates=0;
	int maxCand=opts.maxCandidates;
	vector<bool> keep(clusters.size(), true);
	std::map<long, vector<int> > diagToCluster;
	long targetDiag=0;
	for (a=0; a < clusters.size(); a++) { 
		
		int orderIndex=clusterOrder[a];
		clusters[orderIndex].rank=a;
		float num=1.0;
		float denom=1.0;
		long diag=(long)clusters[orderIndex].tStart - (long)clusters[orderIndex].qStart;
		bool foundDiag=false;
		long clusterDiag, clusterEndDiag;
		//		cerr << "processing cluster on query " << clusters[orderIndex].qStart << "\t" << clusters[orderIndex].qEnd << "\t" << diag << "\t" << orderIndex << "\t" << clusters[orderIndex].tStart << "\t" << clusters[orderIndex].tEnd << endl;
		if (clusters[orderIndex].strand == 0) {
			diagPtr = &forDiagonals;
		}
		else {
			diagPtr = &revDiagonals;
		}
		bool encompassed=false;		
		bool onDiag=false;
		bool nearPoint=false;
		long curClusterDiag=0;
		long diagDist=0;
		long targetClusterDist=0;
		long targetDiagDist=1000;
		for (int d=0; d < diagPtr->size() and encompassed == false; d++) {							
			curClusterDiag=(*diagPtr)[d];
			assert (diagToCluster.find(curClusterDiag) != diagToCluster.end());
			assert (diagToCluster[curClusterDiag].size() > 0);
			for (int di = 0; di < diagToCluster[curClusterDiag].size(); di++) {
				int c=diagToCluster[curClusterDiag][di];
				clusterDiag=(long)clusters[c].tStart - (long) clusters[c].qStart;
				clusterEndDiag=(long)clusters[c].tEnd - (long) clusters[c].qEnd;
				long fey=(long)clusters[c].tStart - (long)clusters[orderIndex].tEnd;
				long fex=(long)clusters[c].qStart - (long)clusters[orderIndex].qEnd;
				long efy=(long)clusters[orderIndex].tStart - (long)clusters[c].tEnd;
				long efx =(long)clusters[orderIndex].tStart - (long)clusters[c].tEnd;
				long fe=(long) sqrt(fex*fex+fey*fey);
				long ef=(long) sqrt(efx*efx+efy*efy);
				diagDist=min(fe,ef);
					
				if (clusters[c].EncompassesInRectangle(clusters[orderIndex],0.5)) {
					//					cout << "cluster " << c << " encompasses " << orderIndex << endl;
					encompassed=true;
					break;
				}
				else {
					//					cout << "cluster " << c << " does not encompass " << orderIndex << "\t" << clusters[c].tEnd-clusters[c].tStart << "\t" << clusters[orderIndex].tEnd - clusters[orderIndex].tStart << endl;
				}					
				if ((abs(clusterDiag - diag) < 1000) or (abs(clusterEndDiag - diag) < 1000)) {
					foundDiag=true;					
					onDiag = true;
					targetDiag=curClusterDiag;
					targetDiagDist = min(abs(clusterDiag - diag), abs(clusterEndDiag - diag));
					//break;
				}
				if (abs(diagDist) < 1000) {
					nearPoint=true;
					targetClusterDist=diagDist;
					targetDiag=curClusterDiag;
					//break;
				}
			}

			if (encompassed == false and (onDiag==true or nearPoint==true)) {
				foundDiag=true;
				break;
			}
			else {
				if (encompassed) {
					break;
				}
			}
		}
		if (foundDiag == false and diagPtr->size() < maxCand and encompassed == false) {
			(*diagPtr).push_back(diag);
			//cerr << "Creating diagonal " << diag << "\t" << clusters[orderIndex].matches.size() << "\t" << clusters[orderIndex].tEnd - clusters[orderIndex].tStart << endl;
			diagToCluster[diag].push_back(orderIndex);
			foundDiag=true;
		}
		else if (foundDiag == true and encompassed == false) {
			/*			cerr << "Keeping match " << clusters[orderIndex].matches.size() << "\t" << orderIndex 
							<< "\ton diag " << diag << "\t" << diagDist << "\t" << (int) nearPoint << "\t" << (int) encompassed << "\t" << targetDiagDist << "\t" << targetClusterDist << "\t" << clusters[orderIndex].qStart << "\t" << clusters[orderIndex].tStart << endl;*/
			assert(targetDiag != 0);
			diagToCluster[targetDiag].push_back(orderIndex);
		}			
		else {
			//			cerr << "Discarding cluster of size " << clusters[orderIndex].matches.size() << " on diag " << diag << endl;
			clusters[orderIndex].matches.resize(0);
		}
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
			// if (alignmentsOrder[first].NumOfAnchors/(float)read.length < 0.00001) { // 0.0005
			// 	alignmentsOrder[first].mapqv = 2;
			// }
			// else {alignmentsOrder[first].mapqv = 60;}
			// alignmentsOrder[first].SetMapqv();
		}
		else {
			// assign mapqv to primary alignment;
			int nmmdiff = alignmentsOrder[first].nmm - alignmentsOrder[first+1].nmm;
			int nsmallgap = (alignmentsOrder[first].ndel + alignmentsOrder[first].nins) - 
								(alignmentsOrder[first+1].ndel + alignmentsOrder[first+1].nins);
			float alnvaluediff = alignmentsOrder[first].totalValue - alignmentsOrder[first+1].totalValue;

			//cerr << "nmmdiff: " << nmmdiff << " nsmallgap: " << nsmallgap << " alnvaluediff: " << alnvaluediff << endl;
			//cerr << "alignmentsOrder[first].NumOfAnchors/(float)read.length: " << alignmentsOrder[first].NumOfAnchors/(float)read.length << endl;
			float denom_1 = 1, denom_2 = 1;
			if (nmmdiff == 0 and nsmallgap == 0 ) {
				alignmentsOrder[first].mapqv = 2;
				alignmentsOrder[first+1].mapqv = 1;
			}
			// else if (alignmentsOrder[first].NumOfAnchors/(float)read.length < 0.00001) {
			// 	alignmentsOrder[first].mapqv = 2;
			// 	alignmentsOrder[first+1].mapqv = 1;
			// }
			else if (alnvaluediff <= 20 and abs(nmmdiff) <= 20 and abs(nsmallgap) <= 20) { //10 & 15 & 20
				alignmentsOrder[first].mapqv = 2;
				alignmentsOrder[first+1].mapqv = 1;
			}
			else if (nmmdiff <= 0 and nsmallgap <= 0) { // nmmdiff=0 and nsmallgap=0 ==> mapqv=52
				if (alignmentsOrder[first].totalValue < alignmentsOrder[first+1].totalValue + 300) { 
					if (nmmdiff < -20) denom_1 = 1;
					else if (nmmdiff >= -20 and nmmdiff < -5) denom_1 = pow(0.9,20-abs(nmmdiff));  
					else denom_1 = pow(0.7,20-abs(nmmdiff));

					if (nsmallgap < -20) denom_2 = 1;
					else if (nsmallgap >= -20 and nsmallgap < -5) denom_2 = pow(0.9,20-abs(nsmallgap));
					else  denom_2 = pow(0.7,20-abs(nsmallgap)); 
					alignmentsOrder[first].mapqv = min(60, (int) (60 + 5*log10(denom_1) + 5*log10(denom_2)));
				}
				else { alignmentsOrder[first].mapqv = 60;}			
			}
			else if (nmmdiff <= 0 and nsmallgap >= 0) { // when nsmallgap=30 and nmmdiff=0; mapqv=0
				if (alignmentsOrder[first].totalValue < alignmentsOrder[first+1].totalValue + 300) { 
					if (nmmdiff < -20) denom_1 = 1;
					else if (nmmdiff >= -20 and nmmdiff < -5) denom_1 = pow(0.9,20-abs(nmmdiff));  
					else denom_1 = pow(0.7,20-abs(nmmdiff));

					denom_2 = pow(0.1, abs(nsmallgap));  
					int cpr=(int) (60 + 5*log10(denom_1) + 20*log10(denom_2)); //60 + 10*log10(denom) 
					if (cpr>=0) alignmentsOrder[first].mapqv = min(60, cpr);
					else alignmentsOrder[first].mapqv = 2;					
				}
				else { alignmentsOrder[first].mapqv = 60;}					
			}
			else if (nmmdiff >= 0 and nsmallgap <= 0) { // when nsmallgap=0 and nmmdiff=30, mapqv=0
				if (alignmentsOrder[first].totalValue < alignmentsOrder[first+1].totalValue + 300) { 
					denom_1 = pow(0.1, abs(nmmdiff));  

					if (nsmallgap < -20) denom_2 = 1;
					else if (nsmallgap >= -20 and nsmallgap < -5) denom_2 = pow(0.9,20-abs(nsmallgap));
					else  denom_2 = pow(0.7,20-abs(nsmallgap)); 
					int cpr=(int) (60 + 20*log10(denom_1) + 5*log10(denom_2)); //60 + 10*log10(denom) 
					if (cpr>=0) alignmentsOrder[first].mapqv = min(60, cpr);
					else alignmentsOrder[first].mapqv = 2;						
				}
				else { alignmentsOrder[first].mapqv = 60;}	
			}
			else { 
				if (alignmentsOrder[first].totalValue < alignmentsOrder[first+1].totalValue + 300) { 
					denom_1 = pow(0.1, abs(nmmdiff));  
					denom_2 = pow(0.1, abs(nsmallgap));  
					int cpr=(int) (60 + 20*log10(denom_1) + 20*log10(denom_2));
					if (cpr>=0) alignmentsOrder[first].mapqv = min(60, cpr);
					else alignmentsOrder[first].mapqv = 2;					
				}
				else {alignmentsOrder[first].mapqv = 60;}	
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
										vector<int> &scoreMat, vector<Arrow> &pathMat, Alignment &aln, Options &options, 
										AffineAlignBuffers &buff) {
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
	int score = AffineOneGapAlign(readSeq, qLen, chromSeq, tLen, 
																options.localMatch, options.localMismatch, options.localIndel, min(drift*2+1,options.localBand), 
																aln, buff);
	/*
	cout << "aligned " << endl
			 << readSeq << endl
			 << chromSeq;
	aln.genomeLen = chromSeq.size();
	aln.read=(char*) readSeq.c_str();
	aln.genome=(char*) chromSeq.c_str();
	aln.CreateAlignmentStrings((char*) readSeq.c_str(), (char*)chromSeq.c_str(), aln.queryString, aln.alignString, aln.refString);
	aln.prepared=true;
	aln.PrintPairwise(cout);
	cout << endl;
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
int nSSE=0;
void RefineSubstrings(char *read, GenomePos readSubStart, GenomePos readSubEnd, char *genome, GenomePos genomeSubStart, GenomePos genomeSubEnd, 
											vector<int> &scoreMat, vector<Arrow> &pathMat, Alignment &aln, Options &opts,
											AffineAlignBuffers &buff) {
	aln.blocks.clear();
	AlignSubstrings(read, readSubStart, readSubEnd, genome, genomeSubStart, genomeSubEnd, scoreMat, pathMat, aln, opts, buff);
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
	Tuple Bi=1;
	Tuple rev_mask = Bi<<63;

	int nForward=0;
	for (int i=0; i < allMatches.size(); i++) {
		int readPos = allMatches[i].first.pos;
		uint64_t refPos = allMatches[i].second.pos;
		if ((allMatches[i].second.t & rev_mask) == 0) {
			// the last bit is 0, the match is in the forward strand
			nForward++;
		} 
		else {
			// the last bit is 1, the match is in the reverse strand
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


//
// This function removes spurious anchors after SDP;
//
void RemoveSpuriousAnchors(FinalChain & chain, Options & opts) {

 	if (chain.size() < 2) return;
	vector<bool> remove(chain.size(), false); // If remove[i] == true, then remove chain[i]
	vector<int> SV;
	vector<int> SVpos;

	//
	// Store SVs in vector SV; Store the anchor just after the SV[c] in SVpos[c];
	//
	for (int c = 1; c < chain.size(); c++) {
		if (chain.strand(c) == chain.strand(c-1)) {
			if (chain.strand(c) == 0) {
				int Gap = (int)(((long int)chain[c].second.pos - (long int)chain[c].first.pos) - 
							    ((long int)chain[c-1].second.pos - (long int)chain[c-1].first.pos));

				if (abs(Gap) >= 2000) { 
					SV.push_back(Gap);	
					SVpos.push_back(c);
				}
			}
			else {
				int Gap = (int)((long int)(chain[c].first.pos + chain[c].second.pos) - 
								(long int)(chain[c-1].first.pos + chain[c-1].second.pos));
				if (abs(Gap) >= 2000) {
					SV.push_back(Gap);	
					SVpos.push_back(c);
				}				
			}
		}
		else {
			SVpos.push_back(c);
			SV.push_back(0);
		}
	}

	for (int c = 1; c < SV.size(); c++) {
		if (SV[c] != 0 and SV[c-1] != 0) {
			bool _check = 0;
			if (SVpos[c] - SVpos[c-1] <= 10) {
				for (int bkc = SVpos[c-1]; bkc < SVpos[c]; bkc++) {
					if (chain.length(bkc) >= 50) {
						_check = 1;
						break;
					}
				}
				if (_check == 0) {
					for (int i = SVpos[c-1]; i < SVpos[c]; i++) {
						if (chain.length(i) < 50) remove[i] = true; // 200
					}				
				}	
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
// This function removes paired indels in the finalchain after SDP;
//
void RemovePairedIndels (FinalChain &chain) {

 	if (chain.size() < 2) return;
	vector<bool> remove(chain.size(), false); // If remove[i] == true, then remove chain[i]
	vector<int> SV;
	vector<int> SVpos;
	vector<long> SVgenome;
	//int s = 0, e = 0;

	//
	// Store SVs in vector SV; Store the anchor just after the SV[c] in SVpos[c];
	//
	for (int c = 1; c < chain.size(); c++) {
		if (chain.strand(c) == chain.strand(c-1)) {
			if (chain.strand(c) == 0) {
				int Gap = (int)(((long int)chain[c].second.pos - (long int)chain[c].first.pos) - 
							    ((long int)chain[c-1].second.pos - (long int)chain[c-1].first.pos));

				if (abs(Gap) > 30) { //30
					SV.push_back(Gap);	
					SVgenome.push_back(chain[c].second.pos);
					SVpos.push_back(c);
				}
			}
			else {
				int Gap = (int)((long int)(chain[c].first.pos + chain[c].second.pos) - 
								(long int)(chain[c-1].first.pos + chain[c-1].second.pos));
				if (abs(Gap) > 30) { //30
					SV.push_back(Gap);	
					SVgenome.push_back(chain[c].second.pos);
					SVpos.push_back(c);
				}				
			}
		}
		else {
			SVgenome.push_back(chain[c].second.pos);
			SVpos.push_back(c);
			SV.push_back(0);
		}
	}

	for (int c = 1; c < SV.size(); c++) {
		//
		// If two adjacent SVs have different types and similar lengths, then delete anchors in between those two SVs.
		// The last condition is to ensure both SV[c] and SV[c-1] are not zeros.
		//
		int blink = max(abs(SV[c]), abs(SV[c-1]));
		if (sign(SV[c]) != sign(SV[c-1]) and SV[c] != 0 and SV[c-1] != 0 and abs(SV[c] + SV[c-1]) < 600) { 
			if ((sign(SV[c]) == true and abs(SVgenome[c] - SVgenome[c-1]) < max(2*blink, 1000)) // SV[c] is ins
				or (sign(SV[c]) == false and abs(SVgenome[c] - SV[c] - SVgenome[c-1]) < max(2*blink, 1000))) { // SV[c] is del
				//
				// remove anchors from SVpos[c-1] to SV[c];
				//
				for (int i = SVpos[c-1]; i < SVpos[c]; i++) {
					if (chain.length(i) < 100) remove[i] = true; // 200
				}
			}			
		} 
		else if (sign(SV[c]) != sign(SV[c-1]) and SV[c] != 0 and SV[c-1] != 0 
					and ((sign(SV[c]) == true and abs(SVgenome[c] - SVgenome[c-1]) < 500) 
						or (sign(SV[c]) == false and abs(SVgenome[c] - SV[c] - SVgenome[c-1]) < 500))) {
				for (int i = SVpos[c-1]; i < SVpos[c]; i++) {
					if (chain.length(i) < 100) remove[i] = true; // 200
				}			
		}
		else if (sign(SV[c]) == sign(SV[c-1]) and SV[c] != 0 and SV[c-1] != 0) { // If two gaps of same typeare too close (<600bp)
			if ((sign(SV[c]) == true and abs(SVgenome[c] - SVgenome[c-1]) < max(2*blink, 1000)) // insertion
				or (sign(SV[c]) == false and abs(SVgenome[c] - SV[c] - SVgenome[c-1]) < max(2*blink, 1000))) { //deletion
				//
				// remove anchors from SVpos[c-1] to SV[c];
				//
				for (int i = SVpos[c-1]; i < SVpos[c]; i++) {
					if (chain.length(i) < 100) remove[i] = true; 
				}
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
// This function removes paired indels in the finalchain after SDP;
//
void RemovePairedIndels (GenomePairs &matches, vector<unsigned int> &chain, vector<int> &lengths) {

 	if (chain.size() < 2) return;
	vector<bool> remove(chain.size(), false); // If remove[i] == true, then remove chain[i]
	vector<int> SV;
	vector<int> SVpos;
	vector<int> SVgenome;
	//int s = 0, e = 0;

	//
	// Store SVs in vector SV; Store the anchor just after the SV[c] in SVpos[c];
	//
	for (int c = 1; c < chain.size(); c++) {
		int Gap = (int)(((long)matches[chain[c]].second.pos - (long)matches[chain[c]].first.pos) - 
					    ((long)matches[chain[c-1]].second.pos - (long)matches[chain[c-1]].first.pos));

		if (abs(Gap) > 30) {
			SV.push_back(Gap);	
			SVgenome.push_back(matches[chain[c]].second.pos);
			SVpos.push_back(c);
		}
	}

	for (int c = 1; c < SV.size(); c++) {
		//
		// If two adjacent SVs have different types and similar lengths, then delete anchors in between those two SVs.
		// The third condition is to ensure both SV[c] and SV[c-1] are not zeros.
		//
		int blink = max(abs(SV[c]), abs(SV[c-1])) ;
		if (sign(SV[c]) != sign(SV[c-1]) and abs(SV[c] + SV[c-1]) < 600 and abs(SV[c]) != 0 and SV[c-1] != 0) { 
			if ((sign(SV[c]) == true and abs(SVgenome[c] - SVgenome[c-1]) < max(2*blink, 1000)) // SV[c] is ins
				or (sign(SV[c]) == false and abs(SVgenome[c] - SV[c] - SVgenome[c-1]) < max(2*blink, 1000))) { // SV[c] is del
				//
				// remove anchors from SVpos[c-1] to SV[c];
				//
				for (int i = SVpos[c-1]; i < SVpos[c]; i++) {
					if (lengths[chain[i]]<100) remove[i] = true;
				}
			}	
		}
		else if (sign(SV[c]) != sign(SV[c-1]) and SV[c] != 0 and SV[c-1] != 0 
					and ((sign(SV[c]) == true and abs(SVgenome[c] - SVgenome[c-1]) < 500) 
						or (sign(SV[c]) == false and abs(SVgenome[c] - SV[c] - SVgenome[c-1]) < 500))) {
				for (int i = SVpos[c-1]; i < SVpos[c]; i++) {
					if (lengths[chain[i]]<100) remove[i] = true; // 200
				}			
		}
		else if (sign(SV[c]) == sign(SV[c-1]) and SV[c] != 0 and SV[c-1] != 0) {
			if ((sign(SV[c]) == true and abs(SVgenome[c] - SVgenome[c-1]) < max(2*blink, 1000)) // insertion
				or (sign(SV[c]) == false and abs(SVgenome[c] - SV[c] - SVgenome[c-1]) < max(2*blink, 1000))) { //deletion
				//
				// remove anchors from SVpos[c-1] to SV[c];
				//
				for (int i = SVpos[c-1]; i < SVpos[c]; i++) {
					if (lengths[chain[i]]<100) remove[i] = true; 
				}
			}			
		}
	}

	int m = 0;
	for (int i = 0; i < chain.size(); i++) {
		if (remove[i] == false) {
			chain[m] = chain[i];
			m++;
		}
	}
	chain.resize(m);
}

//
// This function switches index in splitclusters back 
//
int
switchindex (vector<Cluster> & splitclusters, vector<Primary_chain> & Primary_chains, vector<Cluster> & clusters, Genome &genome) {
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
				Primary_chains[p].chains[h].ch[c] = splitclusters[Primary_chains[p].chains[h].ch[c]].coarse;
			}
		}
	}
	//
	// Change "vector<bool> link" accordingly
	//
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
	//
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
	  		vector<unsigned int>::iterator itp;
	  		itp = unique(Primary_chains[p].chains[h].ch.begin(), Primary_chains[p].chains[h].ch.end()); 
	  		int pdist = distance(Primary_chains[p].chains[h].ch.begin(), itp);
	  		Primary_chains[p].chains[h].ch.resize(pdist);			
		}
	}
	//
	// Remove those Primary_chains which align across chromosome, but overlap a lot on read coordinates;
	//
	// for (int p = 0; p < Primary_chains.size(); p++) {
	// 	//vector<bool> remove(Primary_chains[p].chains.size(), 0);
	// 	for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
	// 		vector<bool> remove(Primary_chains[p].chains[h].ch.size(), 0);
 // 			int c = 1;
 // 			int c_cur, c_prev;
	// 		while (c < Primary_chains[p].chains[h].ch.size()) {
	// 			c_cur = Primary_chains[p].chains[h].ch[c];
	// 			c_prev = Primary_chains[p].chains[h].ch[c-1];
	// 			long clustlen = max(clusters[c_prev].qEnd - clusters[c_prev].qStart, clusters[c_cur].qEnd-clusters[c_cur].qStart);
	// 			long clustdist = 0, overlap_q = 0;
	// 			if (Primary_chains[p].chains[h].link[c-1] == 0) {
	// 				if (clusters[c_prev].tStart > clusters[c_cur].tEnd) clustdist = clusters[c_prev].tStart - clusters[c_cur].tEnd;
	// 			}
	// 			else {
	// 				if (clusters[c_cur].tStart > clusters[c_prev].tEnd) clustdist = clusters[c_cur].tStart - clusters[c_prev].tEnd;
	// 			}
	// 			if (clustdist > min(genome.lengths[clusters[c_cur].chromIndex], genome.lengths[clusters[c_prev].chromIndex])) {
	// 				if (clusters[c_prev].qStart < clusters[c_cur].qEnd) {
	// 					overlap_q = min(clusters[c_prev].qEnd, clusters[c_cur].qEnd)-max(clusters[c_prev].qStart, clusters[c_cur].qStart);
	// 					//cerr << "clustlen: " << clustlen << " overlap_q: " << overlap_q << endl; 
	// 					if (overlap_q > 0 and overlap_q / (float) clustlen > 0.5) {
	// 						remove[c] = 1;
	// 					}
	// 				}
	// 			}
	// 			c++;
	// 		}	
	// 		int cq = 0;
	// 		int lq = 0;
	// 		for (int cr = 0; cr < Primary_chains[p].chains[h].ch.size(); cr++) {	
	// 			if (remove[cr] == 0) {
	// 				Primary_chains[p].chains[h].ch[cq] = Primary_chains[p].chains[h].ch[cr];
	// 				cq++;
	// 				if (cr >= 1){
	// 					Primary_chains[p].chains[h].link[lq] = Primary_chains[p].chains[h].link[cr-1];
	// 					lq++;
	// 				}
	// 			}

	// 		}
	// 		Primary_chains[p].chains[h].ch.resize(cq);
	// 		Primary_chains[p].chains[h].link.resize(lq);
	// 	}
	// }	
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
		maxDN = (long long int) clusters[ph].matches[0].second.pos - (long long int) clusters[ph].matches[0].first.pos;
		minDN = maxDN;
		for (int db = 0; db < clusters[ph].matches.size(); db++) {
			maxDN = max(maxDN, (long long int)clusters[ph].matches[db].second.pos - (long long int)clusters[ph].matches[db].first.pos);
			minDN = min(minDN, (long long int)clusters[ph].matches[db].second.pos - (long long int)clusters[ph].matches[db].first.pos);
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

				CompareLists<LocalTuple, uint32_t>(readIndex->minimizers.begin()+qStartBoundary, 
									readIndex->minimizers.begin()+qEndBoundary, 
									glIndex.minimizers.begin()+ glIndex.tupleBoundaries[lsi], 
									glIndex.minimizers.begin()+ glIndex.tupleBoundaries[lsi+1], smallMatches, smallOpts, 0, 0, false);
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
	long long int minDiagNum, maxDiagNum; 
	long long int diag1, diag2;
	diag1 = 0;
	diag2 = (long long int) (te - (ts - lrts)) - (long long int) (qe - qs); // scale diag1 and diag2 to the local coordinates
	minDiagNum = min(diag1, diag2) - opts.refineSpaceDiag; 
	maxDiagNum = max(diag1, diag2) + opts.refineSpaceDiag; 

	vector<GenomeTuple> EndReadTup, EndGenomeTup;
	StoreMinimizers<GenomeTuple, Tuple>(genome.seqs[ChromIndex]+(ts-lrts), te-ts+lrlength, opts.globalK, opts.globalW, 
											EndGenomeTup, false);
	sort(EndGenomeTup.begin(), EndGenomeTup.end());
	StoreMinimizers<GenomeTuple, Tuple>(strands[st] + qs, qe - qs, opts.globalK, opts.globalW, EndReadTup, false);
	sort(EndReadTup.begin(), EndReadTup.end());
	CompareLists<GenomeTuple, Tuple>(EndReadTup.begin(), EndReadTup.end(), EndGenomeTup.begin(), EndGenomeTup.end(), EndPairs, 
										opts,  maxDiagNum, minDiagNum, false); // By passing maxDiagNum and minDiagNum, this function
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
				or ExtendClusters[cur].tEnd + opts.splitdist < ExtendClusters[prev].tStart) { //100,000
				
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
SeparateChainByStrand(FinalChain & finalchain, vector<vector<int>> & finalSeperateChain, const vector<Cluster> & ExtendClusters) {

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

void 
RefineByLinearAlignment(GenomePos &btc_curReadEnd, GenomePos &btc_curGenomeEnd, GenomePos &btc_nextReadStart, GenomePos &btc_nextGenomeStart, 
						int & str, int & chromIndex, Alignment * alignment, Read & read, Genome & genome, char *strands[2], 
						vector<int> & scoreMat, vector<Arrow> & pathMat, Options & opts, GenomePos & genomeThres,
						AffineAlignBuffers &buff) {
	//
	// find small matches between fragments in gapChain
	int m, rg, gg;
	SetMatchAndGaps(btc_curReadEnd, btc_nextReadStart, btc_curGenomeEnd, btc_nextGenomeStart, m, rg, gg);

	if (m > 0) {
		Alignment betweenAnchorAlignment;
		if (opts.refineLevel & REF_DP) {						
			RefineSubstrings(strands[str], btc_curReadEnd, btc_nextReadStart, genome.seqs[chromIndex], 
											 btc_curGenomeEnd, btc_nextGenomeStart, scoreMat, pathMat, betweenAnchorAlignment, opts, buff);
			int b;
			for (b = 1; b < betweenAnchorAlignment.blocks.size(); b++) {
				assert(betweenAnchorAlignment.blocks[b-1].qPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].qPos);
				assert(betweenAnchorAlignment.blocks[b-1].tPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].tPos);						
			}
			alignment->blocks.insert(alignment->blocks.end(), betweenAnchorAlignment.blocks.begin(), betweenAnchorAlignment.blocks.end());
			if (opts.dotPlot) {
				ofstream Lclust("LinearAlignment.tab", std::ofstream::app);
				for (b = 0; b < betweenAnchorAlignment.blocks.size(); b++) {
					Lclust << betweenAnchorAlignment.blocks[b].qPos << "\t"
						  << betweenAnchorAlignment.blocks[b].tPos << "\t"
						  << betweenAnchorAlignment.blocks[b].qPos + betweenAnchorAlignment.blocks[b].length << "\t"
						  << betweenAnchorAlignment.blocks[b].tPos + betweenAnchorAlignment.blocks[b].length << "\t"
						  << str << endl;
				}
				Lclust.close();
			}	
			//
			// Debug
			//
			// for (int bb = 1; bb < alignment->blocks.size(); bb++) {
			// 	assert(alignment->blocks[bb-1].qPos + alignment->blocks[bb-1].length <= alignment->blocks[bb].qPos);
			// 	assert(alignment->blocks[bb-1].tPos + alignment->blocks[bb-1].length <= alignment->blocks[bb].tPos);	
			// }
			betweenAnchorAlignment.blocks.clear();
		}
	}	
	genomeThres = 0;	
}


//
// This function creates the alignment for a "segment"(one big chunk) on the chain;
//
void
RefinedAlignmentbtwnAnchors(int & cur, int & next, int & str, int & chromIndex, FinalChain & finalchain, Alignment *alignment, 
							Read & read, Genome & genome, char *strands[2], vector<int> & scoreMat, vector<Arrow> & pathMat, 
							Options & tinyOpts, GenomePos & genomeThres, AffineAlignBuffers &buff, const vector<float> & LookUpTable) {

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
	//
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
	assert(curReadEnd <= nextReadStart);
	
	if (curGenomeEnd <= nextGenomeStart) {
		long read_dist = nextReadStart - curReadEnd;
		long genome_dist = nextGenomeStart - curGenomeEnd;
		//
		// The linear alignment is more suitable to find a big gap; 
		//
		//if (tinyOpts.RefineBySDP == true and ((abs(read_dist-genome_dist) <= 1000 and min(read_dist, genome_dist) >= 300) or 
		//												(min(read_dist, genome_dist) > 1000))) { 
		if (tinyOpts.RefineBySDP == true and min(read_dist, genome_dist) >= 300) {
			GenomePairs BtwnPairs;
			//cerr << "abs(read_dist-genome_dist): "  << abs(read_dist-genome_dist)<< endl;
			//cerr << "min(read_dist, genome_dist): " << min(read_dist, genome_dist) << endl;
			//cerr << "curReadEnd: " << curReadEnd << "  curGenomeEnd: " << curGenomeEnd << "  nextReadStart: " << nextReadStart << "  nextGenomeStart: " << nextGenomeStart << endl;
			if (abs(read_dist-genome_dist) < 50) tinyOpts.refineSpaceDiag = 30;
			else tinyOpts.refineSpaceDiag = 80;
			Options nanoOpts = tinyOpts;
			if (min(read_dist, genome_dist) >= 3000 and abs(read_dist-genome_dist) >= 3000) { // big tandem repeat
				nanoOpts.globalK += 3;
			}
			// tinyOpts.refineSpaceDiag = 100; //30
			// if (min(read_dist, genome_dist) >= 30000 and abs(read_dist-genome_dist) < 500) {
			// 	tinyOpts.refineSpaceDiag = 60;//7000
			// }
			RefineSpace(0, BtwnPairs, nanoOpts, genome, read, strands, chromIndex, nextReadStart, curReadEnd, nextGenomeStart, 
							curGenomeEnd, str);
			if (min(read_dist, genome_dist) > 30000 and (BtwnPairs.size() / (float) min(read_dist, genome_dist)) < 0.03) { // this case likely means there is a large SV events. 
				BtwnPairs.clear();
				nanoOpts.refineSpaceDiag = 8000;
				RefineSpace(0, BtwnPairs, nanoOpts, genome, read, strands, chromIndex, nextReadStart, curReadEnd, nextGenomeStart, 
							curGenomeEnd, str);		
			} 
			//cerr << "Refine tiny space done!" << endl;
			GenomePairs ExtendBtwnPairs;
			vector<int> ExtendBtwnPairsMatchesLength;
			vector<unsigned int> BtwnChain;
			LinearExtend(BtwnPairs, ExtendBtwnPairs, ExtendBtwnPairsMatchesLength, nanoOpts, genome, read, chromIndex);
			//
			// insert the next anchor;
			//
			ExtendBtwnPairs.push_back(GenomePair(GenomeTuple(0, nextReadStart), 
												 GenomeTuple(0, nextGenomeStart)));
			ExtendBtwnPairsMatchesLength.push_back(finalchain.length(next));			
			// 
			// insert the previous anchor;
			//
			ExtendBtwnPairs.push_back(GenomePair(GenomeTuple(0, curReadEnd-finalchain.length(cur)), 
												 GenomeTuple(0, curGenomeEnd-finalchain.length(cur))));
			ExtendBtwnPairsMatchesLength.push_back(finalchain.length(cur));
			TrimOverlappedAnchors(ExtendBtwnPairs, ExtendBtwnPairsMatchesLength);
			nanoOpts.coefficient=12;//9 this is calibrately set to 12

			SparseDP_ForwardOnly(ExtendBtwnPairs, ExtendBtwnPairsMatchesLength, BtwnChain, nanoOpts, LookUpTable, 1); //6
			RemovePairedIndels(ExtendBtwnPairs, BtwnChain, ExtendBtwnPairsMatchesLength);
			//
			// Use linear alignment to ligand the gap between local SDP chain;
			//
			GenomePos btc_curReadEnd, btc_nextReadStart, btc_curGenomeEnd, btc_nextGenomeStart;
			btc_curReadEnd = curReadEnd;
			btc_curGenomeEnd = curGenomeEnd;

			if (nanoOpts.dotPlot) {
				ofstream pSclust("BtwnPairs.tab", std::ofstream::app);
				for (int bp = BtwnPairs.size()-1; bp >= 0; bp--) {
					if (str == 0) {
						pSclust << BtwnPairs[bp].first.pos << "\t"
							  << BtwnPairs[bp].second.pos << "\t"
							  << BtwnPairs[bp].first.pos + nanoOpts.globalK << "\t"
							  << BtwnPairs[bp].second.pos + nanoOpts.globalK << "\t"
							  << str << endl;
					}
					else {
						pSclust << read.length - BtwnPairs[bp].first.pos - nanoOpts.globalK << "\t"
							  << BtwnPairs[bp].second.pos + nanoOpts.globalK << "\t"
							  << read.length - BtwnPairs[bp].first.pos << "\t"
							  << BtwnPairs[bp].second.pos << "\t"
							  << str << endl;					
					}					
				}
				pSclust.close();

				ofstream eSclust("ExtendBtwnPairs.tab", std::ofstream::app);
				for (int bp = ExtendBtwnPairs.size()-1; bp >= 0; bp--) {
					if (str == 0) {
						eSclust << ExtendBtwnPairs[bp].first.pos << "\t"
							  << ExtendBtwnPairs[bp].second.pos << "\t"
							  << ExtendBtwnPairs[bp].first.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
							  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
							  << str << endl;
					}
					else {
						eSclust << read.length - ExtendBtwnPairs[bp].first.pos - ExtendBtwnPairsMatchesLength[bp] << "\t"
							  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
							  << read.length - ExtendBtwnPairs[bp].first.pos << "\t"
							  << ExtendBtwnPairs[bp].second.pos << "\t"
							  << str << endl;					
					}					
				}
				eSclust.close();
			}	

			if (nanoOpts.dotPlot) {
				ofstream Sclust("SparseDP_Forward.tab", std::ofstream::app);
					for (int btc = BtwnChain.size()-1; btc >= 0; btc--) {
						if (str == 0) {
							Sclust << ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
								  << ExtendBtwnPairs[BtwnChain[btc]].second.pos << "\t"
								  << ExtendBtwnPairs[BtwnChain[btc]].first.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
								  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
								  << str << endl;
						}
						else {
							Sclust << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos - ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
								  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
								  << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
								  << ExtendBtwnPairs[BtwnChain[btc]].second.pos << "\t"
								  << str << endl;					
						}
					}
				Sclust.close();
			}
			int btc_end = BtwnChain.size()-1; int btc_start = 0;
			if (BtwnChain.back() == ExtendBtwnPairs.size()-1) btc_end = BtwnChain.size()-2;
			if (BtwnChain[0] == ExtendBtwnPairs.size()-2) btc_start = 1;
			
			for (int btc = btc_end; btc >= btc_start; btc--) {
				btc_nextGenomeStart = ExtendBtwnPairs[BtwnChain[btc]].second.pos;
				btc_nextReadStart = ExtendBtwnPairs[BtwnChain[btc]].first.pos;
				RefineByLinearAlignment(btc_curReadEnd, btc_curGenomeEnd, btc_nextReadStart, btc_nextGenomeStart, str, chromIndex, 
										alignment, read, genome, strands, scoreMat, pathMat, nanoOpts, genomeThres, buff);				

				alignment->blocks.push_back(Block(ExtendBtwnPairs[BtwnChain[btc]].first.pos, ExtendBtwnPairs[BtwnChain[btc]].second.pos, 
													ExtendBtwnPairsMatchesLength[BtwnChain[btc]]));

				// This is for debugging
				if (alignment->blocks.size() > 1) {
					int last=alignment->blocks.size();
					assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].qPos);
					assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].tPos);
				}	
				btc_curReadEnd = btc_nextReadStart + ExtendBtwnPairsMatchesLength[BtwnChain[btc]];
				btc_curGenomeEnd = btc_nextGenomeStart + ExtendBtwnPairsMatchesLength[BtwnChain[btc]];			
			}
			//
			// Add the linear alignment after the last anchor on BtwnChain;
			//
			btc_nextGenomeStart = nextGenomeStart;
			btc_nextReadStart = nextReadStart;			
			if (btc_nextGenomeStart > btc_curGenomeEnd and btc_nextReadStart > btc_curReadEnd) {
				RefineByLinearAlignment(btc_curReadEnd, btc_curGenomeEnd, btc_nextReadStart, btc_nextGenomeStart, str, chromIndex, 
										alignment, read, genome, strands, scoreMat, pathMat, nanoOpts, genomeThres, buff);					
			}
		}
		else {
			RefineByLinearAlignment(curReadEnd, curGenomeEnd, nextReadStart, nextGenomeStart, str, chromIndex, 
										alignment, read, genome, strands, scoreMat, pathMat, tinyOpts, genomeThres, buff);				
		}
	}
	else genomeThres = curGenomeEnd;
}

int 
PassgenomeThres(int cur, GenomePos & genomeThres, FinalChain & finalchain) {
	if (genomeThres != 0 and finalchain[cur].second.pos < genomeThres) return 1;
	else return 0;
}

int MapRead(const vector<float> & LookUpTable, Read &read, Genome &genome,
						vector<GenomeTuple> &genomemm, LocalIndex &glIndex, Options &opts, 
						ostream *output, ostream *svsigstrm, Timing &timing, pthread_mutex_t *semaphore=NULL) {
	string baseName = read.name;

	for (int i=0; i < baseName.size(); i++) {	
		if (baseName[i] == '/') baseName[i] = '_';	
		if (baseName[i] == '|') baseName[i] = '_';
	}


	vector<GenomeTuple> readmm; // readmm stores minimizers
	vector<pair<GenomeTuple, GenomeTuple> > allMatches, forMatches, revMatches, matches;
	timing.Start();
	//
	// Add pointers to seq that make code more readable.
	//
	char *readRC;
	CreateRC(read.seq, read.length, readRC);
	char *strands[2] = { read.seq, readRC };

	if (opts.storeAll) {
		Options allOpts = opts;
		allOpts.globalW=1;
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, allOpts.globalK, allOpts.globalW, readmm, false);			
	}
	else {
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, opts.globalK, opts.globalW, readmm, false);
	}
	sort(readmm.begin(), readmm.end()); //sort kmers in readmm(minimizers)

	timing.Tick("Store minimizers");
	//
	// Add matches between the read and the genome.
	//
	CompareLists<GenomeTuple, Tuple>(readmm, genomemm, allMatches, opts);

	// // check duplicate in allMatches
	// for (int am=0; am<allMatches.size(); am++) {
	// 	if (allMatches[am].first.t == 5413 and allMatches[am].first.pos == 47435 and allMatches[am].second.t == 5413 and allMatches[am].second.pos == 36356) {
	// 		cerr << "before sort am: " << am << endl;
	// 	}
	// }
	// Guess that 500 is where the setup/takedown overhead of an array index is equal to sort by value.
	DiagonalSort<GenomeTuple>(allMatches,500); // sort fragments in allMatches by forward diagonal, then by first.pos(read)

	timing.Tick("Sort minimizers");

	//TODO(Jinwen): delete this after debug
	if (opts.dotPlot) {
		ofstream clust("all-matches.dots");
		for (int m=0; m < allMatches.size(); m++) {
			clust << allMatches[m].first.pos << "\t" << allMatches[m].second.pos 
						<< "\t" << allMatches[m].first.pos+ opts.globalK << "\t" 
						<< allMatches[m].second.pos + opts.globalK << "\t0\t0"<<endl;
		}
		clust.close();
	}

	// check duplicate in allMatches
	// for (int am=0; am<allMatches.size(); am++) {
	// 	if (allMatches[am].first.t == 5413 and allMatches[am].first.pos == 47435 and allMatches[am].second.t == 5413 and allMatches[am].second.pos == 36356) {
	// 		cerr << "am: " << am << endl;
	// 	}
	// }

	SeparateMatchesByStrand(read, genome, opts.globalK, allMatches, forMatches, revMatches);
	CleanOffDiagonal(forMatches, opts);
	// if (opts.dotPlot and read.name == "chr1/tig00000001|arrow") {
	// 	ofstream clust("off-diag-all-matches.dots");
	// 	for (int m=0; m < allMatches.size(); m++) {
	// 		clust << allMatches[m].first.pos << "\t" << allMatches[m].second.pos << "\t" << allMatches[m].first.pos+ opts.globalK << "\t" 
	// 			<< allMatches[m].second.pos + opts.globalK << "\t0\t0"<<endl;
	// 	}
	// 	clust.close();
	// }

	//cerr << "opts.globalK " << opts.globalK << endl;
	vector<Cluster> clusters;
	vector<Cluster> roughClusters;
	vector<Cluster> split_roughClusters;
	int forwardStrand=0;
	// maxDiag must be large enough for the following function "StoreDiagonalClusters". 
	// Otherwise fragments on the same line (line with little curve) might end up in several clusters[i], instead of one clusters[i]

	// The strategy we are taking: 
	// first let maxDiag be a small number but not too small(like 500), which also alleviate the cases where anchors are on a curvy line.
	// Then break the clusters[i] into two if any two anchors are father than maxGap. 
	StoreDiagonalClusters(forMatches, roughClusters, opts, 0, forMatches.size(), forwardStrand); // rough == true means only storing "start and end" in every clusters[i]
	for (int c = 0; c < roughClusters.size(); c++) {
		CartesianSort(forMatches, roughClusters[c].start, roughClusters[c].end);
		SplitRoughClustersWithGaps(forMatches, roughClusters[c], split_roughClusters, opts, 0);
	}
	//cerr << "roughClusters.size(): " << roughClusters.size() << " split_roughClusters.size(): " << split_roughClusters.size()<< endl;
	if (opts.dotPlot) {
		ofstream clust("for-matches.dots");
		for (int m=0; m < forMatches.size(); m++) {
			clust << forMatches[m].first.pos << "\t" << forMatches[m].second.pos << "\t" << opts.globalK + forMatches[m].first.pos << "\t"
					<< forMatches[m].second.pos + opts.globalK << "\t" << m << "\t0"<<endl;
		}
		clust.close();
		ofstream rclust("roughClusters.dots");
		for (int m=0; m < roughClusters.size(); m++) {
			int rci = genome.header.Find(roughClusters[m].tStart);
			for (int c = roughClusters[m].start; c < roughClusters[m].end; ++c) {
				rclust << forMatches[c].first.pos << "\t" << forMatches[c].second.pos << "\t" << opts.globalK + forMatches[c].first.pos << "\t"
					<< forMatches[c].second.pos + opts.globalK << "\t" << m << "\t" << genome.header.names[rci]<< "\t0"<<endl;				
			}
		}
		rclust.close();
		ofstream wclust("split_roughClusters.dots");
		for (int m=0; m < split_roughClusters.size(); m++) {
			int rci = genome.header.Find(split_roughClusters[m].tStart);
			for (int c = split_roughClusters[m].start; c < split_roughClusters[m].end; ++c) {
				wclust << forMatches[c].first.pos << "\t" << forMatches[c].second.pos << "\t" << opts.globalK + forMatches[c].first.pos << "\t"
					<< forMatches[c].second.pos + opts.globalK << "\t" << m << "\t" << genome.header.names[rci]<< "\t0"<<endl;				
			}
		}
		wclust.close();

	}

	timing.Tick("Forward-diag-rough-clusters");

	AntiDiagonalSort<GenomeTuple>(revMatches, genome.GetSize(), 500);
	CleanOffDiagonal(revMatches, opts, 1);
	vector<Cluster> revroughClusters;
	vector<Cluster> split_revroughClusters;
	int reverseStrand=1;
	StoreDiagonalClusters(revMatches, revroughClusters, opts, 0, revMatches.size(), reverseStrand);
	for (int c = 0; c < revroughClusters.size(); c++) {
		CartesianSort(revMatches, revroughClusters[c].start, revroughClusters[c].end);
		SplitRoughClustersWithGaps(revMatches, revroughClusters[c], split_revroughClusters, opts, 1);
	}
	//cerr << "revroughClusters.size(): " << revroughClusters.size() << " split_revroughClusters.dots: " << split_revroughClusters.size()<< endl;


	if (opts.dotPlot) {
		ofstream rclust("rev-matches.dots");
		for (int m=0; m < revMatches.size(); m++) {			
			rclust << revMatches[m].first.pos << "\t" << revMatches[m].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[m].first.pos << "\t"
					 << revMatches[m].second.pos << "\t0"<<endl;
		}
		rclust.close();
		ofstream revsclust("revroughClusters.dots");
		for (int m=0; m < revroughClusters.size(); m++) {
			int rci = genome.header.Find(revroughClusters[m].tStart);
			for (int c = revroughClusters[m].start; c < revroughClusters[m].end; ++c) {
				revsclust << revMatches[c].first.pos << "\t" << revMatches[c].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[c].first.pos << "\t"
					 << revMatches[c].second.pos<< "\t" << m << "\t" << genome.header.names[rci]<<"\t0"<<endl;				
			}
		}
		revsclust.close();
		ofstream srevclust("split_revroughClusters.dots");
		for (int m=0; m < split_revroughClusters.size(); m++) {
			int rci = genome.header.Find(split_revroughClusters[m].tStart);
			for (int c = split_revroughClusters[m].start; c < split_revroughClusters[m].end; ++c) {
				srevclust << revMatches[c].first.pos << "\t" << revMatches[c].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[c].first.pos << "\t"
					 << revMatches[c].second.pos<< "\t" << m << "\t" << genome.header.names[rci]<<"\t0"<<endl;				
			}
		}
		srevclust.close();
	}

	//interval_map<GenomePos, int> xIntv;
	//interval_map<GenomePos, int> yIntv;
	timing.Tick("Reverse-diag-rough-clusters");

	for (int c = 0; c < split_roughClusters.size(); c++) {
		int rci = genome.header.Find(split_roughClusters[c].tStart);
		StoreFineClusters(rci, forMatches, clusters, opts, split_roughClusters[c].start, split_roughClusters[c].end, genome, read.length, forwardStrand, c);
	}

	for (int c = 0; c < split_revroughClusters.size(); c++) {
		int rci = genome.header.Find(split_revroughClusters[c].tStart);
		StoreFineClusters(rci, revMatches, clusters, opts, split_revroughClusters[c].start, split_revroughClusters[c].end, genome, read.length, reverseStrand, c);
	}

	//cerr << "clusters.size(): " <<  clusters.size() << endl;
	timing.Tick("Fine-clusters");
	//
	// Split clusters on x and y coordinates, vector<Cluster> splitclusters, add a member for each splitcluster to specify the original cluster it comes from
	//
	// INPUT: vector<Cluster> clusters   OUTPUT: vector<Cluster> splitclusters with member--"coarse" specify the index of the original cluster splitcluster comes from
	if (clusters.size() == 0) {	
		return 0; // This read cannot be mapped to the genome; 
	}

	// if (opts.dotPlot and read.name == "chr1/tig00000001|arrow") {
	// 	ofstream pcclust("clusters-coarse-pre-remove.tab");
	// 	for (int m = 0; m < clusters.size(); m++) {
	// 		if (clusters[m].strand == 0) {
	// 			pcclust << clusters[m].qStart << "\t"
	// 						<< clusters[m].tStart << "\t"
	// 						<< clusters[m].qEnd   << "\t"
	// 						<< clusters[m].tEnd   << "\t"
	// 						<< m << "\t"
	// 					  	<< genome.header.names[clusters[m].chromIndex]<< "\t"
	// 						<< clusters[m].strand << "\t"
	// 						<< clusters[m].outerCluster << "\t"
	// 						<< clusters[m].clusterIndex << "\t"
	// 						<< clusters[m].matches.size() << endl;
	// 		}
	// 		else {
	// 			pcclust << clusters[m].qStart << "\t"
	// 						<< clusters[m].tEnd   << "\t"
	// 						<< clusters[m].qEnd   << "\t"
	// 						<< clusters[m].tStart << "\t"
	// 						<< m << "\t"
	// 					 	<< genome.header.names[clusters[m].chromIndex]<< "\t"
	// 						<< clusters[m].strand << "\t"
	// 						<< clusters[m].outerCluster << "\t"
	// 						<< clusters[m].clusterIndex << "\t"
	// 						<< clusters[m].matches.size() << endl;
	// 		}
	// 	}
	// 	pcclust.close();
	// }

	if (opts.dotPlot) {
		ofstream cpclust("clusters-coarse-pre-remove.tab");
		for (int m = 0; m < clusters.size(); m++) {
			for (int n = 0; n < clusters[m].matches.size(); n++) {
				if (clusters[m].strand == 0) {
					cpclust << clusters[m].matches[n].first.pos << "\t"
						  << clusters[m].matches[n].second.pos << "\t"
						  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].second.pos + opts.globalK << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].strand << endl;
				}
				else {
					cpclust << clusters[m].matches[n].first.pos << "\t"
						  << clusters[m].matches[n].second.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].second.pos << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].strand << endl;
				}				
			}
		}
		cpclust.close();
	}
	
	ClusterOrder fineClusterOrder(&clusters, 1);  // has some bug (delete clusters which should be kept -- cluster9_scaffold_58.fasta and cluster18_contig_234.fasta)
	RemoveOverlappingClusters(clusters, fineClusterOrder.index, opts);
	if (opts.dotPlot) {
		ofstream clust("clusters-post-remove.tab");
		

	}

	/*
	for (int co=0; co < clusters.size(); co++) {
		int idx=fineClusterOrder2.index[co];
		cout << co << "\t" << clusters[idx].size() << "\t" << clusters[idx].qStart << "\t" << clusters[idx].qEnd << endl;
	}
	*/

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
		ofstream clust("clusters_coarse.tab");
		for (int m = 0; m < clusters.size(); m++) {
				if (clusters[m].strand == 0) {
					clust << clusters[m].qStart << "\t"
						  << clusters[m].tStart << "\t"
						  << clusters[m].qEnd << "\t"
						  << clusters[m].tEnd << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].rank << "\t"
						  << clusters[m].strand << endl;
				}
				else {
					clust << clusters[m].qStart << "\t"
						  << clusters[m].tEnd << "\t"
						  << clusters[m].qEnd << "\t"
						  << clusters[m].tStart << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].rank << "\t"
						  << clusters[m].strand << endl;
				}				
		}
		clust.close();
	}

	vector<Cluster> splitclusters;
	SplitClusters(clusters, splitclusters);
	DecideSplitClustersValue(clusters, splitclusters, opts);

	timing.Tick("Split");

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
					  << genome.header.names[clusters[splitclusters[m].coarse].chromIndex]<< "\t"
					  << splitclusters[m].strand << endl;
			}
			else {
				clust << splitclusters[m].qStart << "\t" 
					  << splitclusters[m].tEnd << "\t"
					  << splitclusters[m].qEnd   << "\t"
					  << splitclusters[m].tStart   << "\t"
					  << m << "\t"
					  << splitclusters[m].coarse << "\t"
					  << genome.header.names[clusters[splitclusters[m].coarse].chromIndex]<< "\t"
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
						  << genome.header.names[clusters[splitclusters[m].coarse].chromIndex]<< "\t"
						  << splitclusters[m].strand << "\t"
						  << splitclusters[m].Val << endl;
				}
				else {
					clust << splitclusters[m].qStart << "\t" 
						  << splitclusters[m].tEnd << "\t"
						  << splitclusters[m].qEnd   << "\t"
						  << splitclusters[m].tStart   << "\t"
						  << m << "\t"
						  << splitclusters[m].coarse << "\t"
					 	 << genome.header.names[clusters[splitclusters[m].coarse].chromIndex]<< "\t"
						  << splitclusters[m].strand << "\t"
						  << splitclusters[m].Val << endl;
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
	//cerr << "clusters.size(): " << clusters.size() << endl;
	//cerr << "splitclusters.size(): " << splitclusters.size() << endl;
	vector<Primary_chain> Primary_chains;
	float rate = opts.anchor_rate;
	if (splitclusters.size()/clusters.size() > 20) { // mapping to repetitive region
		rate = rate / 2.0;
	}
	//cerr << "rate: " << rate << endl;

	SparseDP (splitclusters, Primary_chains, opts, LookUpTable, read, rate);
	// for (int p = 0; p < Primary_chains.size(); p++) {
	// 	for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
	// 		for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
	// 			int ph = Primary_chains[p].chains[h].ch[c];
	// 			splitclusters[ph].used = 1;
	// 		}
	// 	}
	// }
	// SparseDP (splitclusters, Primary_chains, opts, LookUpTable, read, rate);

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
							  << splitclusters[ph].Val << "\t"
							  << p << "\t"
							  << h << "\t"
							  << Primary_chains[p].chains[h].ch[c] << "\t"
							  << splitclusters[ph].strand << "\t"
							  << splitclusters[ph].NumofAnchors << "\t"
							  << splitclusters[ph].coarse << endl;
					} 
					else {
						clust << splitclusters[ph].qStart << "\t" 
							  << splitclusters[ph].tEnd << "\t"
							  << splitclusters[ph].qEnd   << "\t"
							  << splitclusters[ph].tStart  << "\t"
							  << splitclusters[ph].Val << "\t"							  
							  << p << "\t"
							  << h << "\t"
							  << Primary_chains[p].chains[h].ch[c] << "\t"
							  << splitclusters[ph].strand << "\t"
							  << splitclusters[ph].NumofAnchors << "\t"						
							  << splitclusters[ph].coarse << endl;
					}
				}
			}
		}
		clust.close();
	}	

	switchindex(splitclusters, Primary_chains, clusters, genome);
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
				//" value: " << Primary_chains[p].chains[h].value << " # of Anchors: " << Primary_chains[p].chains[h].NumOfAnchors << " tStart: " <<  Primary_chains[p].chains[h].tStart << endl;
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
							  << genome.header.names[clusters[ph].chromIndex]<< "\t"
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
							  << genome.header.names[clusters[ph].chromIndex]<< "\t"
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
	// Build local index for refining alignments.
	//
	LocalIndex forwardIndex(glIndex);
	LocalIndex reverseIndex(glIndex);
	LocalIndex *localIndexes[2] = {&forwardIndex, &reverseIndex};
	forwardIndex.IndexSeq(read.seq, read.length);
	reverseIndex.IndexSeq(readRC, read.length); 
	//
	// Set the parameters for merging anchors and 1st SDP
	//
	Options smallOpts=opts;
	smallOpts.coefficient=opts.predefined_coefficient; // used to be 12
	Options tinyOpts=smallOpts;
	tinyOpts.globalMaxFreq=3;
	tinyOpts.maxDiag=5;
	tinyOpts.minDiagCluster=2;
	//
	// Decide whether the number of anchors in each Cluster is enough to skip refining step;
	//
	bool sparse = 0;
	for (int p = 0; p < clusters.size(); p++) {
		//if (((float)(clusters[p].matches.size())/(clusters[p].qEnd - clusters[p].qStart)) < 0.001) sparse = 1;
		if (((float)(clusters[p].Val)/(clusters[p].qEnd - clusters[p].qStart)) < 0.01 and read.length <= 50000) {
			sparse = 1;
		}
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
		smallOpts.coefficient+=3; // used to be 15
		smallOpts.globalMaxFreq=6;
		smallOpts.cleanMaxDiag=10;// used to be 25
		smallOpts.maxDiag=50;
		smallOpts.maxGapBtwnAnchors=100; // used to be 200 // 200 seems a little bit large
		smallOpts.minDiagCluster=3; // used to be 3
		tinyOpts.globalK=smallOpts.globalK-3;

		REFINEclusters(clusters, refinedclusters, genome, read,  glIndex, localIndexes, smallOpts, opts);
		//cerr << "refine cluster done!" << endl;
		// refinedclusters have GenomePos, chromIndex, coarse, matches, strand, refinespace;
		for (int s = 0; s < clusters.size(); s++) {
			RefinedClusters[s] = &refinedclusters[s];
		}
		clusters.clear();
	}
	else {
		tinyOpts.globalK=smallOpts.globalK-8;//-8
		//cerr << "tinyOpts.globalK: " << tinyOpts.globalK << endl;
		//cerr << "smallOpts.globalK: " << smallOpts.globalK << endl;
		//cerr << "opts.globalK: " << opts.globalK << endl;
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
	//cerr << "read.name: " << read.name << endl;

	if (opts.dotPlot) {
		ofstream clust("RefinedClusters.tab", std::ofstream::app);

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
	AffineAlignBuffers buff;
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
			//cerr << "refinement done!" << endl;
			//
			// Do linear extension for each anchors and avoid overlapping locations;
			// INPUT: RefinedClusters; OUTPUT: ExtendClusters;
			// NOTICE: ExtendClusters have members: strand, matches, matchesLengths, GenomePos, chromIndex;
			//
			vector<Cluster> ExtendClusters(Primary_chains[p].chains[h].ch.size());
			LinearExtend(RefinedClusters, ExtendClusters, Primary_chains[p].chains[h].ch, smallOpts, genome, read);
			TrimOverlappedAnchors(ExtendClusters);
			//cerr << "LinearExtend done!" << endl;

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
								  << ep << "\t"
								  << genome.header.names[ExtendClusters[ep].chromIndex]<< "\t"
								  << ExtendClusters[ep].strand << "\t"
								  << eh << endl;
						}
						else {
							clust << ExtendClusters[ep].matches[eh].first.pos << "\t"
								  << ExtendClusters[ep].matches[eh].second.pos + ExtendClusters[ep].matchesLengths[eh] << "\t"
								  << ExtendClusters[ep].matches[eh].first.pos + ExtendClusters[ep].matchesLengths[eh] << "\t"
								  << ExtendClusters[ep].matches[eh].second.pos<< "\t"
								  << ep << "\t"
								  << genome.header.names[ExtendClusters[ep].chromIndex]<< "\t"
								  << ExtendClusters[ep].strand << "\t"
								  << eh << endl;					
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
			//cerr << "splitchains.size(): " << splitchains.size()  << endl;
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
				RemovePairedIndels(finalchain);
				RemoveSpuriousAnchors(finalchain, smallOpts);
				//cerr << "2nd SDP done!" << endl;
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
				SeparateChainByStrand(finalchain, finalSeperateChain, ExtendClusters); 
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
					Alignment *alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, read.name, str, read.qual, genome.seqs[chromIndex],  
																		 genome.lengths[chromIndex], genome.header.names[chromIndex], chromIndex, 0); 
					alignment->NumOfAnchors = Primary_chains[p].chains[h].NumOfAnchors;
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
								//
								// Check the distance between two anchors. If the distance is too large, refine anchors and apply 3rd SDP;
								// Otherwise apply linear alignment. 
								RefinedAlignmentbtwnAnchors(cur, next, str, chromIndex, finalchain, alignment, read, genome,
															strands, scoreMat, pathMat, tinyOpts, genomeThres, buff, LookUpTable);
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
															strands, scoreMat, pathMat, tinyOpts, genomeThres, buff, LookUpTable);
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
			alignments.back().SetFromSegAlignment(smallOpts);
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

	timing.Tick("Final");
	if (opts.storeTiming) {
		for (int a=0; a < (int) min(alignmentsOrder.size(), opts.PrintNumAln); a++){
			for (int s = 0; s < alignmentsOrder[a].SegAlignment.size(); s++) {
				alignmentsOrder[a].SegAlignment[s]->runtime=timing.Elapsed();
				alignmentsOrder[a].SegAlignment[s]->nanchors=allMatches.size();
			}
		}
	}
	if (alignmentsOrder.size() > 0 and alignmentsOrder[0].SegAlignment.size() > 0) {
		for (int a=0; a < (int) min(alignmentsOrder.size(), opts.PrintNumAln); a++){
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
