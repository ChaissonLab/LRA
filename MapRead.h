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
#include "TupleOps.h"
#include "SparseDP.h"
#include "SparseDP_Forward.h"
#include "Chain.h"
#include "overload.h"
#include "LinearExtend.h"
#include "SplitClusters.h"
#include "Timing.h"
#include "ClusterRefine.h"
#include "IndelRefine.h"
#include "LocalRefineAlignment.h"
#include "Map_lowacc.h"
#include "Map_highacc.h"

#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cmath>	// std::log 
#include <sstream>
#include <thread>
#include <climits>
#include <map>

using namespace std;


class SortClusterBySize {
public:
	bool operator()(const Cluster &a, const Cluster &b) {
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
void SwapReadCoordinates(vector<T> &matches, GenomePos readLength, GenomePos kmer){

	for (int i=0; i < matches.size(); i++) {
		matches[i].first.pos = readLength - (matches[i].first.pos+ kmer);
	}
}

// void ReverseClusterStrand(Read &read, Genome &genome, Options &opts, vector<Cluster> &clusters) {
// 	for (int c = 0; c < clusters.size(); c++) {
// 			SwapStrand(read, opts, clusters[c].matches);
// 			clusters[c].strand = 1;
// 	}
// }

// void SetClusterStrand(Read &read, Genome &genome, Options &opts, 
// 											vector<Cluster> &clusters) {
// 	for (int c = 0; c < clusters.size(); c++) {
// 		clusters[c].strand = SetStrand(read, genome, opts, clusters[c].matches);
// 		if (clusters[c].strand == 1) {
// 			SwapStrand(read, opts, clusters[c].matches);
// 		}
// 	}
// }


void 
SeparateMatchesByStrand(Read &read, Genome &genome, int k, vector<GenomePair> &allMatches,  vector<GenomePair> &forMatches,
								vector<GenomePair> &revMatches, string &baseName) {
	//
	// A value of 0 implies forward strand match.
	//
	vector<bool> strand(allMatches.size(), 0);
	int nForward=0;
	for (int i=0; i < allMatches.size(); i++) {
		int readPos = allMatches[i].first.pos; 
		uint64_t refPos = allMatches[i].second.pos;
		char *genomePtr = genome.GlobalIndexToSeq(refPos);
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
	revMatches.resize(allMatches.size() - nForward);
	int i = 0,r = 0,f = 0;
	for (i = 0,r = 0,f = 0; i < allMatches.size(); i++) {
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


int 
MapRead(const vector<float> & LookUpTable, Read &read, Genome &genome, vector<GenomeTuple> &genomemm, LocalIndex &glIndex, Options &opts, 
				ostream *output, ostream *svsigstrm, Timing &timing, IndelRefineBuffers &indelRefineBuffers, pthread_mutex_t *semaphore=NULL) {
	read.unaligned = 0;
	string baseName = read.name;
	for (int i=0; i < baseName.size(); i++) {	
		if (baseName[i] == '/') baseName[i] = '_';	
		if (baseName[i] == '|') baseName[i] = '_';
	}
	vector<GenomeTuple> readmm; // readmm stores minimizers
	vector<pair<GenomeTuple, GenomeTuple> > allMatches, forMatches, revMatches;
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
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, allOpts.globalK, allOpts.globalW, readmm, true);	
		// StoreMinimizers_noncanonical<GenomeTuple, Tuple>(read.seq, read.length, allOpts.globalK, allOpts.globalW, readmm, true);			
	}
	else {
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, opts.globalK, opts.globalW, readmm, true);
		// StoreMinimizers_noncanonical<GenomeTuple, Tuple>(read.seq, read.length, opts.globalK, opts.globalW, readmm, true);
	}
	timing.Tick("Store minimizers");
	sort(readmm.begin(), readmm.end()); //sort kmers in readmm(minimizers)
	timing.Tick("Sort minimizers");
	//
	// Add matches between the read and the genome.
	//
	CompareLists<GenomeTuple, Tuple>(readmm, genomemm, allMatches, opts, true);
	timing.Tick("CompareLists");

	if (opts.dotPlot) {
		ofstream clust("all-matches.dots");
		for (int m = 0; m < allMatches.size(); m++) {
			clust << allMatches[m].first.pos << "\t" << allMatches[m].second.pos
						<< "\t" << allMatches[m].first.pos + opts.globalK << "\t" 
						<< allMatches[m].second.pos+ opts.globalK << endl;
		}
		clust.close();
	}

	SeparateMatchesByStrand(read, genome, opts.globalK, allMatches, forMatches, revMatches, baseName);
	allMatches.clear();
	if (forMatches.size() == 0 and revMatches.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
	} 
	if (opts.dotPlot) {
		ofstream fclust("for-matches_original.dots");
		for (int m = 0; m < forMatches.size(); m++) {
			fclust << forMatches[m].first.pos << "\t" << forMatches[m].second.pos << "\t" << opts.globalK + forMatches[m].first.pos << "\t"
					<< forMatches[m].second.pos + opts.globalK << "\t" << m << endl;
		}
		fclust.close();
		ofstream rclust("rev-matches_original.dots");
		for (int m=0; m < revMatches.size(); m++) {			
			rclust << revMatches[m].first.pos << "\t" << revMatches[m].second.pos + opts.globalK << "\t" << opts.globalK + revMatches[m].first.pos  << "\t"
					 << revMatches[m].second.pos << "\t" << m << endl;
		}
		rclust.close();
	}		

	if (opts.bypassClustering) { 
		return MapRead_lowacc(forMatches, revMatches, LookUpTable, read, genome, genomemm, glIndex, opts, output, svsigstrm, 
					timing, indelRefineBuffers, strands, readRC, semaphore);
	}
	else { 
		return MapRead_highacc(forMatches, revMatches, LookUpTable, read, genome, genomemm, glIndex, opts, output, svsigstrm, 
					timing, indelRefineBuffers, strands, readRC, semaphore);
	}

	// /*
	// if (semaphore != NULL ) {
	// 	pthread_mutex_unlock(semaphore);
	// }
	// */
	// //
	// // Done with one read. Clean memory.
	// //
	// delete[] readRC;
	// for (int a = 0; a < alignments.size(); a++) {
	// 	for (int s = 0; s < alignments[a].SegAlignment.size(); s++) {
	// 		delete alignments[a].SegAlignment[s];
	// 	}
	// }
	// //read.Clear();
	// if (alignments.size() > 0) return 1;
	// return 0;
}

#endif
