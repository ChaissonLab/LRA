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
// #include "GlobalChain.h"
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

void SimpleMapQV(AlignmentsOrder &alignmentsOrder, Read &read) {
	static const float q_coef = 40.0f;
	int len = alignmentsOrder.size();
	for (int r = 0; r < len; r++) {
		if (r == 0 and len == 1) {
			alignmentsOrder[r].mapqv = 60;
		}
		else if (r == 0 and len > 1) {
			// float x = (0.5f * alignmentsOrder[r + 1].value + 0.4f * alignmentsOrder[r + 1].SecondSDPValue + 0.1f * alignmentsOrder[r + 1].FirstSDPValue) / 
			// 			(0.5f * alignmentsOrder[r].value + 0.4f * alignmentsOrder[r].SecondSDPValue + 0.1f * alignmentsOrder[r].FirstSDPValue);
			//float x = (0.8f * alignmentsOrder[r + 1].value + 0.2f * alignmentsOrder[r + 1].FirstSDPValue) / (0.8f * alignmentsOrder[r].value + 0.2f * alignmentsOrder[r].FirstSDPValue);
			float x = (alignmentsOrder[r + 1].value) / (alignmentsOrder[r].value);
			float pen_cm_1 = (alignmentsOrder[r].NumOfAnchors0 > 20? 1.0f : 0.05f ) * alignmentsOrder[r].NumOfAnchors0;
			pen_cm_1 = (alignmentsOrder[r].NumOfAnchors0 >= 5? 1.0f : 0.1f ) * pen_cm_1;
			// float pen_cm_1 = (alignmentsOrder[r].NumOfAnchors0 > 10? 1.0f : 0.1f ) * alignmentsOrder[r].NumOfAnchors0;
			// float pen_cm_2 = (alignmentsOrder[r].NumOfAnchors1 > 100? 1.0f : 0.01f) * alignmentsOrder[r].NumOfAnchors1;
			float mapq = (int)(pen_cm_1 * q_coef * (1.0f - x) * logf(alignmentsOrder[r].value));
			mapq -= (int)(4.343f * logf(len) + .499f);
			mapq = mapq > 0? mapq : 0;
			alignmentsOrder[r].mapqv = mapq < 60? mapq : 60;
			if (r == 0 && len == 2 && alignmentsOrder[r].mapqv == 0) alignmentsOrder[r].mapqv = 1;
		}
		else {
			alignmentsOrder[r].mapqv = 0;
		}
		alignmentsOrder[r].SetMapqv();
		//cerr << "alignmentsOrder[" << r << "].mapq: " << (int) alignmentsOrder[r].mapqv << endl;
	}
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
void SwapReadCoordinates(vector<T> &matches, GenomePos readLength, GenomePos kmer){

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

// template<typename T>
// void 
// UpdateBoundaries(T &matches, GenomePos &qStart, GenomePos &qEnd, GenomePos &tStart, GenomePos &tEnd) {
// 	for (int i =0; i< matches.size(); i++) {
// 		qStart=min(qStart, matches[i].first.pos);
// 		qEnd=max(qEnd, matches[i].first.pos);
// 		tStart=min(tStart, matches[i].second.pos);
// 		tEnd=max(tEnd, matches[i].second.pos);
// 	}
// }
// // int nSSE=0;


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
	// For cases like chain = {22, 125, 19, 125, 16, 17, 125, 57, 125}, compress multiple 125 to one 125;
	//
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			std::map<int, int> appeartimes;
			std::map<int, int> start_pos;
			std::map<int, int> end_pos;
			for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++){
				int ats = Primary_chains[p].chains[h].ch[c];
				if (appeartimes.count(ats) > 0) {
					appeartimes[ats] += 1;
					end_pos[ats] = c + 1;
				}
				else {
					appeartimes[ats] = 1;
					start_pos[ats] = c;
					end_pos[ats] = c + 1;
				}
			}	

			if (start_pos.size() == 0) {continue;}

			vector<tuple<int, int> > start_end;
			for (std::map<int,int>::iterator ait=start_pos.begin(); ait!=start_pos.end(); ++ait) {
				if (end_pos[ait->first] > ait->second + 1) {
					start_end.push_back(make_tuple(ait->second, end_pos[ait->first]));					
				}
			}
			sort(start_end.begin(), start_end.end());

			vector<unsigned int> newch;
			vector<bool> newlink;
			int ste = 0;
			int nc = 0;
			int cur_ste_end = 0;
			while (ste < start_end.size()) {
				while (nc <= get<0>(start_end[ste])) {
					newch.push_back(Primary_chains[p].chains[h].ch[nc]);
					if (newch.size() > 1) {
						newlink.push_back(Primary_chains[p].chains[h].link[nc-1]);
					}
					nc++;					
				}
				nc = get<1>(start_end[ste]);
				//newch.push_back(Primary_chains[p].chains[h].ch[nc]); // add the end
				ste++;
			}
			while (nc < Primary_chains[p].chains[h].ch.size()) {
				newch.push_back(Primary_chains[p].chains[h].ch[nc]);
				if (newch.size() > 1) {
					newlink.push_back(Primary_chains[p].chains[h].link[nc-1]);
				}
				nc++;
			}
			Primary_chains[p].chains[h].ch = newch;
			int snch = newch.size(); int snlink = newlink.size();
			Primary_chains[p].chains[h].ch.resize(snch);
			Primary_chains[p].chains[h].link = newlink;
			Primary_chains[p].chains[h].link.resize(snlink);
		}
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

int MapRead(const vector<float> & LookUpTable, Read &read, Genome &genome, vector<GenomeTuple> &genomemm, LocalIndex &glIndex, Options &opts, 
			ostream *output, ostream *svsigstrm, Timing &timing, pthread_mutex_t *semaphore=NULL) {
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
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, allOpts.globalK, allOpts.globalW, readmm, true, false);			
	}
	else {
		StoreMinimizers<GenomeTuple, Tuple>(read.seq, read.length, opts.globalK, opts.globalW, readmm, true, false);
	}
	timing.Tick("Store minimizers");
	sort(readmm.begin(), readmm.end()); //sort kmers in readmm(minimizers)
	timing.Tick("Sort minimizers");
	// if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
	// 	stringstream outNameStrm;
	// 	outNameStrm << "read.minimizers";
	// 	ofstream baseDots(outNameStrm.str().c_str());
	// 	for (int m=0; m < readmm.size(); m++) {
	// 		baseDots << readmm[m].t << "\t"
	// 				 << readmm[m].pos << "\t" 
	// 				 << readmm[m].pos + opts.globalK << endl;				
	// 	}
	// 	baseDots.close();
	// }
	// cerr << "maxGap: " << opts.maxGap << " maxDiag: " << opts.maxDiag << " cleanMaxDiag: " << opts.cleanMaxDiag << " minDiagCluster: " << opts.minDiagCluster
	// 	<< " minClusterSize: " << opts.minClusterSize << " opts.firstcoefficient: " << opts.firstcoefficient  << " opts.secondcoefficient: " << opts.secondcoefficient << endl;

	//
	// Add matches between the read and the genome.
	//
	CompareLists<GenomeTuple, Tuple>(readmm, genomemm, allMatches, opts, true);
	timing.Tick("CompareLists");

	// if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
	// 	ofstream clust("all-matches.dots");
	// 	for (int m = 0; m < allMatches.size(); m++) {
	// 		clust << allMatches[m].first.pos << "\t" << allMatches[m].second.pos
	// 					<< "\t" << allMatches[m].first.pos + opts.globalK << "\t" 
	// 					<< allMatches[m].second.pos+ opts.globalK << endl;
	// 	}
	// 	clust.close();
	// }

	SeparateMatchesByStrand(read, genome, opts.globalK, allMatches, forMatches, revMatches, baseName);
	vector<Cluster> clusters;
	MatchesToFineClusters(forMatches, clusters, genome, read, opts, timing);
	MatchesToFineClusters(revMatches, clusters, genome, read, opts, timing, 1);

	if (opts.CheckTrueIntervalInFineCluster) {
		CheckTrueIntervalInFineCluster(clusters, read.name, genome);
	}
	//cerr << "clusters.size(): " <<  clusters.size() << endl;

	if (clusters.size() == 0) {
		cerr << "unmapped " << read.name << endl;	
		return 0; // This read cannot be mapped to the genome; 
	}

	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
		ofstream cpclust("clusters-pre-remove.tab");
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

	allMatches.clear();
	forMatches.clear(); 
	revMatches.clear();

	// cerr << "before removing clusters.size(): " << clusters.size() << endl;
	// cerr << "before removing splitclusters.size(): " << splitclusters.size() << endl;

	if (clusters.size() >= 4000) {
		ClusterOrder fineClusterOrder(&clusters, 1);  // has some bug (delete clusters which should be kept -- cluster9_scaffold_58.fasta and cluster18_contig_234.fasta)
		RemoveOverlappingClusters(clusters, fineClusterOrder.index, opts);
	}

	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
		ofstream clust("clusters-post-remove.tab");
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

		ofstream sclust("clusters_coarse.tab");
		for (int m = 0; m < clusters.size(); m++) {
				if (clusters[m].strand == 0) {
					sclust << clusters[m].qStart << "\t"
						  << clusters[m].tStart << "\t"
						  << clusters[m].qEnd << "\t"
						  << clusters[m].tEnd << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].rank << "\t"
						  << clusters[m].strand << endl;
				}
				else {
					sclust << clusters[m].qStart << "\t"
						  << clusters[m].tEnd << "\t"
						  << clusters[m].qEnd << "\t"
						  << clusters[m].tStart << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].rank << "\t"
						  << clusters[m].strand << endl;
				}				
		}
		sclust.close();
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
	if (clusters.size() == 0) {
		cerr << "unmapped " << read.name << endl;	
		return 0; // This read cannot be mapped to the genome;
	}

	//
	// Split clusters on x and y coordinates, vector<Cluster> splitclusters, add a member for each splitcluster to specify the original cluster it comes from.
	// INPUT: vector<Cluster> clusters   OUTPUT: vector<Cluster> splitclusters with member--"coarse" specify the index of the original cluster splitcluster comes from
	//
	vector<Cluster> splitclusters;
	SplitClusters(clusters, splitclusters);
	DecideSplitClustersValue(clusters, splitclusters, opts);
	if (splitclusters.size() == 0) {
		cerr << "unmapped " << read.name << endl;	
		return 0; // This read cannot be mapped to the genome;		
	}
	timing.Tick("SplitClusters");

	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
		ofstream sclust("splitclusters-decideval.tab");
		for (int m = 0; m < splitclusters.size(); m++) {
			if (splitclusters[m].Val !=0) {
				if (splitclusters[m].strand == 0) {
					sclust << splitclusters[m].qStart << "\t" 
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
					sclust << splitclusters[m].qStart << "\t" 
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
		sclust.close();
	}

	//
	// Apply SDP on splitclusters. Based on the chain, clean clusters to make it only contain clusters that are on the chain.   --- vector<Cluster> clusters
	// class: chains: vector<chain> chain: vector<vector<int>>     Need parameters: PrimaryAlgnNum, SecondaryAlnNum
	// NOTICE: chains in Primary_chains do not overlap on Cluster
	//
	vector<Primary_chain> Primary_chains;
	float rate = opts.anchor_rate;
	if (splitclusters.size()/clusters.size() > 20) { // mapping to repetitive region
		rate = rate / 2.0;
	}
	cerr << "clusters.size(): " << clusters.size() << " splitclusters.size(): " << splitclusters.size() << " rate: " << rate << " read.name: " << read.name << endl;

	SparseDP (splitclusters, Primary_chains, opts, LookUpTable, read, rate);
	if (Primary_chains.size() == 0) {
		cerr << "unmapped " << read.name << endl;	
		return 0; // This read cannot be mapped to the genome;			
	}

	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
							  << splitclusters[ph].NumofAnchors0 << "\t"
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
							  << splitclusters[ph].NumofAnchors0 << "\t"						
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

	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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

	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
	//smallOpts.secondcoefficient=opts.predefined_coefficient; // used to be 12
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
		smallOpts.secondcoefficient+=3; // used to be 15
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
	timing.Tick("Refine_clusters");

	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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

	int SizeRefinedClusters = 0;
	for (int p = 0; p < RefinedClusters.size(); p++) {
		SizeRefinedClusters += RefinedClusters[p]->matches.size();
	}
	if (SizeRefinedClusters == 0) {
		cerr << "unmapped " << read.name <<  endl;	
		return 0; // This read is not mapped!
	}

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

			if (Primary_chains[p].chains[h].ch.size() == 0) {
				cerr << "unmapped " << read.name <<  endl;	
				continue;
			}
			alignments.resize(alignments.size() + 1);			
			RefineBtwenClusters_chain(Primary_chains, alignments, RefinedClusters, genome, read, smallOpts, p, h, strands);
			timing.Tick("refine_btwnclusters");
			vector<Cluster> ExtendClusters(Primary_chains[p].chains[h].ch.size());
			LinearExtend_chain(Primary_chains, ExtendClusters, RefinedClusters, smallOpts, genome, read, p, h);
			timing.Tick("LinearExtend");

			int SizeExtendClusters = 0;
			for (int ep = 0; ep < ExtendClusters.size(); ep++) {
				SizeExtendClusters += ExtendClusters[ep].matches.size();
			}	

			if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
								  << ExtendClusters[ep].strand << "\t"
								  << ep << endl;
						}
						else {
							clust << ExtendClusters[ep].matches[eh].first.pos << "\t"
								  << ExtendClusters[ep].matches[eh].second.pos + ExtendClusters[ep].matchesLengths[eh] << "\t"
								  << ExtendClusters[ep].matches[eh].first.pos + ExtendClusters[ep].matchesLengths[eh] << "\t"
								  << ExtendClusters[ep].matches[eh].second.pos<< "\t"
								  << h << "\t"
								  << genome.header.names[ExtendClusters[ep].chromIndex]<< "\t"
								  << ExtendClusters[ep].strand << "\t"
								  << ep << endl;					
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
			SparseDP_and_RefineAlignment_btwn_anchors(Primary_chains, splitchains, ExtendClusters, alignments, smallOpts, 
									LookUpTable, read, strands, p, h, genome, LSC, tinyOpts, buff, svsigstrm);
			timing.Tick("SparseDP_on_anchors and local refinment");


			// for (int st = 0; st < splitchains.size(); st++) {
			// 	//
			// 	// Apply SparseDP on extended anchors on every splitchain;
			// 	// INPUT: vector<unsigned int> splitchain, vector<Cluster> ExtendClusters; OUTPUT: FinalChain finalchain;
			// 	//
			// 	FinalChain finalchain(&ExtendClusters);
			// 	SparseDP(splitchains[st], ExtendClusters, finalchain, smallOpts, LookUpTable, read, timing);
			// 	//
			// 	// RemoveSpuriousAnchors and RemovePairedIndels; 
			// 	//
			// 	// RemovePairedIndels(finalchain); (This is deleting lots of good anchors)
			// 	// RemoveSpuriousAnchors(finalchain, smallOpts);

			// 	//cerr << "2nd SDP done!" << endl;
			// 	if (finalchain.size() == 0) continue; // cannot be mapped to the genome!
			// 	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
			// 		ofstream clust("SparseDP.tab", std::ofstream::app);
			// 		for (int ep = 0; ep < finalchain.chain.size(); ep++) {
			// 			if (finalchain.strand(ep) == 0) {
			// 				clust << finalchain[ep].first.pos << "\t"
			// 					  << finalchain[ep].second.pos << "\t"
			// 					  << finalchain[ep].first.pos + finalchain.length(ep) << "\t"
			// 					  << finalchain[ep].second.pos + finalchain.length(ep) << "\t"
			// 					  << h << "\t"
			// 					  << finalchain.ClusterNum(ep) << "\t"
			// 					  << finalchain.strand(ep) << endl;
			// 			}
			// 			else {
			// 				clust << finalchain[ep].first.pos << "\t"
			// 					  << finalchain[ep].second.pos + finalchain.length(ep) << "\t"
			// 					  << finalchain[ep].first.pos + finalchain.length(ep) << "\t"
			// 					  << finalchain[ep].second.pos << "\t"
			// 					  << h << "\t"
			// 					  << finalchain.ClusterNum(ep) << "\t"
			// 					  << finalchain.strand(ep) << endl;					
			// 			}
			// 		}
			// 		clust.close();
			// 	}	
			// 	//
			// 	// If there are inversions in the path, then seperate finalchain into several parts. 
			// 	// In each part anchors are in the same direction, which is for better manipulation;
			// 	// INPUT: FinalChain finalchain, OUTPUT: vector<int> finalchainSeperateChain;
			// 	//
			// 	vector<vector<int>> finalSeperateChain;
			// 	SeparateChainByStrand(finalchain, finalSeperateChain, ExtendClusters); 
			// 	int LFC = LargestFinalSeperateChain(finalSeperateChain);
			// 	//
			// 	// Refine and store the alignment; NOTICE: when filling in the space between two adjacent anchors, 
			// 	// the process works in forward direction, so we need to flip the small matches
			// 	// found in the spaces before insert them into the alignment if necessary;
			// 	// INPUT: vector<int> finalchainSeperateStrand; OUTPUT: Alignment *alignment;
			// 	//
			// 	int ad = alignments.back().SegAlignment.size();
			// 	for (int fsc = 0; fsc < finalSeperateChain.size(); fsc++) {

			// 		assert(finalSeperateChain[fsc].size() == 2);

			// 		int start = finalSeperateChain[fsc][0];
			// 		int end = finalSeperateChain[fsc][1];
			// 		assert(start < end); 
			// 		bool str = finalchain.strand(start);
			// 		int cln = finalchain.ClusterNum(start);
			// 		int chromIndex = ExtendClusters[cln].chromIndex;	
			// 		Alignment *alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, 
			// 											 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
			// 											 genome.header.names[chromIndex], chromIndex, 0); 
			// 		alignment->NumOfAnchors0 = Primary_chains[p].chains[h].NumOfAnchors0;
			// 		alignment->NumOfAnchors1 = finalchain.size();
			// 		alignment->SecondSDPValue = finalchain.SecondSDPValue;
			// 		alignments.back().SegAlignment.push_back(alignment);
			// 		vector<int> scoreMat;
			// 		vector<Arrow> pathMat;
			// 		//
			// 		// Set the secondary or supplymentary flag for the alignment; 
			// 		// h == 0 ---> primary alignment, otherwise secondary alignment;
			// 		// st > 0 ----> supplymentary alignment; fsc > 0 ----> supplymeantary alignment;
			// 		//
			// 		if (h > 0) alignment->ISsecondary = 1;
			// 		if (splitchains.size() > 1 or finalSeperateChain.size() > 1) alignment->split = 1;					
			// 		if (st != LSC or fsc != LFC) alignment->Supplymentary = 1; //if (st > 0 or fsc > 0) alignment->Supplymentary = 1;	

			// 		GenomePos genomeThres = 0;
			// 		bool inv_str = 0, inversion = 0, prev_inv = 0;
			// 		if (str == 0) {
			// 			int fl = end - 1;
			// 			inv_str = 1; prev_inv = end - 1;
			// 			while (fl > start) {
			// 				int cur = fl;
			// 				int next = fl - 1;
			// 				if (!PassgenomeThres(cur, genomeThres, finalchain)) {
			// 					//
			// 					// Check the distance between two anchors. If the distance is too large, refine anchors and apply 3rd SDP;
			// 					// Otherwise apply linear alignment. 
			// 					//
			// 					RefinedAlignmentbtwnAnchors(cur, next, str, inv_str, chromIndex, finalchain, alignments, alignments.back().SegAlignment.back(), 
			// 												read, genome, strands, scoreMat, pathMat, tinyOpts, genomeThres, buff, 
			// 												LookUpTable, inversion, svsigstrm);

			// 					if (inversion == 1) {
			// 						//alignment.FirstSDPValue /= (float)(prev_inv-fl) / (end-start)
			// 						//alignment->NumOfAnchors = prev_inv-fl;
			// 						alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
			// 						Alignment *next_alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, 
			// 												 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
			// 												 genome.header.names[chromIndex], chromIndex, 0); 
			// 						next_alignment->NumOfAnchors1 = fl-start;
			// 						alignments.back().SegAlignment.push_back(next_alignment);
			// 						prev_inv = fl;
			// 						inversion = 0;
			// 					} 
			// 				}
			// 				fl--;
			// 			}
			// 			if (!PassgenomeThres(start, genomeThres, finalchain)) {
			// 				alignments.back().SegAlignment.back()->blocks.push_back(Block(finalchain[start].first.pos, finalchain[start].second.pos, 
			// 													finalchain.length(start)));
			// 			}
			// 		}
			// 		else {
			// 			int fl = start; 
			// 			inv_str = 0; prev_inv = start;
			// 			while (fl < end - 1) {
			// 				int cur = fl;
			// 				int next = fl + 1;
			// 				if (!PassgenomeThres(cur, genomeThres, finalchain)) {
			// 					RefinedAlignmentbtwnAnchors(cur, next, str, inv_str, chromIndex, finalchain, alignments, alignments.back().SegAlignment.back(), 
			// 												read, genome, strands, scoreMat, pathMat, tinyOpts, genomeThres, buff, 
			// 												LookUpTable, inversion, svsigstrm);
			// 					if (inversion == 1) {
			// 						//alignment.FirstSDPValue /= (float)(prev_inv-fl) / (end-start)
			// 						alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
			// 						//alignment->NumOfAnchors = fl-start;
			// 						Alignment *next_alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, 
			// 												 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
			// 												 genome.header.names[chromIndex], chromIndex, 0); 
			// 						next_alignment->NumOfAnchors1 = end-fl;
			// 						alignments.back().SegAlignment.push_back(next_alignment);
			// 						prev_inv = fl;
			// 						inversion = 0;
			// 					} 
			// 				}
			// 				fl++;
			// 			}	
			// 			if (!PassgenomeThres(end - 1, genomeThres, finalchain)) {
			// 				alignments.back().SegAlignment.back()->blocks.push_back(Block(read.length-finalchain[end - 1].first.pos-finalchain.length(end - 1), 
			// 													finalchain[end - 1].second.pos, finalchain.length(end - 1)));
			// 			}
			// 		}
				
			// 		for (int bb = 1; bb < alignments.back().SegAlignment.back()->blocks.size(); bb++) {
			// 			assert(alignments.back().SegAlignment.back()->blocks[bb-1].qPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= alignments.back().SegAlignment.back()->blocks[bb].qPos);
			// 			assert(alignments.back().SegAlignment.back()->blocks[bb-1].tPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= alignments.back().SegAlignment.back()->blocks[bb].tPos);	
			// 		}
			// 		//
			// 		// Set some parameters in Alignment alignment
			// 		//
			// 		alignments.back().SegAlignment.back()->UpdateParameters(str, smallOpts, LookUpTable, svsigstrm, strands);
			// 	}
			// 	//
			// 	// Decide the inversion flag for each segment alignment; seek +,-,+ or -,+,-
			// 	//
			// 	int js = ad+2;
			// 	while (js < alignments.back().SegAlignment.size()) {
			// 		if ((alignments.back().SegAlignment[js-2]->strand == 0 and alignments.back().SegAlignment[js-1]->strand == 1
			// 			and alignments.back().SegAlignment[js]->strand == 0) or ((alignments.back().SegAlignment[js-2]->strand == 1 
			// 			and alignments.back().SegAlignment[js-1]->strand == 0 and alignments.back().SegAlignment[js]->strand == 1))) {

			// 			//
			// 			// inversion cannot be too far (<10,000) from alignments at both sides;
			// 			//
			// 			if (alignments.back().SegAlignment[js-1]->tStart > alignments.back().SegAlignment[js-2]->tEnd + 10000
			// 				or alignments.back().SegAlignment[js]->tStart > alignments.back().SegAlignment[js-1]->tEnd + 10000) {
			// 				js++;
			// 				continue;
			// 			}
			// 			//
			// 			// length requirements; 
			// 			//
			// 			else if (alignments.back().SegAlignment[js]->nm < 500 or alignments.back().SegAlignment[js-1]->nm < 500
			// 					or alignments.back().SegAlignment[js-2]->nm < 40 or alignments.back().SegAlignment[js-2]->nm > 15000) {
			// 				js++;
			// 				continue;
			// 			}
			// 			else if (alignments.back().SegAlignment[js-2]->typeofaln != 3) {
			// 				alignments.back().SegAlignment[js-1]->typeofaln=3;
			// 			}
			// 		}
			// 		js++;
			// 	}
			// 	// timing.Tick("local refinment");
			// }
			// ///////
			for (int s = 0; s < alignments.back().SegAlignment.size(); s++) {
				//
				// (Mark: add the refining code here)
				// alignments.back() is `class SegAlignmentGroup`. `class SegAlignmentGroup` has member: `vector<Alignment*> SegAlignment`, where 
				// each `SegAlignment[s]` is a segment of alignment when there are multiple segments. 
				// So the function refining can be apply to each of the segment. 
				// function refining(alignments.back().SegAlignment[s])
				//
				if (opts.skipBandedRefine == false) {
					IndelRefineAlignment(read, genome, *alignments.back().SegAlignment[s], opts);
				}
				alignments.back().SegAlignment[s]->CalculateStatistics(smallOpts, svsigstrm, LookUpTable);
			}
			alignments.back().SetFromSegAlignment(smallOpts);
		}
		alignmentsOrder.Update(&alignments);
	}	
	//AlignmentsOrder alignmentsOrder(&alignments);
	SimpleMapQV(alignmentsOrder, read);

	// timing.Tick("Local-SDP");
	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
		for (int a = 0; a < (int) min(alignmentsOrder.size(), opts.PrintNumAln); a++){
			for (int s = 0; s < alignmentsOrder[a].SegAlignment.size(); s++) {
				alignmentsOrder[a].SegAlignment[s]->order = s;
				alignmentsOrder[a].SegAlignment[s]->wholegenomeLen = genome.header.pos[alignmentsOrder[a].SegAlignment[s]->chromIndex];
				
				if (opts.printFormat == "b") {
					alignmentsOrder[a].SegAlignment[s]->PrintBed(*output);
				}
				else if (opts.printFormat == "s") {
					alignmentsOrder[a].SegAlignment[s]->PrintSAM(*output, opts, alignmentsOrder[a].SegAlignment, s, read.passthrough);
				}
				else if (opts.printFormat == "a") {
					alignmentsOrder[a].SegAlignment[s]->PrintPairwise(*output);
				}
				else if (opts.printFormat == "p" or opts.printFormat == "pc") {
					alignmentsOrder[a].SegAlignment[s]->PrintPAF(*output, opts.printFormat == "pc");
				}
			}
		}
	}
	else {
		if (opts.printFormat == "s") {
			Alignment unaligned;
			unaligned.read=read.seq;
			unaligned.readLen=read.length;
			unaligned.SimplePrintSAM(*output, opts, read.passthrough);
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
