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
//#include "NaiveDP.h"
#include "GlobalChain.h"
#include "TupleOps.h"
#include "SparseDP.h"
#include "Merge.h"

using namespace std;



void SwapStrand(Read &read, Options &opts, GenomePairs &matches) {
	for (int m=0; m < matches.size(); m++) {
		matches[m].first.pos = read.length - (matches[m].first.pos + opts.globalK);
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
	vector<bool> strand(allMatches.size());
	int nForward=0;
	for (int i=0; i < allMatches.size(); i++) {
		int readPos = allMatches[i].first.pos;
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


void traceback(vector<int> & clusterchain, int & i, vector<int> & clusters_predecessor, vector<bool> & used, int & m) {

	clusterchain.push_back(i);
	used[i] = 1;
	cerr << "push back " << i << endl;

	if (clusters_predecessor[i] != -1) {
		if (used[clusters_predecessor[i]] != 1) {
			traceback(clusterchain, clusters_predecessor[i], clusters_predecessor, used, m);			
		}
		else {
			for (int lr = 0; lr < clusterchain.size(); ++lr) {
				used[clusterchain[lr]] = 0;
			}			
		}
	}
	else ++m;
}



void MapRead(const vector<float> & LookUpTable, Read &read, Genome &genome, vector<GenomeTuple> &genomemm, LocalIndex &glIndex, Options &opts, ostream *output, pthread_mutex_t *semaphore=NULL) {
	
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

	//DiagonalSort<GenomeTuple>(allMatches); // sort fragments in allMatches by forward diagonal, then by first.pos(read)

	// TODO(Jinwen): delete this after debug
	if (opts.dotPlot) {
		ofstream clust("all-matches.dots");
		for (int m=0; m < allMatches.size(); m++) {
			clust << allMatches[m].first.pos << "\t" << allMatches[m].second.pos << "\t" << opts.globalK << "\t0\t0"<<endl;
		}
		clust.close();
	}


	SeparateMatchesByStrand(read, genome, opts.globalK, allMatches, forMatches, revMatches);
/*
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
*/

	DiagonalSort<GenomeTuple>(forMatches); //sort fragments in forMatches by forward diagonal, then by first.pos(read)
	CleanOffDiagonal(forMatches, opts);

	vector<Cluster> clusters;
	int forwardStrand=0;
	CartesianSort(forMatches); // sort fragments in forMatches by first.pos, then by second.pos
	StoreDiagonalClusters(forMatches, clusters, opts, false, forwardStrand);

	AntiDiagonalSort<GenomeTuple>(revMatches, genome.GetSize());
	// TODO(Jingwen): later improve this part
	//SwapStrand(read, opts, revMatches);
	CleanOffDiagonal(revMatches, opts, 1);
	//SwapStrandBack(read, opts, revMatches);



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





	
	vector<Cluster> revClusters;	
	int reverseStrand=1;
	CartesianTargetSort(revMatches);
	StoreDiagonalClusters(revMatches, revClusters, opts, false, reverseStrand);
	//
	// Add pointers to seq that make code more readable.
	//
	char *readRC;

	CreateRC(read.seq, read.length, readRC);
	char *strands[2] = { read.seq, readRC };

	clusters.insert(clusters.end(), revClusters.begin(), revClusters.end());
	

	if (opts.dotPlot) {
		ofstream clust("clusters-pre-merge.tab");
		for (int c =0; c < clusters.size(); c++) {

			for (int m=0; m < clusters[c].matches.size(); m++) {
				clust << clusters[c].matches[m].first.pos << "\t" << clusters[c].matches[m].second.pos << "\t" << opts.globalK << "\t" << c << "\t" << clusters[c].strand << endl;
			}

			//clust << clusters[c].qStart << "\t" << clusters[c].tStart << "\t" << clusters[c].qEnd << "\t" << clusters[c].tEnd << "\t" << c  << "\t" << clusters[c].strand << endl;
		}
		clust.close();
	}

	ClusterOrder clusterOrder(&clusters);  // clusterOrder is sorted first by strand, then by tStart, then by qStart
	//MergeAdjacentClusters(clusterOrder, genome, opts);


	// Merge clusters
	// use naive way to find top chains
	vector<float> clusters_value(clusters.size(), 0);
	vector<int> clusters_predecessor(clusters.size(), -1);
	clusters_value[0] = (float)(max(clusterOrder[0].qEnd - clusterOrder[0].qStart, clusterOrder[0].tEnd - clusterOrder[0].tStart));
	vector<int> repetitivenum(clusters.size(), 0);

	if (clusters.size() > 1) {
		for (int c = 1; c < clusters.size(); ++c) {

			clusters_value[c] = (float)(max(clusterOrder[c].qEnd - clusterOrder[c].qStart, clusterOrder[c].tEnd - clusterOrder[c].tStart));

			for (int s = 0; s <= c - 1; ++s) {
				if(clusterOrder[c].Overlaps(clusterOrder[s], 0.8) and clusterOrder[c].OverlapsOnRead(read.length, 0.7) 
					and clusterOrder[s].OverlapsOnRead(read.length, 0.7)) {++repetitivenum[c];}

				int gap = abs((int)(clusterOrder[c].qStart - clusterOrder[c].tStart) - (int)(clusterOrder[s].qEnd - clusterOrder[s].tEnd));
				if ( (clusterOrder[s].qEnd - (clusterOrder[s].qEnd - clusterOrder[s].qStart)/4) < clusterOrder[c].qStart 
						and  (clusterOrder[s].tEnd - (clusterOrder[s].tEnd - clusterOrder[s].tStart)/4) < clusterOrder[c].tStart 
						and !clusterOrder[c].Encompasses(clusterOrder[s], 0.4) //and clusterOrder[c].strand == clusterOrder[s].strand
						and gap < opts.maxGap) {

					float objective;
					if (gap < 501)  {
						//objective = clusters_value[c] + clusters_value[s] - log((float)gap) - 0.5*((float)gap);
						objective = clusters_value[c] + clusters_value[s] - 0.5*((float)gap);

					}
					else {
						//objective = clusters_value[c] + clusters_value[s] - LookUpTable[(int)floor((gap-501)/5)] - 0.5*((float)gap);
						objective = clusters_value[c] + clusters_value[s] - 0.5*((float)gap);
					}				
					if (objective >= clusters_value[c]) {
						clusters_value[c] = objective;
						clusters_predecessor[c] = s;
					}			
				}
			} 

		}	
	}


	int repetitive = *max_element(repetitivenum.begin(), repetitivenum.end()); //repetitive+1 shows the number of chains we want to trace back

	//cerr << "the number of chains we want to trace back: " << repetitive+1 << endl;
	vector<vector<int>> clusterchain(repetitive+1);
	Clusters_valueOrder clusters_valueOrder(&clusters_value); // sort clusters_value in descending order
	vector<bool> used(clusters.size(), 0);

	// traceback 
	for (int r= 0, m = 0; m <= repetitive; ++r) {
		//cerr << "r: " << r << "  m: " << m << endl;
		clusterchain[m].clear();
		if (used[clusters_valueOrder.index[r]] != 1) {
			//cerr << " used[" << clusters_valueOrder.index[r] << "] != 1" << endl;
			traceback(clusterchain[m], clusters_valueOrder.index[r], clusters_predecessor, used, m);
		}
	}
	used.clear();

/*
	float max_value = 0;
	int max_ind = -1;

	for (int c = 0; c < clusters.size(); ++c){
		cerr << "c: " << c << "  clusters_value[c]: " << clusters_value[c] << endl;
		if (max_value < clusters_value[c]) {
			max_value = clusters_value[c];
			max_ind = c;
		}
	}

	cerr << "max_ind: " << max_ind << "  max_value: " << max_value << endl;

	traceback(clusterchain, max_ind, clusters_predecessor);
*/

	for (int r = 0; r < clusterchain.size(); ++r) {
		for (int c = 0; c < clusterchain[r].size(); ++c) {
			clusterchain[r][c] = clusterOrder.index[clusterchain[r][c]];			
		}	
	}

	//cerr << "clusterchain.size(): " << clusterchain.size() << endl;
	if (opts.dotPlot) {
		for (int c = 0; c < clusterchain.size(); c++) {
			stringstream outNameStrm;
			outNameStrm << "clusters-sdp." << c << ".dots";
			ofstream baseDots(outNameStrm.str().c_str());
			// clusterchain stores indices which refer to elements in clusters[clusterchain[c]].matches
			for (int s = 0; s < clusterchain[c].size(); ++s) {
				for (int m = 0; m < clusters[clusterchain[c][s]].matches.size(); ++m) {
					baseDots << clusters[clusterchain[c][s]].matches[m].first.pos << "\t" 
							 << clusters[clusterchain[c][s]].matches[m].second.pos << "\t" 
							 << opts.globalK <<"\t" << c << "\t" << clusters[clusterchain[c][s]].strand << endl;					
				}				
			}
			baseDots.close();
		}
	}

	// delete other clusters that are not in top chains 
	vector<bool> ind(clusters.size(), 0);
	for (int c = 0; c < clusterchain.size(); ++c) {
		ind[clusterchain[c][0]] = 1;
		if (clusterchain[c].size() > 1) {
			for (int s = 1; s < clusterchain[c].size(); ++s) {
				clusters[clusterchain[c][0]].matches.insert(clusters[clusterchain[c][0]].matches.end(), clusters[clusterchain[c][s]].matches.begin(), clusters[clusterchain[c][s]].matches.end());
				clusters[clusterchain[c][0]].qEnd = clusters[clusterchain[c][s]].qEnd;
				clusters[clusterchain[c][0]].tEnd = clusters[clusterchain[c][s]].tEnd;
				clusters[clusterchain[c][s]].tEnd=0;
				clusters[clusterchain[c][s]].tStart=0;
				clusters[clusterchain[c][s]].matches.clear();				
			}			
		}
	}		

	// delete clusters with small size
	for (int c= 0; c < clusterchain.size(); ++c) {
		if (((float)(clusters[clusterchain[c][0]].qEnd - clusters[clusterchain[c][0]].qStart))/read.length < 0.5) ind[clusterchain[c][0]] = 0;
	}

	int cl=0;
	int cn;
	for(cl=0, cn=0; cn < clusters.size(); cn++) {
		if (ind[cn] == 1) {
			clusters[cl] = clusters[cn];
			cl++;
		}
	}
	clusters.resize(cl);
	clusterchain.clear();


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
	
	sort(clusters.begin(), clusters.end(), OrderClusterBySize()); // clusters are sorted in descending order
	
	if (opts.dotPlot) {
		ofstream matchfile("long_matches.tab");
		for (int m =0; m < matches.size(); m++) {
			matchfile << matches[m].first.pos << "\t" << matches[m].second.pos << "\t" << opts.globalK << "\t0\t0" << endl;
		}
		matchfile.close();
		ofstream clust("clusters.tab");
		for (int c =0; c < clusters.size(); c++) {
	
			for (int m=0; m < clusters[c].matches.size(); m++) {
				clust << clusters[c].matches[m].first.pos << "\t" << clusters[c].matches[m].second.pos << "\t" << opts.globalK << "\t" << c << "\t" << clusters[c].strand << endl;
			}	
			//clust << clusters[c].qStart << "\t" << clusters[c].tStart << "\t" << clusters[c].qEnd << "\t" << clusters[c].tEnd << "\t" << c << "\t" << clusters[c].strand << endl;
		}
		clust.close();	
	}

	//RemoveOverlappingClusters(clusters, opts); //TODO(Jingwen): check whether need to keep this


	// TODO(Jingwen): only for debug && delete this later 
	if (opts.dotPlot) {
		ofstream clust("clusters-after-remove-overlapping.tab");
		for (int c =0; c < clusters.size(); c++) {

			for (int m=0; m < clusters[c].matches.size(); m++) {
				clust << clusters[c].matches[m].first.pos << "\t" << clusters[c].matches[m].second.pos << "\t" << opts.globalK << "\t" << c << "\t" << clusters[c].strand << endl;
			}
			//clust << clusters[c].qStart << "\t" << clusters[c].tStart << "\t" << clusters[c].qEnd << "\t" << clusters[c].tEnd << "\t" << c  << "\t" << clusters[c].strand << endl;
		}
		clust.close();
	}

	
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
				baseDots << clusters[c].matches[m].first.pos << "\t" 
								 << clusters[c].matches[m].second.pos << "\t" 
								 << smallOpts.globalK << "\t" << c << "\t0" << endl;
			}
			baseDots.close();
		}

		//
		// Get the boundaries of the cluster in both sequences.
		//
		CartesianTargetSort<GenomeTuple>(clusters[c].matches.begin(), clusters[c].matches.end()); // sorted by second.pos and then first.pos
		int nMatch = clusters[c].matches.size();
		GenomePos tPos=clusters[c].matches[0].second.pos;
		int firstChromIndex = genome.header.Find(tPos);
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
		// TODO(Jingwen): clusters[c].matches is cartesianTargetSorted(first by second.pos and then by first.pos)
		// but in this way, readClusterStart and readClusterEnd are not quite accurate. ????????				
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
			GenomePos genomeLocalIndexStart = glIndex.seqOffsets[lsi]  - chromOffset;
			GenomePos genomeLocalIndexEnd   = glIndex.seqOffsets[lsi+1] - 1 - chromOffset;

			int matchStart = CartesianTargetLowerBound<GenomeTuple>(clusters[c].matches.begin(), clusters[c].matches.end(), genomeLocalIndexStart);

			int matchEnd   = CartesianTargetUpperBound<GenomeTuple>(clusters[c].matches.begin(), clusters[c].matches.end(), genomeLocalIndexEnd);

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
			int queryIndexEnd = readIndex->LookupIndex(min(readEnd,   (GenomePos) read.length-1));
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

				DiagonalSort<LocalTuple>(smallMatches);
				CleanOffDiagonal(smallMatches, smallOpts);// TODO(Jingwen): resume this after debug			

				AppendValues<LocalPairs>(refinedClusters[c].matches, smallMatches.begin(), smallMatches.end(), readSegmentStart, genomeLocalIndexStart );
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
				baseDots << refinedClusters[c].matches[m].first.pos << "\t" << refinedClusters[c].matches[m].second.pos << "\t" << smallOpts.globalK << "\t" << c << endl;
			}
			baseDots.close();
		}
			
	}

	//
	// Remove clusters under some minimal number of anchors. By default this is one. 
	//
	//RemoveEmptyClusters(refinedClusters);
	


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
		//int k=glIndex.k;
		if (refinedClusters[r].matches.size() > read.length or 
				opts.refineLevel & REF_LOC == 0) {
			refinedClusters[r].matches= clusters[refinedClusters[r].coarse].matches;
			//k=opts.globalK;
			continue;
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



		vector<Cluster> vt; // vt stores the split result of merged anchors
		if (opts.mergeClusters) {
			vt.clear();

			//MergeClusters(smallOpts, refinedClusters, vt, r);
			mergeClusters (smallOpts, refinedClusters[r].matches, vt, r, baseName);

			//cerr << "vt.size(): " << vt.size() << endl;
			if (opts.dotPlot) {
				stringstream outNameStrm;
				outNameStrm << baseName + "." << r << ".merged.dots";
				ofstream baseDots(outNameStrm.str().c_str());
				for (int m=0; m < vt.size(); m++) {
					baseDots << vt[m].qStart << "\t" 
							 << vt[m].tStart << "\t" 
							 << vt[m].qEnd << "\t" 
							 << vt[m].tEnd << "\t" 
							 << r << endl;
				}
				baseDots.close();
			}


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

		}


		// Perform sparse chaining, uses time O(n log n).
		//
		// Merge anchors that are not overlapping

		//
		// The input are a buntch of fragments which are either stored in vector<Cluster> vt or refinedClusters[r].matches
		// 
		vector<unsigned int> chain; // chain contains the index of the fragments involved in the result of SparseDP.
		if (opts.SparseDP) {
			chain.clear();
			if (opts.mergeClusters and vt.size() < 30000 and vt.size() > 0) {
				// TODO(Jingwen): recover this after fix sdp
				//SparseDP(vt, chain, smallOpts, LookUpTable);
			}
			else if (refinedClusters[r].matches.size() < 30000) {
				SparseDP(refinedClusters[r].matches, chain, smallOpts, LookUpTable);
			}
		}
		else { 
			chain.clear();
			// TODO(Jingwen): add linear sdp here?
		}
		if (opts.dotPlot) {
			stringstream outNameStrm;
			outNameStrm << baseName + "." << r << ".first-sdp.dots";
			ofstream baseDots(outNameStrm.str().c_str());
			for (int c = 0; c < chain.size(); c++) {
				if (opts.mergeClusters) {
					// chain stores indices which refer to elments in vt
					baseDots << vt[chain[c]].qStart << "\t" 
							 << vt[chain[c]].tStart << "\t" 
							 << vt[chain[c]].qEnd << "\t" 
							 << vt[chain[c]].tEnd << "\t"
							 << r << endl;							
				}
				else {
					// chain stores indices which refer to elements in refinedClusters[r].matches
					baseDots << refinedClusters[r].matches[chain[c]].first.pos << "\t" 
							 << refinedClusters[r].matches[chain[c]].second.pos << "\t" 
							 << smallOpts.globalK <<"\t" << r << endl;						
				}		
			}
			baseDots.close();
		}

		GenomePairs tupChain;

		if (opts.mergeClusters and vt.size() > 0) {
			//
			// Add small anchors to tupChain. (Use greedy algorithm to make small anchors not overlap with each other)
			//
			for (int ch=0; ch < chain.size(); ch++) { 

				//cerr << "ch: "<< ch<< endl;
				//cerr << "chain[" << ch <<"]   " << chain[ch] << endl; 

				unsigned int fprev = vt[chain[ch]].start;
				unsigned int fcur = vt[chain[ch]].start;
				tupChain.push_back(GenomePair(GenomeTuple(0, refinedClusters[r].matches[fcur].first.pos), 
																			GenomeTuple(0, refinedClusters[r].matches[fcur].second.pos)));	
				
				//cerr << "tupChain[" << tupChain.size() << "]     push_back: " << "first: " << refinedClusters[r].matches[fcur].first.pos<<"  " << refinedClusters[r].matches[fcur].first.pos + smallOpts.globalK 
				//<< "   second: " << refinedClusters[r].matches[fcur].second.pos << " "<< refinedClusters[r].matches[fcur].second.pos + smallOpts.globalK - 1 <<endl;
				

				++fcur;					
				while (fcur < vt[chain[ch]].end) {
					while (fcur < vt[chain[ch]].end && (refinedClusters[r].matches[fcur].first.pos < refinedClusters[r].matches[fprev].first.pos + smallOpts.globalK ||
												refinedClusters[r].matches[fcur].second.pos < refinedClusters[r].matches[fprev].second.pos + smallOpts.globalK)) {
						++fcur;
					}
					if (fcur != vt[chain[ch]].end) {
						tupChain.push_back(GenomePair(GenomeTuple(0, refinedClusters[r].matches[fcur].first.pos), 
																GenomeTuple(0, refinedClusters[r].matches[fcur].second.pos)));	
				
						//cerr << "tupChain[" << tupChain.size() << "]     push_back: " << "first: " << refinedClusters[r].matches[fcur].first.pos<<"  " << refinedClusters[r].matches[fcur].first.pos + smallOpts.globalK 
						//<< "   second: " << refinedClusters[r].matches[fcur].second.pos << " "<< refinedClusters[r].matches[fcur].second.pos + smallOpts.globalK - 1 <<endl;
				
						assert(refinedClusters[r].matches[fcur].first.pos >= refinedClusters[r].matches[fprev].first.pos + smallOpts.globalK );
						assert(refinedClusters[r].matches[fcur].second.pos >= refinedClusters[r].matches[fprev].second.pos + smallOpts.globalK );								
					}
					fprev = fcur;
					++fcur;
				}
			}
			vt.clear();
		}
		else {
			// chain stores indices which refer to elements in refinedClusters[r].matches
			for (int ch=0; ch < chain.size(); ch++) {
				assert(refinedClusters[r].matches[chain[ch]].first.pos <= read.length); // TODO(Jingwen): delete this after debug
				tupChain.push_back(GenomePair(GenomeTuple(0, refinedClusters[r].matches[chain[ch]].first.pos), 
														GenomeTuple(0, refinedClusters[r].matches[chain[ch]].second.pos)));
			}
		}

		//(TODO)Jingwen: For Debug(remove this later)
		if (opts.dotPlot and opts.mergeClusters) {
			stringstream outNameStrm;
			outNameStrm << baseName + "." << r << ".first-sdp-clean.dots";
			ofstream baseDots(outNameStrm.str().c_str());
			for (int c = 0; c < tupChain.size(); c++) {
				// chain stores indices which refer to elments in vt
				baseDots << tupChain[c].first.pos << "\t" 
						 << tupChain[c].second.pos << "\t" 
						 << tupChain[c].first.pos + smallOpts.globalK << "\t" 
						 << tupChain[c].second.pos + smallOpts.globalK << "\t"
						 << r << endl;							
			}
			baseDots.close();
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

		// remove fragments which are in the middle of an insertion and a deletion OR a deletion and an insertion. 
		//StoreDiagonalClusters(tupChain, chainClust, diagOpts, true); ///Jingwen adds this here, otherwise chainClust is empty
		//RemovePairedIndels(tupChain, chainClust, smallOpts);  


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
		GenomePos chainGenomeStart = tupChain[0].second.pos;
		GenomePos chainGenomeEnd   = tupChain[tupChain.size()-1].second.pos + smallOpts.globalK;

		GenomePos chainReadStart = tupChain[0].first.pos;
		GenomePos chainReadEnd   = tupChain[tupChain.size()-1].first.pos + smallOpts.globalK;


		//GenomePos globalStart = refinedClusters[r].tStart;
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
		vector<GenomePairs> refinedChains(tupChain.size()-1); // Note: refinedChains[i] stores the gap-fragments which locate after chain[i]
		vector<int> refinedChainsLength(tupChain.size()-1, -1); // refinedChainsLength[i] stores the tinyOpts.globalK of the gap-fragments which locate after chain[i]


		vector<int> scoreMat;
		vector<Arrow> pathMat;

		int chainLength = tupChain.size();
		for (int c = 0; chainLength > 0 and c < tupChain.size()-1; c++) {
			GenomePos curGenomeEnd = tupChain[c].second.pos + smallOpts.globalK;
			GenomePos nextGenomeStart = tupChain[c+1].second.pos;

			GenomePos curReadEnd = tupChain[c].first.pos + smallOpts.globalK;
			GenomePos nextReadStart = tupChain[c+1].first.pos;

			assert(nextReadStart >= curReadEnd);
			GenomePos subreadLength = nextReadStart-curReadEnd;
			assert(nextGenomeStart >= curGenomeEnd); // TODO(Jingwen): change this to include inversion
			GenomePos subgenomeLength = nextGenomeStart-curGenomeEnd;

			if (nextReadStart > curReadEnd and nextGenomeStart > curGenomeEnd) {

				if (subreadLength > 50 and 
						subgenomeLength > 50 and
						opts.refineLevel & REF_DYN ) {

					// TODO(Jingwen): should we only consider minLen? because min(subreadLength, subgenomeLength) is the length of possible matches
					//GenomePos maxLen = max(subreadLength, subgenomeLength);
					GenomePos maxLen = min(subreadLength, subgenomeLength);
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
					StoreMinimizers<GenomeTuple, Tuple>( genome.seqs[chromIndex] + curGenomeEnd, subgenomeLength, tinyOpts.globalK, 1, gapGenomeTup, false);

					sort(gapGenomeTup.begin(), gapGenomeTup.end());
					StoreMinimizers<GenomeTuple, Tuple>( strands[refinedClusters[r].strand] + curReadEnd, subreadLength, tinyOpts.globalK, 1, gapReadTup, false);
					sort(gapReadTup.begin(), gapReadTup.end());
					CompareLists(gapReadTup.begin(), gapReadTup.end(), gapGenomeTup.begin(), gapGenomeTup.end(), gapPairs, gapOpts);
					//
					// Remove egregious off diagonal seeds
					//
					DiagonalSort<GenomeTuple>(gapPairs); // sort gapPairs by forward diagonals and then by first.pos
					
					int tinyDiagStart = curReadEnd - curGenomeEnd;
					int tinyDiagEnd =  nextReadStart - nextGenomeStart;
					int diagDiff = abs((int) tinyDiagStart - (int) tinyDiagEnd);

					CleanOffDiagonal(gapPairs, tinyOpts, 0, 0, diagDiff);
					//CartesianTargetSort<GenomeTuple>(gapPairs);

					for(int rm=0; rm < gapPairs.size(); rm++) {
						gapPairs[rm].first.pos  += curReadEnd;
						gapPairs[rm].second.pos += curGenomeEnd;
						assert(gapPairs[rm].first.pos < read.length); // TODO(Jingwen): delete this after debug
					}
					//int gp=0;
					//int nSaved =0;

					//int gpStart;
					//int clusterIndex=0;

					// gapPairs stores all the fragments in the gap. the length of the fragment here is tinyOpts.globalK
					vector<unsigned int> gapChain; // gapChain stores the index of fragments that are involved in the result of sdp
					if (gapPairs.size() > 0) {
						if (opts.SparseDP) {
							gapChain.clear();
							if (gapPairs.size() < 20000) {
								SparseDP(gapPairs, gapChain, tinyOpts, LookUpTable); 
								//cerr << "start 2nd sdp!" << endl;

								if (opts.dotPlot) {
									stringstream outNameStrm;
									outNameStrm << baseName + "." << r << ".second-sdp.dots";
									ofstream baseDots;
									baseDots.open(outNameStrm.str().c_str(), std::ios::app);
									for (int c = 0; c < gapChain.size(); c++) {
										// chain stores indices which refer to elements in refinedClusters[r].matches
										baseDots << gapPairs[gapChain[c]].first.pos << "\t" 
												 << gapPairs[gapChain[c]].second.pos << "\t" 
												 << gapOpts.globalK  << "\t"
												 << r << endl;						
									}
									baseDots.close();
								}

							}
						}
						/*
						else {
							//
							// Test code for global chain. Modify to make it work efficiently if correct.
							//
							vector<Fragment> fragments;
							vector<Endpoint> endpoints;
							vector<int> op;
							seqan::Iterator<seqan::SeedSet<seqan::Seed<seqan::Simple>, seqan::Unordered> >::Type it;
							for ( it = seqan::begin(gapSeedSet); it != seqan::end(gapSeedSet); ++it) {
								fragments.push_back(Fragment(beginPositionH(*it),beginPositionV(*it), endPositionH(*it),endPositionV(*it),endPositionH(*it)-beginPositionH(*it),0));
							}

							GlobalChain(fragments, op, endpoints);
							int i;
							for (int oi=0; oi<op.size(); oi++) {
								IndSeed s(fragments[op[oi]].xl, 
													fragments[op[oi]].yl, 
													fragments[op[oi]].xh, 
													fragments[op[oi]].yh);
								seqan::append(gapChain,s);
							}
						}
						*/
					}

					RemovePairedIndels(curReadEnd, curGenomeEnd, nextReadStart, nextGenomeStart, gapChain, gapPairs, tinyOpts);
					for (unsigned int u = 0; u < gapChain.size(); ++u) {
						refinedChains[c].push_back(gapPairs[gapChain[u]]);
						assert(refinedChains[c].back().first.pos <= read.length);

						if (refinedChains[c].size() > 1) { // TODO(Jingwen): delete the debug code
							int last = refinedChains[c].size();
							assert(refinedChains[c][last -2].first.pos + tinyOpts.globalK <= refinedChains[c][last -1].first.pos);
							assert(refinedChains[c][last -2].second.pos + tinyOpts.globalK <= refinedChains[c][last -1].second.pos);						
						}
					}
					refinedChainsLength[c] = tinyOpts.globalK;	
				}
			}	
		}

		//
		// Refine and store the alignment
		//
		//
		// The alignment is on a substring that starts at the beginning of the first chain.
		//
		GenomePos alnReadPos = tupChain[0].first.pos;
		GenomePos alnRefPos  = tupChain[0].second.pos;
		Alignment *alignment = new Alignment(strands[refinedClusters[r].strand], read.seq, read.length, read.name, refinedClusters[r].strand, genome.seqs[chromIndex],  
											genome.lengths[chromIndex], genome.header.names[chromIndex], chromIndex);

		alignments.push_back(alignment);
		for (int c = 0; chainLength> 0 and  c < chainLength-1; c++) {
			//
			// Chain is with respect to full sequence
			//
			GenomePos curGenomeEnd     = tupChain[c].second.pos + smallOpts.globalK;
			GenomePos nextGenomeStart  = tupChain[c+1].second.pos;

			GenomePos curReadEnd       = tupChain[c].first.pos + smallOpts.globalK;
			GenomePos nextReadStart    = tupChain[c+1].first.pos;
			int curRefinedReadEnd      = curReadEnd;
			int curRefinedGenomeEnd    = curGenomeEnd;
			int nextRefinedReadStart   = nextReadStart;
			int nextRefinedGenomeStart = nextGenomeStart;
			if (opts.dotPlot) {
				dotFile << tupChain[c].first.pos << "\t" 
								<< tupChain[c].second.pos << "\t" 
								<< tupChain[c].first.pos + smallOpts.globalK << "\t"
								<< tupChain[c].second.pos + smallOpts.globalK << "\t" << r << endl;
			}

			alignment->blocks.push_back(Block(tupChain[c].first.pos, tupChain[c].second.pos, smallOpts.globalK)); 
			if (alignment->blocks.size() > 1) {
				int last=alignment->blocks.size();
				assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].qPos);
				assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].tPos);
			}

			//string curAnchor = string(genome.seqs[chromIndex], seqan::beginPositionV(chain[c]), glIndex.k );
			for (int cs = 0; cs < refinedChains[c].size(); cs++) {
				//
				// Refined anchors are with respect to the chained sequence
				nextRefinedReadStart = refinedChains[c][cs].first.pos;
				nextRefinedGenomeStart = refinedChains[c][cs].second.pos;
				
				if (opts.dotPlot) {
					dotFile << refinedChains[c][cs].first.pos << "\t" 
									<< refinedChains[c][cs].second.pos << "\t" 
									<< refinedChains[c][cs].first.pos + refinedChainsLength[c] << "\t"
									<< refinedChains[c][cs].second.pos + refinedChainsLength[c] << "\t" << r << endl;
				}

				// find small matches between fragments in gapChain
				int m, rg, gg;
				SetMatchAndGaps(curRefinedReadEnd, nextRefinedReadStart, curRefinedGenomeEnd, nextRefinedGenomeStart, m, rg, gg);
				if (m > 0) {
					Alignment betweenAnchorAlignment;
					if (opts.refineLevel & REF_DP) {						
						RefineSubstrings(strands[refinedClusters[r].strand], curRefinedReadEnd, nextRefinedReadStart, genome.seqs[chromIndex], 
														 curRefinedGenomeEnd, nextRefinedGenomeStart, scoreMat, pathMat, betweenAnchorAlignment, opts);
						alignment->blocks.insert(alignment->blocks.end(), betweenAnchorAlignment.blocks.begin(), betweenAnchorAlignment.blocks.end());
						int b;
						for (b=1; b < betweenAnchorAlignment.blocks.size(); b++) {
							assert(betweenAnchorAlignment.blocks[b-1].qPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].qPos);
							assert(betweenAnchorAlignment.blocks[b-1].tPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].tPos);						
						}
						betweenAnchorAlignment.blocks.clear();
					}
				}

				curRefinedReadEnd = refinedChains[c][cs].first.pos + refinedChainsLength[c];
				curRefinedGenomeEnd = refinedChains[c][cs].second.pos + refinedChainsLength[c];
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
					RefineSubstrings(strands[refinedClusters[r].strand], curRefinedReadEnd, nextReadStart, genome.seqs[chromIndex], 
													 curRefinedGenomeEnd, nextGenomeStart, scoreMat, pathMat, aln, opts);
					alignment->blocks.insert(alignment->blocks.end(), aln.blocks.begin(), aln.blocks.end());
					aln.blocks.clear();			
				}		
			}
		}
		//		refc.close();
		//			rdc.close();
		alignment->blocks.push_back(Block(tupChain[chainLength-1].first.pos, tupChain[chainLength-1].second.pos, smallOpts.globalK));

		int nm=0;
		for(int b=0; b < alignment->blocks.size(); b++) {
			nm+= alignment->blocks[b].length;
		}
		alignment->nblocks = tupChain.size();
		if (opts.dotPlot) {
			dotFile.close();
		}

		//(TODO)Jingwen: For Debug(remove this later)
		if (opts.dotPlot) {
			stringstream outNameStrm;
			outNameStrm << baseName + "." << r << ".alignment.dots";
			ofstream baseDots(outNameStrm.str().c_str());
			for (int c = 0; c < alignment->blocks.size(); c++) {
				baseDots << alignment->blocks[c].qPos << "\t" 
						 << alignment->blocks[c].tPos << "\t" 
						 << alignment->blocks[c].qPos + alignment->blocks[c].length << "\t" 
						 << alignment->blocks[c].tPos + alignment->blocks[c].length << "\t"
						 << r << endl;							
			}
			baseDots.close();
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


		// for debug. TODO(Jingwen): delete this later
		if (opts.dotPlot) {
			stringstream outNameStrm;
			outNameStrm << baseName << ".finalalignment.dots";
			ofstream baseDots(outNameStrm.str().c_str());
			for (int c = 0; c < alignments[a]->blocks.size(); c++) {
				baseDots << alignments[a]->blocks[c].qPos << "\t" 
						 << alignments[a]->blocks[c].tPos << "\t" 
						 << alignments[a]->blocks[c].qPos + alignments[a]->blocks[c].length << "\t" 
						 << alignments[a]->blocks[c].tPos + alignments[a]->blocks[c].length << "\t"
						 << a << endl;							
			}
			baseDots.close();
		}	

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

	/*
	// get the time for the program
	clock_t end = std::clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cerr << "Time: " << elapsed_secs << endl;
	*/
	
}


#endif
