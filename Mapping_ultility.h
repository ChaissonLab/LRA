#ifndef MAP_ULTILITY_H_
#define MAP_ULTILITY_H_
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

//
// This function switches index in splitclusters back 
//
int
switchindex (vector<Cluster> & splitclusters, vector<Primary_chain> & Primary_chains, vector<Cluster> & clusters, Genome &genome, Read &read) {
	if (read.unaligned) return 0;
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
			map<int, int> appeartimes;
			map<int, int> start_pos;
			map<int, int> end_pos;
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
	//
	// Remove small clusters that are total covered by larger cluster on q coordinates
	//
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			vector<bool> cremove(Primary_chains[p].chains[h].ch.size(), 0);
			for (int c = 1; c < Primary_chains[p].chains[h].ch.size(); c++) {
				int cr = Primary_chains[p].chains[h].ch[c];
				int cp = Primary_chains[p].chains[h].ch[c - 1];
				if (cremove[c - 1] == 0 and clusters[cr].qStart >= clusters[cp].qStart and clusters[cr].qEnd <= clusters[cp].qEnd){
					cremove[c] = 1;
				}
			}
			int sc = 0;
			for (int c = 0; c < Primary_chains[p].chains[h].ch.size(); c++) {
				if (cremove[c] == 0) {
					Primary_chains[p].chains[h].ch[sc] = Primary_chains[p].chains[h].ch[c];
					if (sc >= 1) Primary_chains[p].chains[h].link[sc - 1] = Primary_chains[p].chains[h].link[c - 1];
					sc++;
				}
			}
			Primary_chains[p].chains[h].ch.resize(sc);
			Primary_chains[p].chains[h].link.resize(sc - 1);
		}
	}
	return 0;
}

//
// This function splits the chain if Clusters on the chain are mapped to different chromosomes or different locations (quite far, default: 100000) on the same chromosome;
// Also split the chain when two forward/reverse clusters are chained in reverse/forward direction.
void
SPLITChain(Read &read, vector<Cluster_SameDiag *> &ExtendClusters, vector<SplitChain> & splitchains, vector<bool> & link, const Options & opts) {
	int im = 0;
	vector<int> onec; 
	vector<bool> lk;
	onec.push_back(im);
	int cur = 0, prev = 0;

	while (im < ExtendClusters.size() - 1) {
		cur = im + 1; prev = im;
		if (ExtendClusters[cur]->tStart > ExtendClusters[prev]->tEnd + opts.splitdist // too far
			or ExtendClusters[cur]->tEnd + opts.splitdist < ExtendClusters[prev]->tStart
			or (link[im] == 1 and ExtendClusters[cur]->strand == 0 and ExtendClusters[prev]->strand == 0) // repetitive mapping
			or (link[im] == 0 and ExtendClusters[cur]->strand == 1 and ExtendClusters[prev]->strand == 1) // repetitive mapping
			or (ExtendClusters[cur]->strand == 0 and ExtendClusters[prev]->strand == 1) // inversion
			or (ExtendClusters[cur]->strand == 1 and ExtendClusters[prev]->strand == 0) // inversion
			or (ExtendClusters[prev]->OverlaprateOnGenome(ExtendClusters[cur]) >= 0.3)
			or  (ExtendClusters[prev]->OverlapOnGenome(ExtendClusters[cur]) >= 500 // If two clusters overlap exceeds 0.3, then it is a DUP
				and ExtendClusters[prev]->anchorfreq <= 1.05f and ExtendClusters[cur]->anchorfreq <= 1.05f)) {  // If two clusters overlap exceeds 100bp and they are both linear, then it is a DUP
			splitchains.push_back(SplitChain(onec, lk));
			onec.clear();
			lk.clear();
			onec.push_back(cur);
		}	
		else {
			onec.push_back(cur);
			lk.push_back(link[im]);
		}	
		im++;
	}
	if (!onec.empty()) splitchains.push_back(SplitChain(onec, lk));

	for (int m = 0; m < splitchains.size(); m++) {
		splitchains[m].QStart = ExtendClusters[splitchains[m][0]]->qStart;
		splitchains[m].QEnd = ExtendClusters[splitchains[m][0]]->qEnd;
		for (int n = 1; n < splitchains[m].size(); n++) {
			splitchains[m].QStart = min(splitchains[m].QStart, ExtendClusters[splitchains[m][n]]->qStart);
			splitchains[m].QEnd = max(splitchains[m].QEnd, ExtendClusters[splitchains[m][n]]->qEnd);
		}
	}
}

//
// This function splits the chain if Clusters on the chain are mapped to different chromosomes or different locations (quite far, default: 100000) on the same chromosome;
// Also split the chain when two forward/reverse clusters are chained in reverse/forward direction.
void
SPLITChain(Read &read, UltimateChain &chain, vector<SplitChain> &splitchains, vector<bool> &splitchains_link, 
				vector<pair<GenomePos, GenomePos>> &splitchains_qpos, const Options &opts) {
	int im = 0;
	vector<int> onec; 
	vector<bool> lk;
	onec.push_back(im);
	int cur = 0, prev = 0;

	while (im < chain.size() - 1) {
		cur = im + 1; prev = im;
		if (chain.tStart(cur) > chain.tEnd(prev) + opts.splitdist // too far
			or chain.tEnd(cur) + opts.splitdist < chain.tStart(prev)
			or (chain.link[im] == 1 and chain.strand(cur) == 0 and chain.strand(prev) == 0) // repetitive mapping and DUP
			or (chain.link[im] == 0 and chain.strand(cur)== 1 and chain.strand(prev) == 1) // repetitive mapping and DUP
			or (chain.strand(cur) == 0 and chain.strand(prev) == 1) // inversion
			or (chain.strand(cur) == 1 and chain.strand(prev) == 0)) { //inversion

			splitchains.push_back(SplitChain(onec, lk));
			splitchains_qpos.push_back(make_pair(chain.qStart(onec.back()), chain.qEnd(onec[0])));
			if ((chain.strand(cur) == 0 and chain.strand(prev) == 1) or (chain.strand(cur) == 1 and chain.strand(prev) == 0)) {splitchains_link.push_back(1);}
			else {splitchains_link.push_back(0);}
			onec.clear();
			lk.clear();
			onec.push_back(cur);
		}	
		else {
			onec.push_back(cur);
			lk.push_back(chain.link[im]);
		}	
		im++;
	}
	if (!onec.empty()) {
		splitchains.push_back(SplitChain(onec, lk));
		splitchains_qpos.push_back(make_pair(chain.qStart(onec.back()), chain.qStart(onec[0])));
	}
}

bool push_new(Genome &genome, vector<int> &onec, vector<bool> &lk, vector<SplitChain> &splitchains, vector<bool> &splitchains_link, UltimateChain &chain, int cur) {

	splitchains.push_back(SplitChain(onec, lk, &chain, chain.strand(onec[0])));
	splitchains.back().ClusterIndex.push_back(chain.ClusterNum(onec[0]));
	for (int c = 1; c < onec.size(); c++) {
		if (chain.ClusterNum(onec[c]) != splitchains.back().ClusterIndex.back()) {
			splitchains.back().ClusterIndex.push_back(chain.ClusterNum(onec[c]));
		}
	}
	splitchains.back().QStart = chain.qStart(onec.back()); splitchains.back().QEnd = chain.qEnd(onec[0]);
	if (chain.strand(onec[0]) == 0) {
		splitchains.back().TStart = chain.tStart(onec.back()); splitchains.back().TEnd = chain.tEnd(onec[0]);
	}
	else {
		splitchains.back().TStart = chain.tStart(onec[0]); splitchains.back().TEnd = chain.tEnd(onec.back());
	}

	if (splitchains.back().CHROMIndex(genome)) {
		splitchains.pop_back();
		onec.clear();
		lk.clear();
		onec.push_back(cur);
		return 0;
	}
	onec.clear();
	lk.clear();
	onec.push_back(cur);
	return 1;
}

void
SPLITChain(Genome &genome, Read &read, UltimateChain &chain, vector<SplitChain> &splitchains, vector<bool> &splitchains_link, const Options &opts) {
	int im = 0;
	vector<int> onec; 
	vector<bool> lk;
	onec.push_back(im);
	int cur = 0, prev = 0;

	while (im < chain.size() - 1) {
		cur = im + 1; prev = im;
		int dist = 0;
		assert(chain.qStart(prev) >= chain.qEnd(cur));
		int qdist = chain.qStart(prev) - chain.qEnd(cur);
		int tdist = (chain.tStart(prev) > chain.tEnd(cur))? chain.tStart(prev) - chain.tEnd(cur) : chain.tEnd(cur) - chain.tStart(prev);
		dist = min(qdist, tdist);
		if (chain.strand(cur) == chain.strand(prev) and dist >= 1000 and abs(chain.diag(cur) - chain.diag(prev)) <= ceil(0.15 * dist)) {  // missing TRA and INV
			if (push_new(genome, onec, lk, splitchains, splitchains_link, chain, cur)) {
				splitchains_link.push_back(0);
			}	
		}
		else if (chain.tStart(cur) > chain.tEnd(prev) + opts.splitdist or chain.tEnd(cur) + opts.splitdist < chain.tStart(prev)) {// TRA
			if (push_new(genome, onec, lk, splitchains, splitchains_link, chain, cur)) {
				splitchains_link.push_back(0);
			}
		}
		else if ((chain.link[im] == 1 and chain.strand(cur) == 0 and chain.strand(prev) == 0) or (chain.link[im] == 0 and chain.strand(cur)== 1 and chain.strand(prev) == 1)) { // DUP
			if (push_new(genome, onec, lk, splitchains, splitchains_link, chain, cur)) {
				splitchains_link.push_back(1);
				// if (chain.strand(cur) == 0) splitchains_link.push_back(0);
				// else splitchains_link.push_back(1);
			}			
		}
		else if ((chain.strand(cur) == 0 and chain.strand(prev) == 1) or (chain.strand(cur) == 1 and chain.strand(prev) == 0)) { // INV
			if (push_new(genome, onec, lk, splitchains, splitchains_link, chain, cur)) {
				splitchains_link.push_back(1);
			}				
		}
		else {
			onec.push_back(cur);
			lk.push_back(chain.link[im]);
		}	
		im++;
	}
	if (!onec.empty()) {
		push_new(genome, onec, lk, splitchains, splitchains_link, chain, cur);
	}

	for (int im = 0; im < splitchains.size(); im++) { // reverse the order for forward chain for refining chain
		if (splitchains[im].Strand == 0) { 
			reverse(splitchains[im].sptc.begin(), splitchains[im].sptc.end()); 
			reverse(splitchains[im].link.begin(), splitchains[im].link.end());
		}	
	}
}

// int 
// LargestUltimateChain(vector<UltimateChain> &ultimatechains) {
// 	int maxi = 0;
// 	for (int mi = 1; mi < ultimatechains.size(); mi++) {
// 		if (ultimatechains[mi].size() > ultimatechains[maxi].size()) {
// 			maxi = mi;
// 		}
// 	}
// 	return maxi;
// }

void 
output_unaligned(Read &read, const Options &opts, ostream &output) {
	// cerr << "unmapped: " << read.name << endl;
	if (opts.printFormat == "s") {
		Alignment unaligned = Alignment(read.seq, read.length, read.name, read.qual);
		unaligned.SimplePrintSAM(output, opts, read.passthrough);
	}
}

void 
OUTPUT(AlignmentsOrder &alignmentsOrder, Read &read, const Options &opts, Genome &genome, ostream *output){

	if (alignmentsOrder.size() > 0 and alignmentsOrder[0].SegAlignment.size() > 0) {
		int primary_num = 0;
		for (int a = 0; a < (int) min(alignmentsOrder.size(), opts.PrintNumAln); a++){
			for (int s = alignmentsOrder[a].SegAlignment.size() - 1; s >= 0; s--) {
				if (alignmentsOrder[a].SegAlignment[s]->Supplymentary == 0) primary_num++;
				alignmentsOrder[a].SegAlignment[s]->order = alignmentsOrder[a].SegAlignment.size() - 1 - s;
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
		assert(primary_num <= opts.PrintNumAln);
	}
	else if (read.unaligned == 1) {
	output_unaligned(read, opts, *output);
	}
}


void SimpleMapQV(AlignmentsOrder &alignmentsOrder, Read &read, const Options &opts) {
	if (read.unaligned) return;
	float q_coef;
	if (opts.bypassClustering) q_coef = 1.0f; // 40
	else q_coef = 22.0f; // 40
	int len = alignmentsOrder.size(); // number of primary aln and secondary aln
	for (int r = 0; r < len; r++) {
		if (r == 0 and len == 1) {
			// set mapq for each segment
			for (int s = alignmentsOrder[r].SegAlignment.size() - 1; s >= 0 ; s--) {
				float pen_cm_1;
				if (!opts.bypassClustering) {
					pen_cm_1 = (alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0 > 20? 1.0f : 0.05f ) * alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0;
					pen_cm_1 = (alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0 >= 5? 1.0f : 0.1f ) * pen_cm_1;					
				}
				else {
					pen_cm_1 = (alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0 > 10? 1.0f : 0.05f ) * alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0;
					pen_cm_1 = (alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0 >= 5? 1.0f : 0.02f ) * pen_cm_1;	// punish more
				}
				float identity = ((float) alignmentsOrder[r].SegAlignment[s]->nm ) / (//alignmentsOrder[r].SegAlignment[s]->nm + 
												      alignmentsOrder[r].SegAlignment[s]->nmm + 
												   alignmentsOrder[r].SegAlignment[s]->ndel +
												   alignmentsOrder[r].SegAlignment[s]->nins);
				identity = (identity < 1? identity : 1);
				float l = ( alignmentsOrder[r].SegAlignment[s]->value > 3? logf(alignmentsOrder[r].SegAlignment[s]->value / opts.globalK) : 0);
				long mapq;
				if (!opts.bypassClustering) mapq = (int)(pen_cm_1 * q_coef * l * identity);
				else mapq = (int)(pen_cm_1 * q_coef * identity);			
				// long mapq = (int)(pen_cm_1 * q_coef * l * identity);
				mapq = mapq > 0? mapq : 0;
				// if (1/identity >= 0.95f) {mapq = mapq < 60? mapq : 10;}
				alignmentsOrder[r].SegAlignment[s]->mapqv = mapq < 60? mapq : 60;
				if (r == 0 && len == 2 && alignmentsOrder[r].SegAlignment[s]->mapqv == 0) alignmentsOrder[r].SegAlignment[s]->mapqv = 1;
			}			
		}
		else if (r == 0 and len > 1) {
			// set mapq for each segment
			float x = (alignmentsOrder[r + 1].value) / (alignmentsOrder[r].value);
			for (int s = alignmentsOrder[r].SegAlignment.size() - 1; s >= 0 ; s--) {
				float pen_cm_1;
				if (!opts.bypassClustering) {
					pen_cm_1 = (alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0 > 20? 1.0f : 0.05f ) * alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0;
					pen_cm_1 = (alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0 >= 5? 1.0f : 0.1f ) * pen_cm_1;					
				}
				else {
					pen_cm_1 = (alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0 > 10? 1.0f : 0.05f ) * alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0;
					pen_cm_1 = (alignmentsOrder[r].SegAlignment[s]->NumOfAnchors0 >= 5? 1.0f : 0.02f ) * pen_cm_1;
				}

				float identity = (float) alignmentsOrder[r].SegAlignment[s]->nm / ( //alignmentsOrder[r].SegAlignment[s]->nm +
												   alignmentsOrder[r].SegAlignment[s]->nmm + 
												   alignmentsOrder[r].SegAlignment[s]->ndel +
												   alignmentsOrder[r].SegAlignment[s]->nins);
				float l = ( alignmentsOrder[r].SegAlignment[s]->value > 3? logf(alignmentsOrder[r].SegAlignment[s]->value / opts.globalK) : 0);
				identity = (identity < 1? identity : 1);
				long mapq;
				if (x >= 0.980f) {
					mapq = (int)(pen_cm_1 * (1.0f - x) * identity);
				}
				else if (!opts.bypassClustering) {
					mapq = (int)(pen_cm_1 * q_coef * (1.0f - x) * l * identity);
				}
				else {
					mapq = (int)(pen_cm_1 * q_coef * (1.0f - x) * identity);
				}
				// long mapq = (int)(pen_cm_1 * q_coef * (1.0f - x) * l);
				mapq -= (int)(4.343f * logf(len) + .499f);
				mapq = mapq > 0? mapq : 0;
				// if (1/identity >= 0.95f) {mapq = mapq < 60? mapq : 10;}
				alignmentsOrder[r].SegAlignment[s]->mapqv = mapq < 60? mapq : 60;
				if (r == 0 && len == 2 && alignmentsOrder[r].SegAlignment[s]->mapqv == 0) alignmentsOrder[r].SegAlignment[s]->mapqv = 1;
			}
		}
		else {
			for (int s = alignmentsOrder[r].SegAlignment.size() - 1; s >= 0 ; s--) {
				alignmentsOrder[r].SegAlignment[s]->mapqv = 0;
			}
		}
	}
}

void RemoveOverlappingClusters(vector<Cluster> &clusters, vector<int> &clusterOrder, const Options &opts) {
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

#endif
