#ifndef MAP_HIGHACC_H_
#define MAP_HIGHACC_H_
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
#include "Mapping_ultility.h"

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

int MapRead_highacc(GenomePairs &forMatches, GenomePairs &revMatches, const vector<float> & LookUpTable, Read &read, Genome &genome, 
								vector<GenomeTuple> &genomemm, LocalIndex &glIndex, const Options &opts, ostream *output, ostream *svsigstrm,
								Timing &timing, IndelRefineBuffers &indelRefineBuffers, char *strands[2], char* readRC, pthread_mutex_t *semaphore=NULL) {
	vector<Cluster> clusters;
	MatchesToFineClusters(forMatches, clusters, genome, read, opts, timing);
	MatchesToFineClusters(revMatches, clusters, genome, read, opts, timing, 1);
	// if (opts.debug and opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
	// 	ofstream cpclust("clusters-pre-remove.tab");
	// 	for (int m = 0; m < clusters.size(); m++) {
	// 		for (int n = 0; n < clusters[m].matches.size(); n++) {
	// 			if (clusters[m].strand == 0) {
	// 				cpclust << clusters[m].matches[n].first.pos << "\t"
	// 					  << clusters[m].matches[n].second.pos << "\t"
	// 					  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
	// 					  << clusters[m].matches[n].second.pos + opts.globalK << "\t"
	// 					  << m << "\t"
	// 					  << genome.header.names[clusters[m].chromIndex]<< "\t"
	// 					  << clusters[m].strand << endl;
	// 			}
	// 			else {
	// 				cpclust << clusters[m].matches[n].first.pos << "\t"
	// 					  << clusters[m].matches[n].second.pos + opts.globalK << "\t"
	// 					  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
	// 					  << clusters[m].matches[n].second.pos << "\t"
	// 					  << m << "\t"
	// 					  << genome.header.names[clusters[m].chromIndex]<< "\t"
	// 					  << clusters[m].strand << endl;
	// 			}				
	// 		}
	// 	}
	// 	cpclust.close();
	// }
	//
	// Continue work on Clusters
	//
	if (opts.CheckTrueIntervalInFineCluster) {
		// CheckTrueIntervalInFineCluster(clusters, read.name, genome, read);
	}
	if (clusters.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
	}
	//cerr << "clusters.size(): " <<  clusters.size() << endl;

	forMatches.clear(); 
	revMatches.clear();

	// cerr << "before removing clusters.size(): " << clusters.size() << endl;
	// cerr << "before removing splitclusters.size(): " << splitclusters.size() << endl;
	// if (clusters.size() >= 4000) {
	// 	ClusterOrder fineClusterOrder(&clusters, 1);  // has some bug (delete clusters which should be kept -- cluster9_scaffold_58.fasta and cluster18_contig_234.fasta)
	// 	RemoveOverlappingClusters(clusters, fineClusterOrder.index, opts);
	// }
	// if (clusters.size() == 0) {
	// 	cerr << "unmapped " << read.name << endl;	
	// 	return 0; // This read cannot be mapped to the genome;
	// }

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
	// Split clusters on x and y coordinates, vector<Cluster> splitclusters, add a member for each splitcluster to specify the original cluster it comes from.
	// INPUT: vector<Cluster> clusters   OUTPUT: vector<Cluster> splitclusters with member--"coarse" specify the index of the original cluster splitcluster comes from
	//
	vector<Cluster> splitclusters;
	SplitClusters(clusters, splitclusters, read);
	DecideSplitClustersValue(clusters, splitclusters, opts, read);
	if (splitclusters.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
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
	if (splitclusters.size()/clusters.size() > 20) rate = rate / 2.0;// mapping to repetitive region

	// cerr << "clusters.size(): " << clusters.size() << " splitclusters.size(): " << splitclusters.size() << " rate: " << rate << " read.name: " << read.name << endl;

	SparseDP (splitclusters, Primary_chains, opts, LookUpTable, read, rate);
	// cerr << "1st SDP is done" << endl;
	if (Primary_chains.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
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
	switchindex(splitclusters, Primary_chains, clusters, genome, read);
	timing.Tick("Sparse DP - clusters");
	splitclusters.clear();
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
	if (lm == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;		
	}
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
	if (Primary_chains.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
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
	Options tinyOpts=smallOpts;
	tinyOpts.globalMaxFreq=3;
	tinyOpts.maxDiag=5;
	tinyOpts.minDiagCluster=2;

	//
	// Decide whether the number of anchors in each Cluster is enough to skip refining step;
	//
	bool sparse = 0;
	for (int p = 0; p < clusters.size(); p++) {
		if (((float)(clusters[p].matches.size())/(clusters[p].qEnd - clusters[p].qStart)) <= 0.05 and read.length <= 50000) sparse = 1;
	}
	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) cerr << "sparse: " << sparse << endl;
	//
	// Refining each cluster in "clusters" needed if CLR reads are aligned OR CCS reads with very few anchors. Otherwise, skip this step
	// After this step, the t coordinates in clusters and refinedclusters have been substract chromOffSet. 
	//
	vector<Cluster> refinedclusters(clusters.size());
 	vector<Cluster*> RefinedClusters(clusters.size());
 	//
 	// 1) read is not highly accurate; 2) read is highly accrate but not this read
 	//
	if (!opts.SkipLocalMinimizer and (opts.HighlyAccurate == false or (opts.HighlyAccurate == true and sparse == 1))) {
			
		smallOpts.globalK=glIndex.k;
		smallOpts.globalW=glIndex.w;
		smallOpts.secondcoefficient+=3; // used to be 15
		smallOpts.globalMaxFreq=10;
		smallOpts.cleanMaxDiag=10;// used to be 25
		smallOpts.maxDiag=50;
		smallOpts.maxGapBtwnAnchors=100; // used to be 200 // 200 seems a little bit large
		smallOpts.minDiagCluster=3; // used to be 3

		REFINEclusters(clusters, refinedclusters, genome, read, glIndex, localIndexes, smallOpts, opts);
		// cerr << "refine cluster done!" << endl;
		// refinedclusters have GenomePos, chromIndex, coarse, matches, strand, refinespace;
		for (int s = 0; s < clusters.size(); s++) {
			refinedclusters[s].anchorfreq = clusters[s].anchorfreq; // inherit anchorfreq
			RefinedClusters[s] = &refinedclusters[s];
		}
		clusters.clear();
	}
	else {
		tinyOpts.globalK=smallOpts.globalK-3;
		for (int s = 0; s < clusters.size(); s++) {	
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
	
	tinyOpts.globalK=smallOpts.globalK-3;
	tinyOpts.globalW=tinyOpts.localW;
	int K, W;
	if (opts.HighlyAccurate and sparse == 0) {K = opts.globalK; W = opts.globalW;}
	else {K = smallOpts.globalK; W = smallOpts.globalW;}

	if (RefinedClusters.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
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
	if (Primary_chains.size() == 0 or Primary_chains[0].chains.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
	}	
	//
	// For each chain, check the two ends and spaces between adjacent clusters. If the spaces are too wide, go to find anchors in the banded region.
	// For each chain, we have vector<Cluster> btwnClusters to store anchors;
	// RefineBtwnClusters and Do Liear extension for anchors
	//
	vector<Cluster> RevBtwnCluster;
	vector<tuple<int, int, int> > tracerev;
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			if (Primary_chains[p].chains[h].ch.size() == 0) continue;
				RefineBtwnClusters_chain(K, W, Primary_chains, RefinedClusters, RevBtwnCluster, tracerev, genome, read, smallOpts, p, h, strands);
		}
	}	

	if (opts.debug and opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
		ofstream clust("RefinedClusters.tab", std::ofstream::app);
		for (int p = 0; p < RefinedClusters.size(); p++) {
			for (int h = 0; h < RefinedClusters[p]->matches.size(); h++) {
				assert(RefinedClusters[p]->matches[h].first.pos + K <= read.length);
				if (RefinedClusters[p]->strand == 0) {
					clust << RefinedClusters[p]->matches[h].first.pos << "\t"
						  << RefinedClusters[p]->matches[h].second.pos + genome.header.pos[RefinedClusters[p]->chromIndex] << "\t"
						  << RefinedClusters[p]->matches[h].first.pos + K << "\t"
						  << RefinedClusters[p]->matches[h].second.pos + K + genome.header.pos[RefinedClusters[p]->chromIndex] << "\t"
						  << p << "\t"
						  << genome.header.names[RefinedClusters[p]->chromIndex] <<"\t"
						  << RefinedClusters[p]->strand << endl;
				}
				else {
					clust << RefinedClusters[p]->matches[h].first.pos << "\t"
						  << RefinedClusters[p]->matches[h].second.pos + K + genome.header.pos[RefinedClusters[p]->chromIndex] << "\t"
						  << RefinedClusters[p]->matches[h].first.pos + K << "\t"
						  << RefinedClusters[p]->matches[h].second.pos + genome.header.pos[RefinedClusters[p]->chromIndex]<< "\t"
						  << p << "\t"
						  << genome.header.names[RefinedClusters[p]->chromIndex] <<"\t"
						  << RefinedClusters[p]->strand << endl;					
				}
			}
		}
		clust.close();
	}	
	//
	// Add back RevBtwnCluster; Edit Primary_chains based on tracerev;
	//
	for (int p = 0; p < tracerev.size() ; p++) {
		int h = get<0>(tracerev[p]); 
		int I = get<1>(tracerev[p]); 
		Primary_chains[0].chains[h].ch.insert(Primary_chains[0].chains[h].ch.begin() + I, get<2>(tracerev[p]) + RefinedClusters.size());
		Primary_chains[0].chains[h].link.insert(Primary_chains[0].chains[h].link.begin() + (I - 1), 1);
		Primary_chains[0].chains[h].link[I] = 1;
	}
	int a = RefinedClusters.size();
	RefinedClusters.resize(a + RevBtwnCluster.size());
	for (int p = 0; p < RevBtwnCluster.size(); p++) {
		RefinedClusters[a + p] = &RevBtwnCluster[p];
	}

	timing.Tick("Refine_btwnclusters");

	vector<Cluster> extend_clusters;
	int overlap = 0;
	for (int p = 0; p < Primary_chains.size(); p++) {
		for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
			if (Primary_chains[p].chains[h].ch.size() == 0) continue;
				int cur_s = extend_clusters.size();
				extend_clusters.resize(cur_s + Primary_chains[p].chains[h].ch.size());
				LinearExtend_chain(Primary_chains[p].chains[h].ch, extend_clusters, RefinedClusters, smallOpts, genome, read, cur_s, overlap, 1, K);
		}
	}	
	timing.Tick("LinearExtend");

	int SizeRefinedClusters = 0, SizeExtendClusters = 0;
	for (int p = 0; p < RefinedClusters.size(); p++) {
		SizeRefinedClusters += RefinedClusters[p]->matches.size();
	}
	for (int ep = 0; ep < extend_clusters.size(); ep++) {
		SizeExtendClusters += extend_clusters[ep].matches.size();
	}	
	if (SizeRefinedClusters == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
	}	
	// cerr << "SizeRefinedClusters: " << SizeRefinedClusters << "   SizeExtendClusters: " << SizeExtendClusters << endl;
	// 	   << "  read.name:"<< read.name <<  endl;
	// cerr << "LinearExtend efficiency: " << (float)SizeExtendClusters/(float)SizeRefinedClusters << endl;
	// cerr << "overlapped anchors: " << overlap << " total:" << SizeExtendClusters <<endl;
	
	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
		ofstream clust("ExtendClusters.tab", ofstream::app);
		for (int ep = 0; ep < extend_clusters.size(); ep++) {
			for (int eh = 0; eh < extend_clusters[ep].matches.size(); eh++) {
				assert(extend_clusters[ep].matches[eh].first.pos + extend_clusters[ep].matchesLengths[eh] <= read.length);
				if (extend_clusters[ep].strand == 0) {
					clust << extend_clusters[ep].matches[eh].first.pos << "\t"
						  << extend_clusters[ep].matches[eh].second.pos + genome.header.pos[extend_clusters[ep].chromIndex]<< "\t"
						  << extend_clusters[ep].matches[eh].first.pos + extend_clusters[ep].matchesLengths[eh] << "\t"
						  << extend_clusters[ep].matches[eh].second.pos + extend_clusters[ep].matchesLengths[eh] + genome.header.pos[extend_clusters[ep].chromIndex]<< "\t"
						  << genome.header.names[extend_clusters[ep].chromIndex]<< "\t"
						  << extend_clusters[ep].strand << "\t"
						  << ep << endl;
				}
				else {
					clust << extend_clusters[ep].matches[eh].first.pos << "\t"
						  << extend_clusters[ep].matches[eh].second.pos + extend_clusters[ep].matchesLengths[eh] + genome.header.pos[extend_clusters[ep].chromIndex] << "\t"
						  << extend_clusters[ep].matches[eh].first.pos + extend_clusters[ep].matchesLengths[eh]<< "\t"
						  << extend_clusters[ep].matches[eh].second.pos + genome.header.pos[extend_clusters[ep].chromIndex]<< "\t"
						  << genome.header.names[extend_clusters[ep].chromIndex]<< "\t"
						  << extend_clusters[ep].strand << "\t"
						  << ep << endl;					
				}
			}
		}
		clust.close();
	}	

	clusters.clear();
	refinedclusters.clear();
	RefinedClusters.clear();
	RevBtwnCluster.size();

	vector<Cluster_SameDiag> samediag_clusters;
	MergeMatchesSameDiag(extend_clusters, samediag_clusters, opts); // There are a lot matches on the same diagonal, especially for contig;
	timing.Tick("Merged ExtendClusters");

	vector<SegAlignmentGroup> alignments;
	AlignmentsOrder alignmentsOrder(&alignments);
	AffineAlignBuffers buff;
	if (read.unaligned == 0) {
		int cur_cluster = 0;
		for (int p = 0; p < Primary_chains.size(); p++) {
			for (int h = 0; h < Primary_chains[p].chains.size(); h++) {
				if (Primary_chains[p].chains[h].ch.size() == 0) continue;
				alignments.resize(alignments.size() + 1);	

				vector<Cluster_SameDiag *> ExtendClusters(Primary_chains[p].chains[h].ch.size());
				for (int v = 0; v < Primary_chains[p].chains[h].ch.size(); v++) {
					// ExtendClusters[v] = &(extend_clusters[cur_cluster + v]);
					ExtendClusters[v] = &(samediag_clusters[cur_cluster + v]);
					assert(ExtendClusters[v]->coarse == cur_cluster + v);
				}

				//
				// Split the chain Primary_chains[p].chains[h] if clusters are aligned to different chromosomes; 
				// SplitAlignment is class that vector<* vector<Cluster>>
				// INPUT: vector<Cluster> ExtendClusters; OUTPUT:  vector<vector<unsigned int>> splitchain;
				//
				vector<SplitChain> splitchains;
				SPLITChain(read, ExtendClusters, splitchains, Primary_chains[p].chains[h].link, smallOpts);
				// cerr << "splitchains.size(): " << splitchains.size() << endl;
				int LSC = LargestSplitChain_dist(splitchains);
				//
				// Apply SDP on all splitchains to get the final rough alignment path;
				// store the result in GenomePairs tupChain; 
				// We need vector<Cluster> tupClusters for tackling anchors of different strands
				// NOTICE: Insert 4 points for anchors in the overlapping regions between Clusters;
				//
				//cerr << "splitchains.size(): " << splitchains.size()  << endl;
				LocalRefineAlignment(Primary_chains, splitchains, ExtendClusters, alignments, smallOpts, 
						     LookUpTable, read, strands, p, h, genome, LSC, tinyOpts, buff, svsigstrm, extend_clusters, false);
				ExtendClusters.clear();

				for (int s = 0; s < alignments.back().SegAlignment.size(); s++) {
					if (opts.skipBandedRefine == false) {
						IndelRefineAlignment(read, genome, *alignments.back().SegAlignment[s], opts, indelRefineBuffers, true);
					}
					//					if (s == 0 or s == alignments.back().SegAlignment.size() - 1) alignments.back().SegAlignment[s]->RetrieveEnd(s);// retrieve the ends
					alignments.back().SegAlignment[s]->CalculateStatistics(smallOpts, svsigstrm, LookUpTable); // final value gets compuated here
				}
				alignments.back().SetFromSegAlignment(smallOpts);
				cur_cluster += Primary_chains[p].chains[h].ch.size();
			}
			alignmentsOrder.Update(&alignments);
		}	
	}
	extend_clusters.clear();
	SimpleMapQV(alignmentsOrder, read, smallOpts);		

	timing.Tick("2nd SDP + local alignment");
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

	if (opts.storeTiming) {
		for (int a=0; a < (int) min(alignmentsOrder.size(), opts.PrintNumAln); a++){
			for (int s = 0; s < alignmentsOrder[a].SegAlignment.size(); s++) {
				alignmentsOrder[a].SegAlignment[s]->runtime=timing.Elapsed();
			}
		}
	}
	if (read.unaligned or alignments.size() == 0) {
		output_unaligned(read, opts, *output);
		return 0;
	} 	
	OUTPUT(alignmentsOrder, read, opts, genome, output);
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
	if (alignments.size() > 0) return 1;
	return 0;
}

#endif
