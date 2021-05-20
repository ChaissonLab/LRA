#ifndef MAP_LOWACC_H_
#define MAP_LOWACC_H_
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
#include "ChainRefine.h"
#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cmath>	// std::log 
#include <sstream>
#include <thread>
#include <climits>
#include <map>
#include "RefineBreakpoint.h"

using namespace std;

void RemoveSpuriousSplitChain(vector<SplitChain> &chains, vector<bool> &spchain_link) {
	int total = 0;
	for (int i = 0; i < chains.size(); i++) {total += chains[i].size();}
	int filter = max((int) floor(0.02f * (float) total), 2);
	int filter_suspiciousDUPINV = max((int) floor(0.03f * (float) total), 2);
	// cerr << min(filter_suspiciousDUPINV, 4) << " " << min(filter, 2) << endl;

	vector<bool> remove(chains.size(), 0);
	for (int i = 0; i < chains.size(); i++) {
		if (chains[i].size() < min(filter, 2)) { // generally require >= 2 anchors
			remove[i] = 1;
		}
		if (i > 0 and spchain_link[i - 1] == 1 and chains[i].size() < min(filter_suspiciousDUPINV, 4)) { // for DUP and INV, requires more anchors
			remove[i] = 1;
		} 
	}
	int c = 0;
	for (int i = 0; i < remove.size(); i++) {
		if (!remove[i]) {
			chains[c] = chains[i];
			if (c > 1) spchain_link[c - 1] = spchain_link[i - 1];
			c++;
		}
	}
	chains.resize(c);
	if (c > 1) spchain_link.resize(c - 1);
	else spchain_link.clear();
	return;

}

int MapRead_lowacc(GenomePairs &forMatches, GenomePairs &revMatches, const vector<float> & LookUpTable, Read &read, Genome &genome, 
								vector<GenomeTuple> &genomemm, LocalIndex &glIndex, const Options &opts, ostream *output, ostream *svsigstrm,
								Timing &timing, IndelRefineBuffers &indelRefineBuffers, char *strands[2], char* readRC, pthread_mutex_t *semaphore=NULL) {
	//
	// bypass clustering splitting
	//
	vector<Cluster> clusters;
	//	cerr << "Cleaning " << endl;
	CleanMatches(forMatches, clusters, genome, read, opts, timing);
	CleanMatches(revMatches, clusters, genome, read, opts, timing, 1);
	//	cerr << "Cleaning done " << endl;
	forMatches.clear(); revMatches.clear();
	if (clusters.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
	}
	bool repetitivecluster = false;
	for (int s = 0; s < clusters.size(); s++) {	
		if (clusters[s].anchorfreq > 1.0f and clusters[s].anchorfreq <= 2.0f and clusters[s].matches.size() >= 500) {repetitivecluster = true;}
	}		
	if (opts.dotPlot and opts.readname == read.name) {
		ofstream cpclust("clusters-pre-remove.tab");
		for (int m = 0; m < clusters.size(); m++) {
			for (int n = 0; n < clusters[m].matches.size(); n++) {
				if (clusters[m].strand == 0) {
					cpclust << clusters[m].matches[n].first.pos << "\t"
						  << clusters[m].matches[n].second.pos + genome.header.pos[clusters[m].chromIndex]<< "\t"
						  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].second.pos + opts.globalK + genome.header.pos[clusters[m].chromIndex] << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].strand << "\t"
						  << clusters[m].anchorfreq<< endl;
				}
				else {
					cpclust << clusters[m].matches[n].first.pos << "\t"
						  << clusters[m].matches[n].second.pos + opts.globalK + genome.header.pos[clusters[m].chromIndex] << "\t"
						  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].second.pos + genome.header.pos[clusters[m].chromIndex] << "\t"
						  << m << "\t"
						  << genome.header.names[clusters[m].chromIndex]<< "\t"
						  << clusters[m].strand << "\t"
						  << clusters[m].anchorfreq<< endl;
				}				
			}
		}
		cpclust.close();
	}
	for (int s = 0; s < clusters.size(); s++) {	
		// Subtract chromOffSet from t coord.
		GenomePos chromOffset = genome.header.pos[clusters[s].chromIndex];
		for (int m = 0; m < clusters[s].matches.size(); m++) {
			clusters[s].matches[m].second.pos -= chromOffset;
		}
		clusters[s].tStart -= chromOffset;
		clusters[s].tEnd -= chromOffset;
	}
	vector<Cluster> ext_clusters(clusters.size());
	vector<UltimateChain> chains;
	//
	// Linear Extend on pure matches
	//
	for (int d = 0; d < clusters.size(); d++) {
		LinearExtend(&clusters[d].matches, ext_clusters[d].matches, ext_clusters[d].matchesLengths, opts, genome, read, clusters[d].chromIndex, clusters[d].strand, 1, opts.globalK);
		// ext_clusters[d].strand = clusters[d].strand;
		// ext_clusters[d].chromIndex = clusters[d].chromIndex;
		DecideCoordinates(ext_clusters[d], clusters[d].strand, clusters[d].chromIndex, clusters[d].anchorfreq);
	}
	//
	// Linear Extend efficiency
	//
	int o = 0, a = 0;
	for (int s = 0; s < clusters.size(); s++) {o += clusters[s].matches.size();}
	for (int s = 0; s < clusters.size(); s++) {a += ext_clusters[s].matches.size();}
	//cerr << "Linear Extend efficiency: " << (float) a/o << endl;
	for (int s = 0; s < ext_clusters.size(); s++) {	
		// Subtract chromOffSet from t coord.
		GenomePos chromOffset = genome.header.pos[clusters[s].chromIndex];
		for (int m = 0; m < ext_clusters[s].matches.size(); m++) {
			ext_clusters[s].matches[m].second.pos += chromOffset;
		}
		ext_clusters[s].tStart += chromOffset;
		ext_clusters[s].tEnd += chromOffset;
	}	
	clusters.clear();
	if (opts.dotPlot and opts.readname == read.name) {
		ofstream clust("Initial_ExtendClusters.tab", ofstream::app);
		for (int ep = 0; ep < ext_clusters.size(); ep++) {
			for (int eh = 0; eh < ext_clusters[ep].matches.size(); eh++) {	
				if (ext_clusters[ep].strand == 0) {
					clust << ext_clusters[ep].matches[eh].first.pos << "\t"
						  << ext_clusters[ep].matches[eh].second.pos << "\t"
						  << ext_clusters[ep].matches[eh].first.pos + ext_clusters[ep].matchesLengths[eh] << "\t"
						  << ext_clusters[ep].matches[eh].second.pos + ext_clusters[ep].matchesLengths[eh] << "\t"
						  << genome.header.names[ext_clusters[ep].chromIndex]<< "\t"
						  << ext_clusters[ep].strand << "\t"
						  << ep << endl;
				}
				else {
					clust << ext_clusters[ep].matches[eh].first.pos << "\t"
						  << ext_clusters[ep].matches[eh].second.pos + ext_clusters[ep].matchesLengths[eh] << "\t"
						  << ext_clusters[ep].matches[eh].first.pos + ext_clusters[ep].matchesLengths[eh] << "\t"
						  << ext_clusters[ep].matches[eh].second.pos<< "\t"
						  << genome.header.names[ext_clusters[ep].chromIndex]<< "\t"
						  << ext_clusters[ep].strand << "\t"
						  << ep << endl;					
				}
			}
		}
		clust.close();
	}
	//
	// SDP on matches
	//
	float match_rate = opts.initial_anchorbonus;
	if (repetitivecluster) {match_rate = 3;}
	Options optsSDP=opts;
	optsSDP.freeGap=4;	
	SparseDP(ext_clusters, chains, optsSDP, LookUpTable, read, match_rate);
	for (int p = 0; p < chains.size(); p++) { 
		RemoveSpuriousJump<UltimateChain>(chains[p]); 
		// chains[p].CleanSpurious();
	}

	if (chains.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
	} 		
	if (opts.dotPlot and opts.readname == read.name) {
		ofstream clust("Initial_SparseDP.tab", ofstream::app);
		for (int s = 0; s < chains.size(); s++) {
			for (int ep = 0; ep < chains[s].chain.size(); ep++) {
				if (chains[s].strand(ep) == 0) {
					clust << chains[s].qStart(ep) << "\t"
						  << chains[s].tStart(ep) << "\t"
						  << chains[s].qEnd(ep) << "\t"
						  << chains[s].tEnd(ep) << "\t"
						  << s << "\t"
						  << chains[s].ClusterNum(ep) << "\t"
						  << chains[s].strand(ep) << endl;
				}
				else {
					clust << chains[s].qStart(ep) << "\t"
						  << chains[s].tEnd(ep) << "\t"
						  << chains[s].qEnd(ep) << "\t"
						  << chains[s].tStart(ep) << "\t"
						  << s << "\t"
						  << chains[s].ClusterNum(ep) << "\t"
						  << chains[s].strand(ep) << endl;					
				}
			}				
		}
		clust.close();
	}	
	//
	// Close the space btwn matches
	//
	Options smallOpts=opts;
	Options tinyOpts=smallOpts;
	tinyOpts.globalMaxFreq=3;
	tinyOpts.maxDiag=5;
	tinyOpts.minDiagCluster=2;
	smallOpts.globalK=glIndex.k;
	smallOpts.globalW=glIndex.w;
	smallOpts.secondcoefficient+=3; // used to be 15
	smallOpts.globalMaxFreq=6;
	smallOpts.cleanMaxDiag=10;// used to be 25
	smallOpts.maxDiag=50;
	smallOpts.maxGapBtwnAnchors=100; // used to be 200 // 200 seems a little bit large
	smallOpts.minDiagCluster=3; // used to be 3
	tinyOpts.globalK=smallOpts.globalK-3;
	tinyOpts.globalW=tinyOpts.localW;
	//
	// Build local index for refining alignments.
	//
	LocalIndex forwardIndex(glIndex);
	LocalIndex reverseIndex(glIndex);
	LocalIndex *localIndexes[2] = {&forwardIndex, &reverseIndex};
	forwardIndex.IndexSeq(read.seq, read.length);
	reverseIndex.IndexSeq(readRC, read.length); 

	vector<SegAlignmentGroup> alignments;
	AlignmentsOrder alignmentsOrder(&alignments);
	AffineAlignBuffers buff;

	for (int p = 0; p < chains.size(); p++) {
		//
		// parse chain to splitchain (INV, DUP, TRA, Missing TRA and INV)
		//
		vector<SplitChain> spchain; vector<bool> spchain_link;
		SPLITChain(genome, read, chains[p], spchain, spchain_link, opts);
		RemoveSpuriousSplitChain(spchain, spchain_link); // remove spurious splitchain of <= 3 anchors
		if (p == 0 and spchain.size() == 0) {
			read.unaligned = 1; output_unaligned(read, opts, *output);
			return 0;
		} 	
		else if (p > 0 and spchain.size() == 0) break;
		if (opts.dotPlot and opts.readname == read.name) {
			ofstream clust("Initial_splitchain.tab", ofstream::app);
			for (int s= 0; s < spchain.size(); s++) {
			 	for (int ep= 0; ep < spchain[s].size(); ep++) {
					if (spchain[s].Strand == 0) {
						clust << spchain[s].genomepair(ep).first.pos << "\t"
							  << spchain[s].genomepair(ep).second.pos + genome.header.pos[spchain[s].CHROMIndex(genome)] << "\t"
							  << spchain[s].genomepair(ep).first.pos + spchain[s].length(ep) << "\t"
							  << spchain[s].genomepair(ep).second.pos + spchain[s].length(ep) + genome.header.pos[spchain[s].CHROMIndex(genome)] << "\t"
							  << s << "\t"
							  << p << "\t"
							  << spchain[s].Strand  << endl;
					}
					else {
						clust << spchain[s].genomepair(ep).first.pos  << "\t"
							  << spchain[s].genomepair(ep).second.pos + spchain[s].length(ep) + genome.header.pos[spchain[s].CHROMIndex(genome)] << "\t"
							  << spchain[s].genomepair(ep).first.pos + spchain[s].length(ep) << "\t"
							  << spchain[s].genomepair(ep).second.pos + genome.header.pos[spchain[s].CHROMIndex(genome)]<< "\t"
							  << s << "\t"
							  << p << "\t"
							  << spchain[s].Strand  << endl;					
					}
				}
			}				
			clust.close();
		}	
		//
		// refine splitchain
		//
		vector<Cluster> refined_clusters(spchain.size());
		Refine_splitchain(spchain, chains[p], refined_clusters, ext_clusters, genome, read, glIndex, localIndexes, smallOpts, opts);
		if (opts.dotPlot and opts.readname == read.name) {
			ofstream clust("RefinedClusters.tab", std::ofstream::app);
			for (int t = 0; t < refined_clusters.size(); t++) {
				for (int h = 0; h < refined_clusters[t].matches.size(); h++) {
					if (refined_clusters[t].strand  == 0) {
						clust << refined_clusters[t].matches[h].first.pos << "\t"
							  << refined_clusters[t].matches[h].second.pos + genome.header.pos[refined_clusters[t].chromIndex] << "\t"
							  << refined_clusters[t].matches[h].first.pos + smallOpts.globalK << "\t"
							  << refined_clusters[t].matches[h].second.pos + smallOpts.globalK  + genome.header.pos[refined_clusters[t].chromIndex] << "\t"
							  << t << "\t"
							  << p << "\t"
							  << genome.header.names[refined_clusters[t].chromIndex] <<"\t"
							  << refined_clusters[t].strand << endl;
					}
					else {
						clust << refined_clusters[t].matches[h].first.pos << "\t"
							  << refined_clusters[t].matches[h].second.pos + smallOpts.globalK  + genome.header.pos[refined_clusters[t].chromIndex] << "\t"
							  << refined_clusters[t].matches[h].first.pos + smallOpts.globalK << "\t"
							  << refined_clusters[t].matches[h].second.pos  + genome.header.pos[refined_clusters[t].chromIndex] << "\t"
							  << t << "\t"
							  << p << "\t"
							  << genome.header.names[refined_clusters[t].chromIndex] <<"\t"
							  << refined_clusters[t].strand << endl;					
					}
				}
			}
			clust.close();
		}	

		if (smallOpts.trimrefine) {
			SetGenomeOffset(spchain, genome);
			SubtractGenomeOffset(spchain);		
			TrimSplitChainDiagonal(spchain, refined_clusters);
			AddGenomeOffset(spchain);			
		}

		if (opts.dotPlot and opts.readname == read.name) {
			ofstream clust("RefinedClustersPostTrim.tab", std::ofstream::app);
			for (int t = 0; t < refined_clusters.size(); t++) {
				for (int h = 0; h < refined_clusters[t].matches.size(); h++) {
					assert(refined_clusters[t].matches[h].second.pos + smallOpts.globalK <= genome.lengths[refined_clusters[t].chromIndex]);
					if (refined_clusters[t].strand  == 0) {
						clust << refined_clusters[t].matches[h].first.pos << "\t"
							  << refined_clusters[t].matches[h].second.pos + genome.header.pos[refined_clusters[t].chromIndex] << "\t"
							  << refined_clusters[t].matches[h].first.pos + smallOpts.globalK << "\t"
							  << refined_clusters[t].matches[h].second.pos + smallOpts.globalK + genome.header.pos[refined_clusters[t].chromIndex]<< "\t"
							  << t << "\t"
							  << p << "\t"
							  << genome.header.names[refined_clusters[t].chromIndex] <<"\t"
							  << refined_clusters[t].strand << endl;
					}
					else {
						clust << refined_clusters[t].matches[h].first.pos << "\t"
							  << refined_clusters[t].matches[h].second.pos + smallOpts.globalK + genome.header.pos[refined_clusters[t].chromIndex]<< "\t"
							  << refined_clusters[t].matches[h].first.pos + smallOpts.globalK << "\t"
							  << refined_clusters[t].matches[h].second.pos + genome.header.pos[refined_clusters[t].chromIndex]<< "\t"
							  << t << "\t"
							  << p << "\t"
							  << genome.header.names[refined_clusters[t].chromIndex] <<"\t"
							  << refined_clusters[t].strand << endl;					
					}
				}
			}
			clust.close();
		}	
		
		//
		// refine between splitchain -- uncomment the next block!
		//
		vector<Cluster> RevBtwnCluster; // in case INV happens
		
		vector<tuple<int, int, int> > tracerev;
		Refine_Btwnsplitchain(spchain, refined_clusters, RevBtwnCluster, tracerev, genome, read, smallOpts, strands, spchain_link);

		/*
		//
		// Add back RevBtwnCluster; Edit spchain based on tracerev;
		//
		for (int t = 0; t < tracerev.size(); t++) {
			int I = get<1>(tracerev[t]);
			spchain.insert(spchain.begin() + I, SplitChain(get<2>(tracerev[t]) + refined_clusters.size()));
			spchain_link.insert(spchain_link.begin() + (I - 1), 1);
			spchain_link[I] = 1;
		}
		*/
		timing.Tick("FirstRefine");
		//
		// SplitChain spcluster -- each element points to a refined_cluster on the chain
		//
		SplitChain spcluster;
		spcluster.sptc.resize(spchain.size());
		spcluster.link = spchain_link;
		for (int t = 0; t < spchain.size(); t++) { spcluster.sptc[t] = spchain[t].clusterIndex;}
		assert(spcluster.sptc.size() == spcluster.link.size() + 1);
		spchain.clear(); spchain_link.clear();

		int a = refined_clusters.size();
		vector<Cluster *> Refined_Clusters(a);
		//vector<Cluster *> Refined_Clusters(a + RevBtwnCluster.size());
		for (int t = 0; t < refined_clusters.size(); t++) {Refined_Clusters[t] = &refined_clusters[t];}
		/*for (int t = 0; t < RevBtwnCluster.size(); t++) {Refined_Clusters[a + t] = &RevBtwnCluster[t];}*/
		if (p == 0 and Refined_Clusters.size() == 0) {
			read.unaligned = 1;
			output_unaligned(read, opts, *output);
			return 0;
		}	
		else if (p > 0 and Refined_Clusters.size() == 0) break;
		if (opts.dotPlot and opts.readname == read.name) {
			ofstream clust("WholeRefinedClusters.tab", std::ofstream::app);
			for (int t = 0; t < Refined_Clusters.size(); t++) {
				for (int h = 0; h < Refined_Clusters[t]->matches.size(); h++) {
					assert(Refined_Clusters[t]->matches[h].second.pos + smallOpts.globalK <= genome.lengths[Refined_Clusters[t]->chromIndex]);
					if (Refined_Clusters[t]->strand == 0) {
						clust << Refined_Clusters[t]->matches[h].first.pos << "\t"
							  << Refined_Clusters[t]->matches[h].second.pos + genome.header.pos[Refined_Clusters[t]->chromIndex] << "\t"
							  << Refined_Clusters[t]->matches[h].first.pos + smallOpts.globalK << "\t"
							  << Refined_Clusters[t]->matches[h].second.pos + smallOpts.globalK + genome.header.pos[Refined_Clusters[t]->chromIndex]<< "\t"
							  << t << "\t"
							  << p << "\t"
							  << genome.header.names[Refined_Clusters[t]->chromIndex] <<"\t"
							  << Refined_Clusters[t]->strand << endl;
					}
					else {
						clust << Refined_Clusters[t]->matches[h].first.pos << "\t"
							  << Refined_Clusters[t]->matches[h].second.pos + smallOpts.globalK + genome.header.pos[Refined_Clusters[t]->chromIndex] << "\t"
							  << Refined_Clusters[t]->matches[h].first.pos + smallOpts.globalK << "\t"
							  << Refined_Clusters[t]->matches[h].second.pos + genome.header.pos[Refined_Clusters[t]->chromIndex] << "\t"
							  << t << "\t"
							  << p << "\t"
							  << genome.header.names[Refined_Clusters[t]->chromIndex] <<"\t"
							  << Refined_Clusters[t]->strand << endl;					
					}
				}
			}
			clust.close();
		}	
		//
		// Merge adjacent Refined_Clusters
		//
		vector<Merge_SplitChain> mergeinfo; 
		SplitChain merge_spcluster;
		MergeChain(Refined_Clusters, mergeinfo, merge_spcluster, spcluster);
		assert(mergeinfo.size() == merge_spcluster.size());

		//
		// Linear extend matches 
		//
		int overlap = 0;
		vector<Cluster> extend_clusters;
		extend_clusters.resize(merge_spcluster.size());
		// for (int r = 0; r < Refined_Clusters.size(); r++) {
		// 	if (Refined_Clusters[r]->strand == 1) {
		// 		SwapStrand(read, opts, (*Refined_Clusters[r])); // LinearExtend is only for forward direction, so need to flip the reversed strand
		// 		Refined_Clusters[r]->flip = 1;
		// 	}
		// 	else {
		// 		Refined_Clusters[r]->flip = 0;
		// 	}
		// }
		for (int r = 0; r < mergeinfo.size(); r++) {
			bool st; int chromIndex; float anchorfreq;
			assert(mergeinfo[r].merged_clusterIndex.size() > 0);
			int cI;
			for (int t = 0; t < mergeinfo[r].merged_clusterIndex.size(); t++) {
				cI = mergeinfo[r].merged_clusterIndex[t];
				st = Refined_Clusters[cI]->strand; 
				chromIndex = Refined_Clusters[cI]->chromIndex;
				anchorfreq = Refined_Clusters[cI]->anchorfreq;
				LinearExtend(&(Refined_Clusters[cI]->matches), extend_clusters[r].matches, extend_clusters[r].matchesLengths, smallOpts, genome, read, chromIndex, st, 0, smallOpts.globalK);
			}
			DecideCoordinates(extend_clusters[r], st, chromIndex, anchorfreq);
			// if (st == 1) {
			// 	assert(Refined_Clusters[cI]->flip == 1);
			// 	SwapStrand(read, opts, extend_clusters[r]);
			// }
		}
		TrimOverlappedAnchors(extend_clusters, 0);

		int SizeRefinedClusters = 0, SizeExtendClusters = 0;
		for (int r = 0; r < Refined_Clusters.size(); r++) {
			SizeRefinedClusters += Refined_Clusters[r]->matches.size();
		}
		for (int ep = 0; ep < extend_clusters.size(); ep++) {
			SizeExtendClusters += extend_clusters[ep].matches.size();
		}	
		//		cerr << "Refined versus extend " << SizeRefinedClusters << "\t" << SizeExtendClusters << endl;
		
		if (p == 0 and SizeRefinedClusters == 0) {
			read.unaligned = 1;
			output_unaligned(read, opts, *output);
			return 0;
		}	
		else if (p > 0 and SizeRefinedClusters == 0) break;

		if (opts.dotPlot and opts.readname == read.name) {
			ofstream Eclust("ExtendClusters.tab", ofstream::app);
			for (int ep = 0; ep < extend_clusters.size(); ep++) {
				for (int eh = 0; eh < extend_clusters[ep].matches.size(); eh++) {	
					assert(extend_clusters[ep].matches[eh].second.pos + extend_clusters[ep].matchesLengths[eh] <= genome.lengths[extend_clusters[ep].chromIndex]);
					if (extend_clusters[ep].strand == 0) {
						Eclust << extend_clusters[ep].matches[eh].first.pos << "\t"
							  << extend_clusters[ep].matches[eh].second.pos + genome.header.pos[extend_clusters[ep].chromIndex] << "\t"
							  << extend_clusters[ep].matches[eh].first.pos + extend_clusters[ep].matchesLengths[eh] << "\t"
							  << extend_clusters[ep].matches[eh].second.pos + extend_clusters[ep].matchesLengths[eh] + genome.header.pos[extend_clusters[ep].chromIndex] << "\t"
							  << genome.header.names[extend_clusters[ep].chromIndex]<< "\t"
							  << extend_clusters[ep].strand << "\t"
							  << ep << "\t"
							  << p <<  endl;
					}
					else {
						Eclust << extend_clusters[ep].matches[eh].first.pos << "\t"
							  << extend_clusters[ep].matches[eh].second.pos + extend_clusters[ep].matchesLengths[eh] + genome.header.pos[extend_clusters[ep].chromIndex] << "\t"
							  << extend_clusters[ep].matches[eh].first.pos + extend_clusters[ep].matchesLengths[eh] << "\t"
							  << extend_clusters[ep].matches[eh].second.pos + genome.header.pos[extend_clusters[ep].chromIndex] << "\t"
							  << genome.header.names[extend_clusters[ep].chromIndex]<< "\t"
							  << extend_clusters[ep].strand << "\t"
							  << ep << "\t"
							  << p <<endl;					
					}
				}
			}
			Eclust.close();
		}	
		Refined_Clusters.clear();
		refined_clusters.clear();
		RevBtwnCluster.clear();
		//
		// SDP on clusters
		//
		timing.Tick("Extend");
		vector<UltimateChain> ultimatechains(merge_spcluster.size());
		for (int t = 0; t < merge_spcluster.size(); t++) {
			ultimatechains[t].clusters = &extend_clusters;
			ultimatechains[t].NumOfAnchors0 = chains[p].NumOfAnchors0;
			if (extend_clusters.size() > 0) { // and extend_clusters[0].matches.size() < 3*read.length

				SparseDP(merge_spcluster[t], extend_clusters, ultimatechains[t], smallOpts, LookUpTable, read);

				// ultimatechains[t].DebugCheck(read.length, genome);
				RemovePairedIndels<UltimateChain>(ultimatechains[t]); 
				RemoveSpuriousAnchors(ultimatechains[t]);
			}
			// ultimatechains[t].CleanSpurious();
		}

		if (opts.dotPlot and opts.readname == read.name) {
			ofstream Rclust("Refined_SparseDP.tab", ofstream::app);
			for (int t = 0; t < ultimatechains.size(); t++) {
				for (int s = 0; s < ultimatechains[t].chain.size(); s++) {
					if (ultimatechains[t].strand(s) == 0) {
						Rclust << ultimatechains[t].qStart(s) << "\t"
							  << ultimatechains[t].tStart(s) + genome.header.pos[ultimatechains[t].chromIndex(s)] << "\t"
							  << ultimatechains[t].qEnd(s) << "\t"
							  << ultimatechains[t].tEnd(s)+ genome.header.pos[ultimatechains[t].chromIndex(s)] << "\t"
							  << p << "\t" // chain_num
							  << t << "\t" // cluster_num
							  << ultimatechains[t].ClusterNum(s) << "\t"
							  << ultimatechains[t].strand(s) << endl;
					}
					else {
						Rclust << ultimatechains[t].qStart(s) << "\t"
							  << ultimatechains[t].tEnd(s) + genome.header.pos[ultimatechains[t].chromIndex(s)]<< "\t"
							  << ultimatechains[t].qEnd(s) << "\t"
							  << ultimatechains[t].tStart(s) + genome.header.pos[ultimatechains[t].chromIndex(s)]<< "\t"
							  << p << "\t" // chain_num
							  << t << "\t" // cluster_num
							  << ultimatechains[t].ClusterNum(s) << "\t"
							  << ultimatechains[t].strand(s) << endl;					
					}
				}				
			}
			Rclust.close();
		}	


		alignments.resize(alignments.size() + 1); 
		int LSC = LargestSplitChain(ultimatechains);
		LocalRefineAlignment(ultimatechains, extend_clusters, alignments, smallOpts, LookUpTable, read, strands, p, genome, LSC, tinyOpts, buff, svsigstrm);
		if (p == 0 and alignments.back().SegAlignment.size() == 0) {
			read.unaligned = 1;
			break;
		}
		for (int s = 0; s < alignments.back().SegAlignment.size(); s++) {
			if (opts.skipBandedRefine == false) { IndelRefineAlignment(read, genome, *alignments.back().SegAlignment[s], smallOpts, indelRefineBuffers); }
			if (s == 0 or s == alignments.back().SegAlignment.size() - 1) alignments.back().SegAlignment[s]->RetrieveEnd(s);// retrieve the ends
		}
		if (opts.refineBreakpoint) {
		  for (int s = 0; s < alignments.back().SegAlignment.size(); s++) {
		    if (s > 0) {
		    //
		    // segments are guaranteed to be in order from right to left on the read.
		    //
		      
		      RefineBreakpoint(read, genome, *alignments.back().SegAlignment[s], *alignments.back().SegAlignment[s-1], opts);
		    }
		  }
		}

		for (int s = 0; s < alignments.back().SegAlignment.size(); s++) {
		  alignments.back().SegAlignment[s]->CalculateStatistics(smallOpts, svsigstrm, LookUpTable);
		}
		alignments.back().SetFromSegAlignment(smallOpts);
		extend_clusters.clear();
		ultimatechains.clear();
		timing.Tick("SDP2");
	}
	if (read.unaligned or alignments.size() == 0) {
		output_unaligned(read, opts, *output);
		return 0;
	} 	
	alignmentsOrder.Update(&alignments);
	SimpleMapQV(alignmentsOrder, read, smallOpts);
	int e=timing.Elapsed();

	for (int a=0; a < alignments.size(); a++) {
	  for (int s=0; s< alignments[a].SegAlignment.size(); s++) {
	    alignments[a].SegAlignment[s]->runtime=e;
	  }
	}
	OUTPUT(alignmentsOrder, read, opts, genome, output);
	//
	// Done with one read. Clean memory.
	//
	// delete[] readRC;
	for (int a = 0; a < alignments.size(); a++) {
		for (int s = 0; s < alignments[a].SegAlignment.size(); s++) {
			delete alignments[a].SegAlignment[s];
		}
	}
	timing.Tick("Done");
	 
	if (alignments.size() > 0) return 1;
	return 0;
}

#endif
