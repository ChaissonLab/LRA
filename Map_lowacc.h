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

using namespace std;

void RemoveSpuriousSplitChain(vector<SplitChain> &chains, vector<bool> &spchain_link) {
	vector<bool> remove(chains.size(), 0);
	for (int i = 0; i < chains.size(); i++) {
		if (chains[i].size() < 3) {
			remove[i] == 1;
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
	spchain_link.resize(c - 1);
	return;

}

int MapRead_lowacc(GenomePairs &forMatches, GenomePairs &revMatches, const vector<float> & LookUpTable, Read &read, Genome &genome, 
								vector<GenomeTuple> &genomemm, LocalIndex &glIndex, Options &opts, ostream *output, ostream *svsigstrm,
								Timing &timing, IndelRefineBuffers &indelRefineBuffers, char *strands[2], char* readRC, pthread_mutex_t *semaphore=NULL) {
	//
	// bypass clustering splitting
	//
	vector<Cluster> clusters;
	CleanMatches(forMatches, clusters, genome, read, opts, timing);
	CleanMatches(revMatches, clusters, genome, read, opts, timing, 1);
	forMatches.clear(); revMatches.clear();
	if (clusters.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
	}
	bool NolinearCluster = true;
	for (int s = 0; s < clusters.size(); s++) {	
		if (clusters[s].anchorfreq <= 1.5f and clusters[s].matches.size() >= 10) {NolinearCluster = false; break;}
	}		
	if (NolinearCluster) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
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
						  << clusters[m].strand << "\t"
						  << clusters[m].anchorfreq<< endl;
				}
				else {
					cpclust << clusters[m].matches[n].first.pos << "\t"
						  << clusters[m].matches[n].second.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].first.pos + opts.globalK << "\t"
						  << clusters[m].matches[n].second.pos << "\t"
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
		LinearExtend(&clusters[d].matches, ext_clusters[d].matches, ext_clusters[d].matchesLengths, opts, genome, read, clusters[d].chromIndex, clusters[d].strand, 1);
		ext_clusters[d].strand = clusters[d].strand;
		ext_clusters[d].chromIndex = clusters[d].chromIndex;
		DecideCoordinates(ext_clusters[d]);
	}
	//
	// Linear Extend efficiency
	//
	int o = 0, a = 0;
	for (int s = 0; s < clusters.size(); s++) {o += clusters[s].matches.size();}
	for (int s = 0; s < clusters.size(); s++) {a += ext_clusters[s].matches.size();}
	//cerr << "Linear Extend efficiency: " << (float) a/o << endl;
	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
	//
	// SDP on matches
	//
	SparseDP(ext_clusters, chains, opts, LookUpTable, read);
	// for (int p = 0; p < chains.size(); p++) { 
	// 	RemovePairedIndels<UltimateChain>(chains[p]); 
	// 	chains[p].CleanSpurious();
	// }

	if (chains.size() == 0) {
		read.unaligned = 1;
		output_unaligned(read, opts, *output);
		return 0;
	} 		
	for (int s = 0; s < ext_clusters.size(); s++) {	
		// Subtract chromOffSet from t coord.
		GenomePos chromOffset = genome.header.pos[clusters[s].chromIndex];
		for (int m = 0; m < ext_clusters[s].matches.size(); m++) {
			ext_clusters[s].matches[m].second.pos -= chromOffset;
		}
		ext_clusters[s].tStart -= chromOffset;
		ext_clusters[s].tEnd -= chromOffset;
	}
	if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
		//
		// debug
		//
		for (int t = 0; t < spchain.size(); t++) {
			for (int l = 0; l < spchain[t].size(); l++) {
				int sl = spchain[t].chain->ClusterNum(spchain[t][0]);
				assert(spchain[t].chain->ClusterNum(spchain[t][l]) == sl);
			}
		}
		RemoveSpuriousSplitChain(spchain, spchain_link); // remove spurious splitchain of <= 3 anchors
		if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
			ofstream clust("Initial_splitchain.tab", ofstream::app);
			for (int s= 0; s < spchain.size(); s++) {
			 	for (int ep= 0; ep < spchain[s].size(); ep++) {
					if (spchain[s].Strand == 0) {
						clust << spchain[s].genomepair(ep).first.pos << "\t"
							  << spchain[s].genomepair(ep).second.pos << "\t"
							  << spchain[s].genomepair(ep).first.pos + spchain[s].length(ep) << "\t"
							  << spchain[s].genomepair(ep).second.pos + spchain[s].length(ep) << "\t"
							  << s << "\t"
							  << p << "\t"
							  << spchain[s].Strand  << endl;
					}
					else {
						clust << spchain[s].genomepair(ep).first.pos  << "\t"
							  << spchain[s].genomepair(ep).second.pos + spchain[s].length(ep) << "\t"
							  << spchain[s].genomepair(ep).first.pos + spchain[s].length(ep) << "\t"
							  << spchain[s].genomepair(ep).second.pos << "\t"
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
		// int o = ext_clusters.size(); ext_clusters.resize(o + 1); // dummy cluster
		// ext_clusters.back().matches.resize(1); 
		// chains[p].chain.push_back(0); chains[p].ClusterIndex.push_back(o); // dummy
		vector<Cluster> refined_clusters(spchain.size());
		Refine_splitchain(spchain, chains[p], refined_clusters, ext_clusters, genome, read, glIndex, localIndexes, smallOpts, opts);
		// chains[p].chain.pop_back(); chains[p].ClusterIndex.pop_back(); ext_clusters.pop_back(); // remove dummy
		//
		// refine between splitchain -- uncomment the next block!
		//
		
		vector<Cluster> RevBtwnCluster; // in case INV happens
		/*
		vector<tuple<int, int, int> > tracerev;
		Refine_Btwnsplitchain(spchain, refined_clusters, RevBtwnCluster, tracerev, genome, read, smallOpts, strands, spchain_link);
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

		//
		// SplitChain spcluster -- each element points to a refined_cluster on the chain
		//
		SplitChain spcluster;
		spcluster.sptc.resize(spchain.size());
		spcluster.link = spchain_link;
		for (int t = 0; t < spchain.size(); t++) { spcluster.sptc[t] = spchain[t].clusterIndex;}
		assert(spcluster.sptc.size() == spcluster.link.size() + 1);
		// spchain.clear(); spchain_link.clear();

		int a = refined_clusters.size();
		vector<Cluster *> Refined_Clusters(a + RevBtwnCluster.size());
		for (int t = 0; t < refined_clusters.size(); t++) {Refined_Clusters[t] = &refined_clusters[t];}
		/*for (int t = 0; t < RevBtwnCluster.size(); t++) {Refined_Clusters[a + t] = &RevBtwnCluster[t];}*/
		if (Refined_Clusters.size() == 0) {
			read.unaligned = 1;
			output_unaligned(read, opts, *output);
			return 0;
		}	
		if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
			ofstream clust("RefinedClusters.tab", std::ofstream::app);
			for (int t = 0; t < Refined_Clusters.size(); t++) {
				for (int h = 0; h < Refined_Clusters[t]->matches.size(); h++) {
					if (Refined_Clusters[t]->strand == 0) {
						clust << Refined_Clusters[t]->matches[h].first.pos << "\t"
							  << Refined_Clusters[t]->matches[h].second.pos << "\t"
							  << Refined_Clusters[t]->matches[h].first.pos + smallOpts.globalK << "\t"
							  << Refined_Clusters[t]->matches[h].second.pos + smallOpts.globalK << "\t"
							  << t << "\t"
							  << p << "\t"
							  << genome.header.names[Refined_Clusters[t]->chromIndex] <<"\t"
							  << Refined_Clusters[t]->strand << endl;
					}
					else {
						clust << Refined_Clusters[t]->matches[h].first.pos << "\t"
							  << Refined_Clusters[t]->matches[h].second.pos + smallOpts.globalK << "\t"
							  << Refined_Clusters[t]->matches[h].first.pos + smallOpts.globalK << "\t"
							  << Refined_Clusters[t]->matches[h].second.pos<< "\t"
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
		// Linear extend matches 
		//
		int overlap = 0;
		vector<Cluster> extend_clusters;
		extend_clusters.resize(spcluster.size());
		LinearExtend_chain(spcluster.sptc, extend_clusters, Refined_Clusters, smallOpts, genome, read, 0, overlap, 0);

		int SizeRefinedClusters = 0, SizeExtendClusters = 0;
		for (int r = 0; r < Refined_Clusters.size(); r++) {
			SizeRefinedClusters += Refined_Clusters[r]->matches.size();
		}
		for (int ep = 0; ep < extend_clusters.size(); ep++) {
			SizeExtendClusters += extend_clusters[ep].matches.size();
		}	
		if (SizeRefinedClusters == 0) {
			read.unaligned = 1;
			output_unaligned(read, opts, *output);
			return 0;
		}	

		if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
			ofstream Eclust("ExtendClusters.tab", ofstream::app);
			for (int ep = 0; ep < extend_clusters.size(); ep++) {
				for (int eh = 0; eh < extend_clusters[ep].matches.size(); eh++) {	
					if (extend_clusters[ep].strand == 0) {
						Eclust << extend_clusters[ep].matches[eh].first.pos << "\t"
							  << extend_clusters[ep].matches[eh].second.pos << "\t"
							  << extend_clusters[ep].matches[eh].first.pos + extend_clusters[ep].matchesLengths[eh] << "\t"
							  << extend_clusters[ep].matches[eh].second.pos + extend_clusters[ep].matchesLengths[eh] << "\t"
							  << genome.header.names[extend_clusters[ep].chromIndex]<< "\t"
							  << extend_clusters[ep].strand << "\t"
							  << ep << "\t"
							  << p <<  endl;
					}
					else {
						Eclust << extend_clusters[ep].matches[eh].first.pos << "\t"
							  << extend_clusters[ep].matches[eh].second.pos + extend_clusters[ep].matchesLengths[eh] << "\t"
							  << extend_clusters[ep].matches[eh].first.pos + extend_clusters[ep].matchesLengths[eh] << "\t"
							  << extend_clusters[ep].matches[eh].second.pos<< "\t"
							  << genome.header.names[extend_clusters[ep].chromIndex]<< "\t"
							  << extend_clusters[ep].strand << "\t"
							  << ep << "\t"
							  << p <<endl;					
					}
				}
			}
			Eclust.close();
		}	
		// Refined_Clusters.clear();
		// refined_clusters.clear();
		// RevBtwnCluster.clear();
		//
		// SDP on clusters
		//
		vector<UltimateChain> ultimatechains(spcluster.size());
		for (int t = 0; t < spcluster.size(); t++) {
			ultimatechains[t].clusters = &extend_clusters;
			SparseDP(spcluster[t], extend_clusters, ultimatechains[t], opts, LookUpTable, read);
			RemovePairedIndels<UltimateChain>(ultimatechains[t]); 
			// ultimatechains[t].CleanSpurious();
		}
		// 
		// SparseDP(spcluster, extend_clusters, ultimatechain, opts, LookUpTable, read);

		if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
			ofstream Rclust("Refined_SparseDP.tab", ofstream::app);
			for (int t = 0; t < ultimatechains.size(); t++) {
				for (int s = 0; s < ultimatechains[t].chain.size(); s++) {
					if (ultimatechains[t].strand(s) == 0) {
						Rclust << ultimatechains[t].qStart(s) << "\t"
							  << ultimatechains[t].tStart(s) << "\t"
							  << ultimatechains[t].qEnd(s) << "\t"
							  << ultimatechains[t].tEnd(s) << "\t"
							  << p << "\t" // chain_num
							  << s << "\t" // cluster_num
							  << ultimatechains[t].ClusterNum(s) << "\t"
							  << ultimatechains[t].strand(s) << endl;
					}
					else {
						Rclust << ultimatechains[t].qStart(s) << "\t"
							  << ultimatechains[t].tEnd(s) << "\t"
							  << ultimatechains[t].qEnd(s) << "\t"
							  << ultimatechains[t].tStart(s) << "\t"
							  << p << "\t" // chain_num
							  << s << "\t" // cluster_num
							  << ultimatechains[t].ClusterNum(s) << "\t"
							  << ultimatechains[t].strand(s) << endl;					
					}
				}				
			}
			Rclust.close();
		}	


		alignments.resize(alignments.size() + 1); 
		// vector<SplitChain> sp_ulchain; vector<bool> sp_ulink; vector<pair<GenomePos, GenomePos>> sp_ulpos;
		// SPLITChain(read, ultimatechain, sp_ulchain, sp_ulink, sp_ulpos, opts);
		int LSC = LargestUltimateChain(ultimatechains);
		LocalRefineAlignment(ultimatechains, extend_clusters, alignments, smallOpts, LookUpTable, read, strands, p, genome, LSC, tinyOpts, buff, svsigstrm);
		if (alignments.back().SegAlignment.size() == 0) {
			read.unaligned = 1;
			break;
		}
		for (int s = 0; s < alignments.back().SegAlignment.size(); s++) {
			if (opts.skipBandedRefine == false) { IndelRefineAlignment(read, genome, *alignments.back().SegAlignment[s], opts, indelRefineBuffers); }
			alignments.back().SegAlignment[s]->CalculateStatistics(smallOpts, svsigstrm, LookUpTable);
		}
		alignments.back().SetFromSegAlignment(smallOpts);
	}
	if (alignments.size() == 0) {
		output_unaligned(read, opts, *output);
		return 0;
	} 	
	alignmentsOrder.Update(&alignments);
	SimpleMapQV(alignmentsOrder, read, smallOpts);	
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
	if (alignments.size() > 0) {
		return 1;
	}
	return 0;
	// vector<SegAlignmentGroup> alignments;
	// AlignmentsOrder alignmentsOrder(&alignments);
	// AffineAlignBuffers buff;
	// for (int p = 0; p < chains.size(); p++) {
	// 	assert(chains[p].chain.size() > 0);
	// 	alignments.resize(alignments.size() + 1);	
	// 	vector<SplitChain> spchain; vector<bool> spchain_link; vector<pair<GenomePos, GenomePos>> spchain_qpos;
	// 	SPLITChain(read, chains[p], spchain, spchain_link, spchain_qpos, chains[p].link, opts);
	// 	int LSC = LargestSplitChain(spchain);
	// 	SparseDP_and_RefineAlignment_btwn_anchors(chains[p], spchain, spchain_link, spchain_qpos, ext_clusters, alignments, smallOpts, 
	// 							LookUpTable, read, strands, p, genome, LSC, tinyOpts, buff, svsigstrm);
	// 	if (alignments.back().SegAlignment.size() == 0) {
	// 		read.unaligned = 1;
	// 		break;
	// 	}
	// 	for (int s = 0; s < alignments.back().SegAlignment.size(); s++) {
	// 		if (opts.skipBandedRefine == false) { IndelRefineAlignment(read, genome, *alignments.back().SegAlignment[s], opts, indelRefineBuffers); }
	// 		alignments.back().SegAlignment[s]->CalculateStatistics(smallOpts, svsigstrm, LookUpTable);
	// 	}
	// 	alignments.back().SetFromSegAlignment(smallOpts);
	// }
	// if (read.unaligned) {
	// 	output_unaligned(read, opts, *output);
	// 	return 0;
	// } 	
	// alignmentsOrder.Update(&alignments);
	// alignmentsOrder.SimpleMapQV(read, smallOpts);	
	// OUTPUT(alignmentsOrder, read, opts, genome, output);
	// //
	// // Done with one read. Clean memory.
	// //
	// delete[] readRC;
	// for (int a = 0; a < alignments.size(); a++) {
	// 	for (int s = 0; s < alignments[a].SegAlignment.size(); s++) {
	// 		delete alignments[a].SegAlignment[s];
	// 	}
	// }
	// if (alignments.size() > 0) {
	// 	return 1;
	// }
	// return 0;
}

#endif