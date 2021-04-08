#ifndef LOCAL_REFINE_ALIGNMENT_H
#define LOCAL_REFINE_ALIGNMENT_H
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

#include <iostream>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <cmath>	// std::log 
#include <sstream>
#include <thread>
#include <climits>
#include <map>

//
// This function seperates finalchain by different strands; 
//
void debug_alignment (Alignment *alignment, int readlength, int genomeLen) {
	// This is for debugging
	if (alignment->blocks.size() > 1) {
		int last=alignment->blocks.size();
		assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length < readlength);
		assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length < genomeLen);
		assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].qPos);
		assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].tPos);
	}	
}

void
SeparateChainByStrand(UltimateChain & chain, vector<vector<int>> & finalSeperateChain) {
	vector<int> sep;
	int fl = 0;
	sep.push_back(fl);
	while (fl < chain.size() - 1) {

		if (chain.strand(fl) == chain.strand(fl+1)) {
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
	if (fl == chain.size() - 1) {
		sep.push_back(fl+1);
		finalSeperateChain.push_back(sep);
	}	
}

int 
LargestFinalSeperateChain(vector<vector<int>> &finalSeperateChain) {
	int maxi = 0;
	for (int mi = 1; mi < finalSeperateChain.size(); mi++) {
		if (finalSeperateChain[mi].back() - finalSeperateChain[mi][0] > finalSeperateChain[maxi].back() - finalSeperateChain[maxi][0]){
			maxi = mi;
		}
	}
	return maxi;	
}

int 
PassgenomeThres(int cur, GenomePos & genomeThres, UltimateChain & ultimatechain) {
	if (genomeThres != 0 and ultimatechain[cur].second.pos < genomeThres) return 1;
	else return 0;
}

int 
Matched(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te) {
	return min(qe-qs+1, te-ts+1); // TODO(Jingwen): check whether should add 1
}

void 
SetMatchAndGaps(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te, int &m, int &qg, int &tg) {
	m=Matched(qs, qe, ts, te);
	qg=qe-qs+1-m; //TODO(Jingwen): check whether should add 1
	tg=te-ts+1-m;
}

int 
AlignSubstrings(char *qSeq, GenomePos &qStart, GenomePos &qEnd, char *tSeq, GenomePos &tStart, GenomePos &tEnd,
						vector<int> &scoreMat, vector<Arrow> &pathMat, Alignment &aln, const Options &options, 
						AffineAlignBuffers &buff) {
	int qLen = qEnd-qStart;
	int tLen = tEnd-tStart;
	int drift = abs(qLen - tLen);
	int k = max(7, drift+1);
	
	/*
	int score = KBandAlign(&qSeq[qStart], qEnd-qStart, &tSeq[tStart], tEnd-tStart, -5,3,2,2, k, scoreMat, pathMat, aln);*/
	string readSeq(&qSeq[qStart], qEnd-qStart);
	string chromSeq(&tSeq[tStart],tEnd-tStart);
	int score = AffineOneGapAlign(readSeq, qLen, chromSeq, tLen, options.localMatch, options.localMismatch, options.localIndel, 
								min(drift*2+1,options.localBand), aln, buff);
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

void 
RefineSubstrings(char *read, GenomePos readSubStart, GenomePos readSubEnd, char *genome, GenomePos genomeSubStart, 
						GenomePos genomeSubEnd, vector<int> &scoreMat, vector<Arrow> &pathMat, Alignment &aln, const Options &opts,
						AffineAlignBuffers &buff) {
	aln.blocks.clear();
	AlignSubstrings(read, readSubStart, readSubEnd, genome, genomeSubStart, genomeSubEnd, scoreMat, pathMat, aln, opts, buff);
	for (int b = 0; b < aln.blocks.size(); b++) {
		aln.blocks[b].qPos += readSubStart;
		aln.blocks[b].tPos += genomeSubStart;

	}
}

void 
RefineByLinearAlignment(GenomePos &btc_curReadEnd, GenomePos &btc_curGenomeEnd, GenomePos &btc_nextReadStart, GenomePos &btc_nextGenomeStart, 
						bool str, int chromIndex, Alignment * alignment, Read & read, Genome & genome, char *strands[2], 
						vector<int> & scoreMat, vector<Arrow> & pathMat, const Options & opts, AffineAlignBuffers &buff) {
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
			// for (b = 1; b < betweenAnchorAlignment.blocks.size(); b++) {
			// 	assert(betweenAnchorAlignment.blocks[b-1].qPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].qPos);
			// 	assert(betweenAnchorAlignment.blocks[b-1].tPos + betweenAnchorAlignment.blocks[b-1].length <= betweenAnchorAlignment.blocks[b].tPos);						
			// }
			alignment->blocks.insert(alignment->blocks.end(), betweenAnchorAlignment.blocks.begin(), betweenAnchorAlignment.blocks.end());
			if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
				ofstream Lclust("LinearAlignment.tab", std::ofstream::app);
				for (b = 0; b < betweenAnchorAlignment.blocks.size(); b++) {
					Lclust << betweenAnchorAlignment.blocks[b].qPos << "\t"
						  << betweenAnchorAlignment.blocks[b].tPos + genome.header.pos[chromIndex]<< "\t"
						  << betweenAnchorAlignment.blocks[b].qPos + betweenAnchorAlignment.blocks[b].length << "\t"
						  << betweenAnchorAlignment.blocks[b].tPos + betweenAnchorAlignment.blocks[b].length + genome.header.pos[chromIndex]<< "\t"
						  << str << endl;
				}
				Lclust.close();
			}	
			// for (int bb = 1; bb < alignment->blocks.size(); bb++) {
			// 	assert(alignment->blocks[bb-1].qPos + alignment->blocks[bb-1].length <= read.length);
			// 	assert(alignment->blocks[bb-1].tPos + alignment->blocks[bb-1].length <= genome.lengths[chromIndex]);
			// 	assert(alignment->blocks[bb-1].qPos + alignment->blocks[bb-1].length <= alignment->blocks[bb].qPos);
			// 	assert(alignment->blocks[bb-1].tPos + alignment->blocks[bb-1].length <= alignment->blocks[bb].tPos);	
			// }
			betweenAnchorAlignment.blocks.clear();
		}
	}	
	// genomeThres = 0;	
}

void
SwitchToOriginalAnchors(FinalChain &prev, UltimateChain &cur, vector<Cluster_SameDiag *> &ExtendClusters, vector<Cluster> &extend_clusters) {
	cur = UltimateChain(&extend_clusters, prev.SecondSDPValue);
	for (int i = 0; i < prev.size(); i++) {
		int c = prev.ClusterNum(i);
		int k = prev.chain[i];
		for (int j = ExtendClusters[c]->end[k] - 1; j >= ExtendClusters[c]->start[k]; j--) {
			cur.chain.push_back(j);
			cur.ClusterIndex.push_back(ExtendClusters[c]->coarse);
		}			
	}
}

//
// This function creates the alignment for a "segment"(one big chunk) in the chain;
//
void
RefinedAlignmentbtwnAnchors(int &cur, int &next, bool &str, bool &inv_str, int &chromIndex, UltimateChain &chain, vector<SegAlignmentGroup> &alignments,
							Alignment *alignment, Read &read, Genome &genome, char *strands[2], vector<int> &scoreMat, 
							vector<Arrow> &pathMat, Options tinyOpts, AffineAlignBuffers &buff, 
							const vector<float> & LookUpTable, bool &inversion, bool &breakalignment, ostream *svsigstrm) {

	if (str == 0) alignment->blocks.push_back(Block(chain[cur].first.pos, chain[cur].second.pos, chain.length(cur)));
	else {
		alignment->blocks.push_back(Block(read.length - chain[cur].first.pos - chain.length(cur), chain[cur].second.pos, chain.length(cur)));
	}
	assert(chain[cur].first.pos + chain.length(cur) <= read.length);
	assert(chain[cur].second.pos + chain.length(cur) <= genome.lengths[chromIndex]);
	debug_alignment(alignment, read.length, genome.lengths[chromIndex]);
	//
	// Refine the alignment in the space between two adjacent anchors;
	//
	GenomePos curGenomeEnd, curReadEnd, nextGenomeStart, nextReadStart;
	if (str == 0) {
		curReadEnd = chain[cur].first.pos + chain.length(cur);
		nextReadStart = chain[next].first.pos;
		curGenomeEnd = chain[cur].second.pos + chain.length(cur);
		nextGenomeStart = chain[next].second.pos; 
	}
	else { // flip the space if it is reversed stranded
		curReadEnd = read.length - chain[cur].first.pos;
		nextReadStart = read.length - chain[next].first.pos - chain.length(next);
		curGenomeEnd = chain[cur].second.pos + chain.length(cur);
		nextGenomeStart = chain[next].second.pos;
	}
	assert(curReadEnd <= nextReadStart);
	assert(nextReadStart <= read.length);
	assert(nextGenomeStart <= genome.lengths[chromIndex]);
	
	if (curGenomeEnd <= nextGenomeStart) {
		long read_dist = nextReadStart - curReadEnd;
		long genome_dist = nextGenomeStart - curGenomeEnd;
		//
		// The linear alignment is more suitable to find a big gap; 
		//
		if (tinyOpts.RefineBySDP == true and min(read_dist, genome_dist) >= 300) {
			GenomePairs *BtwnPairs;
			GenomePairs for_BtwnPairs;
			GenomePairs rev_BtwnPairs;
			// cerr << "abs(read_dist-genome_dist): "  << abs(read_dist-genome_dist)<< endl;
			// cerr << "min(read_dist, genome_dist): " << min(read_dist, genome_dist) << endl;
			// cerr << "curReadEnd: " << curReadEnd << "  curGenomeEnd: " << curGenomeEnd << "  nextReadStart: " << nextReadStart << "  nextGenomeStart: " << nextGenomeStart << endl;
			int refineSpaceDiag = 0;
			//
			// Create a diagonal band that is not too big, not too small
			int sv_diag = max(read_dist, genome_dist) - min(read_dist, genome_dist);
			if (tinyOpts.readType == Options::contig or tinyOpts.readType == Options::ccs ) {
			  refineSpaceDiag = min((int) floor(max(80.f, 0.01f * read_dist)), 500);
			}
			else if (tinyOpts.readType == Options::clr or tinyOpts.readType == Options::ont) {
			  refineSpaceDiag = min((int) floor(max(100.f, 0.15f * read_dist)), 2000);	
			}
			refineSpaceDiag = max(2*sv_diag, refineSpaceDiag);

			// int refineSpaceDiag = (int) (0.15f * read_dist);	
			// if (refineSpaceDiag >= 100) cerr << "refineSpaceDiag: " << refineSpaceDiag << " read.name: " << read.name << endl;

			if (max(read_dist,genome_dist) < 100) {
				tinyOpts.globalK = 6;
				tinyOpts.localW  = 5;
			}
			else if (max(read_dist,genome_dist) < 500) {
				tinyOpts.globalK = 9;
				tinyOpts.localW  = 7;
			}
			else {
				tinyOpts.globalK = 13;
				tinyOpts.localW  = 7;
			}
				
 			RefineSpace(tinyOpts.globalK, tinyOpts.localW, refineSpaceDiag, 0, for_BtwnPairs, tinyOpts, genome, read, 
 					strands, chromIndex, nextReadStart, curReadEnd, nextGenomeStart, curGenomeEnd, str);
			//
			// If anchor is still too sparse, try finding seeds in the inversed direction
			//
			if ((for_BtwnPairs.size() / (float) min(read_dist, genome_dist)) < tinyOpts.anchorstoosparse and alignments.back().SegAlignment.back()->blocks.size() >=5 ) { // try inversion
				GenomePos temp = curReadEnd;
				curReadEnd = read.length - nextReadStart;
				nextReadStart = read.length - temp;
				RefineSpace(tinyOpts.globalK, tinyOpts.globalW, refineSpaceDiag, 0, rev_BtwnPairs, tinyOpts, genome, read, strands, chromIndex, nextReadStart, curReadEnd, nextGenomeStart, 
							curGenomeEnd, inv_str);	

				// if ((rev_BtwnPairs.size() / (float) min(read_dist, genome_dist)) < tinyOpts.anchorstoosparse and min(read_dist, genome_dist) >= 1000) {
				// 	for_BtwnPairs.clear();
				// 	// try refine in a larger band 1000 (sometimes INS and DEL)
 			// 		RefineSpace(tinyOpts.globalK, tinyOpts.localW, 1000, 0, for_BtwnPairs, tinyOpts, genome, read, 
 			// 				strands, chromIndex, nextReadStart, curReadEnd, nextGenomeStart, curGenomeEnd, str);	
 			// 		// break the alignment if still no anchors
 			// 		if ((for_BtwnPairs.size() / (float) min(read_dist, genome_dist)) < tinyOpts.anchorstoosparse) {
				// 		// break the alignment;
				// 		for_BtwnPairs.clear();
				// 		rev_BtwnPairs.clear();
				// 		breakalignment = 1;
				// 		inversion = 0;
				// 		return;		 						
 			// 		}			
				// }

				if ((rev_BtwnPairs.size() / (float) min(read_dist, genome_dist)) < tinyOpts.anchorstoosparse and min(read_dist, genome_dist) >= 1000) {
					// break the alignment;
					for_BtwnPairs.clear();
					rev_BtwnPairs.clear();
					breakalignment = 1;
					inversion = 0;
					return;					
				}
				if (for_BtwnPairs.size() >= rev_BtwnPairs.size()) {
					BtwnPairs = &for_BtwnPairs;
					rev_BtwnPairs.clear();
					inversion = 0;
					GenomePos temp = curReadEnd;
					curReadEnd = read.length - nextReadStart;
					nextReadStart = read.length - temp;
				}
				else {
					BtwnPairs = &rev_BtwnPairs;
					for_BtwnPairs.clear();
					inversion = 1;
				}				
			}
			else {
				BtwnPairs = &for_BtwnPairs;
			}

			if (BtwnPairs->size() > 0) {
				GenomePairs ExtendBtwnPairs;
				vector<int> ExtendBtwnPairsMatchesLength;
				vector<unsigned int> BtwnChain;
				LinearExtend(BtwnPairs, ExtendBtwnPairs, ExtendBtwnPairsMatchesLength, tinyOpts, genome, read, chromIndex, 0, 0, tinyOpts.globalK);

				if (inversion == 0) {
					//
					// insert the next anchor;
					//
					ExtendBtwnPairs.push_back(GenomePair(GenomeTuple(0, nextReadStart), GenomeTuple(0, nextGenomeStart)));
					ExtendBtwnPairsMatchesLength.push_back(chain.length(next));
					// 
					// insert the previous anchor;
					//
					ExtendBtwnPairs.push_back(GenomePair(GenomeTuple(0, curReadEnd-chain.length(cur)), GenomeTuple(0, curGenomeEnd-chain.length(cur))));
					ExtendBtwnPairsMatchesLength.push_back(chain.length(cur));			
				}

				TrimOverlappedAnchors(ExtendBtwnPairs, ExtendBtwnPairsMatchesLength);
				// nanoOpts.coefficient=12;//9 this is calibrately set to 12

				float inv_value = 0; int inv_NumofAnchors = 0;
				SparseDP_ForwardOnly(ExtendBtwnPairs, ExtendBtwnPairsMatchesLength, BtwnChain, tinyOpts, LookUpTable, inv_value, inv_NumofAnchors, 2); //1
				RemovePairedIndels(ExtendBtwnPairs, BtwnChain, ExtendBtwnPairsMatchesLength);

				if (tinyOpts.dotPlot and read.name == tinyOpts.readname) {
					ofstream pSclust("BtwnPairs.tab", std::ofstream::app);
					for (int bp = BtwnPairs->size()-1; bp >= 0; bp--) {
						if (inversion == 0) {
							if (str == 0) {
								pSclust << (*BtwnPairs)[bp].first.pos << "\t"
									  << (*BtwnPairs)[bp].second.pos + genome.header.pos[chromIndex] << "\t"
									  << (*BtwnPairs)[bp].first.pos + tinyOpts.globalK << "\t"
									  << (*BtwnPairs)[bp].second.pos + tinyOpts.globalK + genome.header.pos[chromIndex] << "\t"
									  << 0 << endl;
							}
							else if (str == 1) {
								pSclust << read.length - (*BtwnPairs)[bp].first.pos - tinyOpts.globalK << "\t"
									  << (*BtwnPairs)[bp].second.pos + tinyOpts.globalK + genome.header.pos[chromIndex]<< "\t"
									  << read.length - (*BtwnPairs)[bp].first.pos << "\t"
									  << (*BtwnPairs)[bp].second.pos + genome.header.pos[chromIndex] << "\t"
									  << 1 << endl;					
							}						
						} 
						else {
							if (inv_str == 0) {
								pSclust << (*BtwnPairs)[bp].first.pos << "\t"
									  << (*BtwnPairs)[bp].second.pos + genome.header.pos[chromIndex] << "\t"
									  << (*BtwnPairs)[bp].first.pos + tinyOpts.globalK << "\t"
									  << (*BtwnPairs)[bp].second.pos + tinyOpts.globalK + genome.header.pos[chromIndex]<< "\t"
									  << 0 << endl;
							}
							else if (inv_str == 1) {
								pSclust << read.length - (*BtwnPairs)[bp].first.pos - tinyOpts.globalK << "\t"
									  << (*BtwnPairs)[bp].second.pos + tinyOpts.globalK + genome.header.pos[chromIndex] << "\t"
									  << read.length - (*BtwnPairs)[bp].first.pos << "\t"
									  << (*BtwnPairs)[bp].second.pos + genome.header.pos[chromIndex] << "\t"
									  << 1 << endl;					
							}							
						}
					}
					pSclust.close();

					ofstream eSclust("ExtendBtwnPairs.tab", std::ofstream::app);
					for (int bp = ExtendBtwnPairs.size()-1; bp >= 0; bp--) {
						if (inversion == 0) {
							if (str == 0) {
								eSclust << ExtendBtwnPairs[bp].first.pos << "\t"
									  << ExtendBtwnPairs[bp].second.pos + genome.header.pos[chromIndex]<< "\t"
									  << ExtendBtwnPairs[bp].first.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] + genome.header.pos[chromIndex]<< "\t"
									  << 0 << endl;
							}
							else if (str == 1) {
								eSclust << read.length - ExtendBtwnPairs[bp].first.pos - ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] + genome.header.pos[chromIndex] << "\t"
									  << read.length - ExtendBtwnPairs[bp].first.pos << "\t"
									  << ExtendBtwnPairs[bp].second.pos + genome.header.pos[chromIndex] << "\t"
									  << 1 << endl;					
							}							
						}
						else {
							if (inv_str == 0) {
								eSclust << ExtendBtwnPairs[bp].first.pos << "\t"
									  << ExtendBtwnPairs[bp].second.pos + genome.header.pos[chromIndex] << "\t"
									  << ExtendBtwnPairs[bp].first.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] + genome.header.pos[chromIndex] << "\t"
									  << 0 << endl;
							}
							else if (inv_str == 1) {
								eSclust << read.length - ExtendBtwnPairs[bp].first.pos - ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] + genome.header.pos[chromIndex] << "\t"
									  << read.length - ExtendBtwnPairs[bp].first.pos << "\t"
									  << ExtendBtwnPairs[bp].second.pos + genome.header.pos[chromIndex]<< "\t"
									  << 1 << endl;					
							}							
						}
					}
					eSclust.close();
				}	

				if (tinyOpts.dotPlot and tinyOpts.readname == read.name) {
					ofstream Sclust("SparseDP_Forward.tab", std::ofstream::app);
						for (int btc = BtwnChain.size()-1; btc >= 0; btc--) {
							if (inversion == 0) {
								if (str == 0) {
									Sclust << ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + genome.header.pos[chromIndex]<< "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].first.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] + genome.header.pos[chromIndex]<< "\t"
										  << 0 << endl;
								}
								else if (str == 1) {
									Sclust << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos - ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] + genome.header.pos[chromIndex]<< "\t"
										  << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + genome.header.pos[chromIndex]<< "\t"
										  << 1 << endl;					
								}		
							}
							else {
								if (inv_str == 0) {
									Sclust << ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + genome.header.pos[chromIndex] << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].first.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] + genome.header.pos[chromIndex]<< "\t"
										  << 0 << endl;
								}
								else if (inv_str == 1) {
									Sclust << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos - ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] + genome.header.pos[chromIndex] << "\t"
										  << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + genome.header.pos[chromIndex] << "\t"
										  << 1 << endl;					
								}							
							}

						}
					Sclust.close();
				}
				//
				// Use linear alignment to ligand the gap between local SDP chain;
				//
				GenomePos btc_curReadEnd, btc_nextReadStart, btc_curGenomeEnd, btc_nextGenomeStart;
				btc_curReadEnd = curReadEnd;
				btc_curGenomeEnd = curGenomeEnd;
				int btc_end = BtwnChain.size()-1; int btc_start = 0;
				if (BtwnChain.back() == ExtendBtwnPairs.size()-1) btc_end = BtwnChain.size()-2;
				if (BtwnChain[0] == ExtendBtwnPairs.size()-2) btc_start = 1;

				if (inversion == 1) {
					alignment->UpdateParameters(str, tinyOpts, LookUpTable, svsigstrm, strands);
					Alignment *inv_alignment = new Alignment(inv_value, strands[str], read.seq, read.length, read.name, inv_str, read.qual, genome.seqs[chromIndex],  
																		 genome.lengths[chromIndex], genome.header.names[chromIndex], chromIndex); 
					inv_alignment->NumOfAnchors1 = BtwnChain.size();
					inv_alignment->NumOfAnchors0 = BtwnChain.size();
					inv_alignment->Supplymentary = 1;
					alignments.back().SegAlignment.push_back(inv_alignment);	
				}
				
				for (int btc = btc_end; btc >= btc_start; btc--) {
					btc_nextGenomeStart = ExtendBtwnPairs[BtwnChain[btc]].second.pos;
					btc_nextReadStart = ExtendBtwnPairs[BtwnChain[btc]].first.pos;
					RefineByLinearAlignment(btc_curReadEnd, btc_curGenomeEnd, btc_nextReadStart, btc_nextGenomeStart, str, chromIndex, 
											alignments.back().SegAlignment.back(), read, genome, strands, scoreMat, pathMat, tinyOpts, buff);				

					alignments.back().SegAlignment.back()->blocks.push_back(Block(ExtendBtwnPairs[BtwnChain[btc]].first.pos, ExtendBtwnPairs[BtwnChain[btc]].second.pos, 
														ExtendBtwnPairsMatchesLength[BtwnChain[btc]]));
					debug_alignment(alignments.back().SegAlignment.back(), read.length, genome.lengths[chromIndex]);
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
											alignments.back().SegAlignment.back(), read, genome, strands, scoreMat, pathMat, tinyOpts, buff);					
				}				
			}
			else {
				RefineByLinearAlignment(curReadEnd, curGenomeEnd, nextReadStart, nextGenomeStart, str, chromIndex, 
									alignment, read, genome, strands, scoreMat, pathMat, tinyOpts, buff);
			}
		}
		else {
			RefineByLinearAlignment(curReadEnd, curGenomeEnd, nextReadStart, nextGenomeStart, str, chromIndex, 
									alignment, read, genome, strands, scoreMat, pathMat, tinyOpts, buff);				
		}
	}
	// else genomeThres = curGenomeEnd;
	return;
}

void
LocalRefineAlignment(vector<Primary_chain> &Primary_chains, vector<SplitChain> &splitchains, vector<Cluster_SameDiag *> &ExtendClusters, 
		vector<SegAlignmentGroup> &alignments, const Options &smallOpts, const vector<float> & LookUpTable, Read &read, char *strands[2], int &p, int &h, 
		     Genome &genome, int &LSC, const Options &tinyOpts, AffineAlignBuffers &buff, ostream *svsigstrm, vector<Cluster> &extend_clusters, bool refineEnd=true) {
	for (int st = 0; st < splitchains.size(); st++) {
		//
		// Apply SparseDP on extended anchors on every splitchain;
		// INPUT: vector<unsigned int> splitchain, vector<Cluster> ExtendClusters; OUTPUT: FinalChain finalchain;
		//
	  FinalChain finalchain(&ExtendClusters);
		SparseDP(splitchains[st], ExtendClusters, finalchain, smallOpts, LookUpTable, read);
		//timing.Tick("2nd SDP");
		
		if (tinyOpts.RemovePairedIndels) {
		  RemoveSmallPairedIndels<FinalChain> (finalchain);
		  RemovePairedIndels<FinalChain>(finalchain, refineEnd); 
		}
		if (tinyOpts.RemoveSpuriousAnchors) RemoveSpuriousAnchors(finalchain);
		//
		// switch back to the original anchors;
		//
		UltimateChain ultimatechain; 
		SwitchToOriginalAnchors(finalchain, ultimatechain, ExtendClusters, extend_clusters);
		finalchain.clear();

		//cerr << "2nd SDP done!" << endl;
		if (ultimatechain.size() == 0) continue; // cannot be mapped to the genome!
		//
		// Refine and store the alignment; NOTICE: when filling in the space between two adjacent anchors, 
		// the process works in forward direction, so we need to flip the small matches
		// found in the spaces before insert them into the alignment if necessary;
		// INPUT: vector<int> finalchainSeperateStrand; OUTPUT: Alignment *alignment;
		//
		int ad = alignments.back().SegAlignment.size();
		int start = 0, end = ultimatechain.size();
		assert(start < end); 
		bool str = ultimatechain.strand(start);
		int cln = ultimatechain.ClusterNum(start);
		int chromIndex = extend_clusters[cln].chromIndex;

		if (smallOpts.dotPlot and read.name == smallOpts.readname ) {
			ofstream clust("SparseDP.tab", ofstream::app);
			for (int ep = 0; ep < ultimatechain.size(); ep++) {
				if (ultimatechain.strand(ep) == 0) {
					clust << ultimatechain[ep].first.pos << "\t"
						  << ultimatechain[ep].second.pos + genome.header.pos[chromIndex] << "\t"
						  << ultimatechain[ep].first.pos + ultimatechain.length(ep) << "\t"
						  << ultimatechain[ep].second.pos + ultimatechain.length(ep) + genome.header.pos[chromIndex]<< "\t"
						  << st << "\t"
						  << ultimatechain.ClusterNum(ep) << "\t"
						  << ultimatechain.strand(ep) << endl;
				}
				else {
					clust << ultimatechain[ep].first.pos << "\t"
						  << ultimatechain[ep].second.pos + ultimatechain.length(ep) + genome.header.pos[chromIndex]<< "\t"
						  << ultimatechain[ep].first.pos + ultimatechain.length(ep) << "\t"
						  << ultimatechain[ep].second.pos + genome.header.pos[chromIndex] << "\t"
						  << st << "\t"
						  << ultimatechain.ClusterNum(ep) << "\t"
						  << ultimatechain.strand(ep) << endl;					
				}
			}
			clust.close();
		}	
		Alignment *alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, 
											 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
											 genome.header.names[chromIndex], chromIndex); 
		alignment->NumOfAnchors0 = Primary_chains[p].chains[h].NumOfAnchors0;
		alignment->NumOfAnchors1 = end - start;
		alignment->SecondSDPValue = ultimatechain.SecondSDPValue;
		alignments.back().SegAlignment.push_back(alignment);
		vector<int> scoreMat;
		vector<Arrow> pathMat;
		//
		// Set the secondary or supplymentary flag for the alignment; 
		//
		if (h > 0) alignment->ISsecondary = 1;
		if (st != LSC) alignment->Supplymentary = 1;
		bool inv_str = 0, inversion = 0, breakalignment = 0, prev_inv = 0;
		if (str == 0) {
			int last = end - 1;
			int fl = end - 1;
			inv_str = 1; prev_inv = end - 1;
			while (fl > start) {
				int cur = fl;
				int next = fl - 1;
				//
				// Check the distance between two anchors. If the distance is too large, refine anchors and apply 3rd SDP;
				// Otherwise apply linear alignment. 
				//
				RefinedAlignmentbtwnAnchors(cur, next, str, inv_str, chromIndex, ultimatechain, alignments, alignments.back().SegAlignment.back(), 
											read, genome, strands, scoreMat, pathMat, tinyOpts, buff, 
											LookUpTable, inversion, breakalignment, svsigstrm);
				if (inversion == 1 or breakalignment == 1) { // an inversion found in between two anchors; Need to make a new alignment for the next
					// alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
					if (inversion == 1) alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
					else alignments.back().SegAlignment.back()->UpdateParameters(str, smallOpts, LookUpTable, svsigstrm, strands);
					alignments.back().SegAlignment.back()->NumOfAnchors0 = Primary_chains[p].chains[h].NumOfAnchors0;
					alignments.back().SegAlignment.back()->NumOfAnchors1 = last - fl;
					Alignment *next_alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, 
											 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
											 genome.header.names[chromIndex], chromIndex); 
					// next_alignment->NumOfAnchors0 = Primary_chains[p].chains[h].NumOfAnchors0;
					// next_alignment->NumOfAnchors1 = fl - start;
					next_alignment->Supplymentary = 1;
					alignments.back().SegAlignment.push_back(next_alignment);
					last = fl;
					prev_inv = fl;
					inversion = 0;
					breakalignment = 0;
				} 
				fl--;
			}
			alignments.back().SegAlignment.back()->NumOfAnchors0 = Primary_chains[p].chains[h].NumOfAnchors0;
			alignments.back().SegAlignment.back()->NumOfAnchors1 = last - fl;
			alignments.back().SegAlignment.back()->blocks.push_back(Block(ultimatechain[start].first.pos, ultimatechain[start].second.pos, 
												ultimatechain.length(start)));
		}
		else {
			int last = start;
			int fl = start; 
			inv_str = 0; prev_inv = start;
			while (fl < end - 1) {
				int cur = fl;
				int next = fl + 1;
				RefinedAlignmentbtwnAnchors(cur, next, str, inv_str, chromIndex, ultimatechain, alignments, alignments.back().SegAlignment.back(), 
											read, genome, strands, scoreMat, pathMat, tinyOpts, buff, 
											LookUpTable, inversion, breakalignment, svsigstrm);
				if (inversion == 1 or breakalignment == 1) {
					if (inversion == 1) alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
					else alignments.back().SegAlignment.back()->UpdateParameters(str, smallOpts, LookUpTable, svsigstrm, strands);
					// alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
					alignments.back().SegAlignment.back()->NumOfAnchors1 = fl - last;
					alignments.back().SegAlignment.back()->NumOfAnchors0 = Primary_chains[p].chains[h].NumOfAnchors0;
					Alignment *next_alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, 
											 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
											 genome.header.names[chromIndex], chromIndex); 
					// next_alignment->NumOfAnchors0 = Primary_chains[p].chains[h].NumOfAnchors0;
					// next_alignment->NumOfAnchors1 = end-fl;
					next_alignment->Supplymentary = 1;
					alignments.back().SegAlignment.push_back(next_alignment);
					last = fl;
					prev_inv = fl;
					inversion = 0;
					breakalignment = 0;
				} 
				fl++;
			}	
			alignments.back().SegAlignment.back()->NumOfAnchors0 = Primary_chains[p].chains[h].NumOfAnchors0;
			alignments.back().SegAlignment.back()->NumOfAnchors1 = fl - last;
			alignments.back().SegAlignment.back()->blocks.push_back(Block(read.length - ultimatechain[end - 1].first.pos - ultimatechain.length(end - 1), 
												ultimatechain[end - 1].second.pos, ultimatechain.length(end - 1)));
		}
		// // debug checking
		// for (int bb = 1; bb < alignments.back().SegAlignment.back()->blocks.size(); bb++) {
		// 	assert(alignments.back().SegAlignment.back()->blocks[bb-1].qPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= read.length);
		// 	assert(alignments.back().SegAlignment.back()->blocks[bb-1].tPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= genome.lengths[chromIndex]);
		// 	assert(alignments.back().SegAlignment.back()->blocks[bb-1].qPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= alignments.back().SegAlignment.back()->blocks[bb].qPos);
		// 	assert(alignments.back().SegAlignment.back()->blocks[bb-1].tPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= alignments.back().SegAlignment.back()->blocks[bb].tPos);	
		// }
		//
		// Set some parameters in Alignment alignment
		//
		alignments.back().SegAlignment.back()->UpdateParameters(str, smallOpts, LookUpTable, svsigstrm, strands);
		ultimatechain.clear();
		//
		// Decide the inversion flag for each segment alignment (PAF format); seek +,-,+ or -,+,-
		//
		int js = ad + 2;
		while (js < alignments.back().SegAlignment.size()) {
			if ((alignments.back().SegAlignment[js-2]->strand == 0 and alignments.back().SegAlignment[js-1]->strand == 1
				and alignments.back().SegAlignment[js]->strand == 0) or ((alignments.back().SegAlignment[js-2]->strand == 1 
				and alignments.back().SegAlignment[js-1]->strand == 0 and alignments.back().SegAlignment[js]->strand == 1))) {
				//
				// inversion cannot be too far (<10,000) from alignments at both sides;
				//
				if (alignments.back().SegAlignment[js-1]->tStart > alignments.back().SegAlignment[js-2]->tEnd + 10000
					or alignments.back().SegAlignment[js]->tStart > alignments.back().SegAlignment[js-1]->tEnd + 10000) {
					js++;
					continue;
				}
				//
				// length requirements; 
				//
				else if (alignments.back().SegAlignment[js]->nm < 500 or alignments.back().SegAlignment[js-1]->nm < 500
						or alignments.back().SegAlignment[js-2]->nm < 40 or alignments.back().SegAlignment[js-2]->nm > 15000) {
					js++;
					continue;
				}
				else if (alignments.back().SegAlignment[js-2]->typeofaln != 3) {
					alignments.back().SegAlignment[js-1]->typeofaln=3;
				}
			}
			js++;
		}

	}	
}

// void RefindEnds(GenomePos &qPos, int cur, UltimateChain &chain, bool str, Alignment *alignment, const Options &opts, 
// 				Genome &genome, Read &read, int &chromIndex, char *strands[2], vector<int> &scoreMat, 
// 				vector<Arrow> &pathMat, AffineAlignBuffers &buff, const vector<float> & LookUpTable, bool block = 1) { // block = 1: right, block = 0: left
	
// 	GenomePos curGenomeEnd, curReadEnd, nextGenomeStart, nextReadStart;
// 	if (block) { // right
// 		if (qPos <= chain.qEnd(cur) + 100) return;
// 		if (str == 0) { // forward
// 			curReadEnd = chain.qEnd(cur); nextReadStart = qPos - 100; 
// 			curGenomeEnd = chain.tEnd(cur); nextGenomeStart = curGenomeEnd + (nextReadStart - curReadEnd);
// 		}
// 		else { // flip the space if it is reversed stranded
// 			curReadEnd = read.length - (qPos - 100); nextReadStart = read.length - chain.qEnd(cur); 
// 			nextGenomeStart = chain.tStart(cur); curGenomeEnd = nextGenomeStart - (nextReadStart - curReadEnd);
// 		}		
// 	}
// 	else { // left
// 		if (chain.qStart(cur) <= qPos + 100) return;
// 		if (str == 0) {
// 			curReadEnd = qPos + 100; nextReadStart = chain.qStart(cur); 
// 			curGenomeEnd = chain.tStart(cur) - (nextReadStart - curReadEnd); nextGenomeStart = chain.tStart(cur);
// 		}
// 		else {
// 			curReadEnd = read.length - chain.qStart(cur); nextReadStart = read.length - (qPos + 100); 
// 			curGenomeEnd = chain.tStart(cur) - (nextReadStart - curReadEnd); nextGenomeStart = chain.tStart(cur);
// 		}

// 	}
// 	assert(curGenomeEnd <= nextGenomeStart); 
// 	assert(curReadEnd <= nextReadStart);
// 	if (curGenomeEnd >= nextGenomeStart + 100 or curReadEnd >= nextReadStart + 100) return;
// 	GenomePairs EndPairs;
// 	RefineSpace(opts.globalK, opts.globalW, 20, 0, EndPairs, opts, genome, read, strands, chromIndex, nextReadStart, curReadEnd, nextGenomeStart, curGenomeEnd, str);	
// 	if (EndPairs.size() > 0 and (EndPairs.size() / (float) (nextReadStart - curReadEnd)) >= opts.anchorstoosparse ) {
// 		GenomePairs ExtendEndPairs; vector<int> ExtendEndPairsMatchesLength; vector<unsigned int> EndChain;
// 		LinearExtend(&EndPairs, ExtendEndPairs, ExtendEndPairsMatchesLength, opts, genome, read, chromIndex, 0, 0, opts.globalK);	
// 		TrimOverlappedAnchors(ExtendEndPairs, ExtendEndPairsMatchesLength);
// 		float value = 0; int NumofAnchors = 0;
// 		SparseDP_ForwardOnly(ExtendEndPairs, ExtendEndPairsMatchesLength, EndChain, opts, LookUpTable, value, NumofAnchors, 1/2); //1
// 		RemovePairedIndels(ExtendEndPairs, EndChain, ExtendEndPairsMatchesLength);
// 		if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
// 			ofstream eSclust("ExtendBtwnPairs.tab", std::ofstream::app);
// 			for (int bp = ExtendEndPairs.size()-1; bp >= 0; bp--) {
// 				if (str == 0) {
// 					eSclust << ExtendEndPairs[bp].first.pos << "\t"
// 						  << ExtendEndPairs[bp].second.pos << "\t"
// 						  << ExtendEndPairs[bp].first.pos + ExtendEndPairsMatchesLength[bp] << "\t"
// 						  << ExtendEndPairs[bp].second.pos + ExtendEndPairsMatchesLength[bp] << "\t"
// 						  << 0 << endl;
// 				}
// 				else if (str == 1) {
// 					eSclust << read.length - ExtendEndPairs[bp].first.pos - ExtendEndPairsMatchesLength[bp] << "\t"
// 						  << ExtendEndPairs[bp].second.pos + ExtendEndPairsMatchesLength[bp] << "\t"
// 						  << read.length - ExtendEndPairs[bp].first.pos << "\t"
// 						  << ExtendEndPairs[bp].second.pos << "\t"
// 						  << 1 << endl;					
// 				}							
// 			}
// 			eSclust.close();
// 			if (opts.dotPlot and read.name == opts.readname) {
// 				ofstream Sclust("SparseDP_Forward.tab", std::ofstream::app);
// 					for (int btc = EndChain.size()-1; btc >= 0; btc--) {
// 						if (str == 0) {
// 							Sclust << ExtendEndPairs[EndChain[btc]].first.pos << "\t"
// 								  << ExtendEndPairs[EndChain[btc]].second.pos << "\t"
// 								  << ExtendEndPairs[EndChain[btc]].first.pos + ExtendEndPairsMatchesLength[EndChain[btc]] << "\t"
// 								  << ExtendEndPairs[EndChain[btc]].second.pos + ExtendEndPairsMatchesLength[EndChain[btc]] << "\t"
// 								  << 0 << endl;
// 						}
// 						else if (str == 1) {
// 							Sclust << read.length - ExtendEndPairs[EndChain[btc]].first.pos - ExtendEndPairsMatchesLength[EndChain[btc]] << "\t"
// 								  << ExtendEndPairs[EndChain[btc]].second.pos + ExtendEndPairsMatchesLength[EndChain[btc]] << "\t"
// 								  << read.length - ExtendEndPairs[EndChain[btc]].first.pos << "\t"
// 								  << ExtendEndPairs[EndChain[btc]].second.pos << "\t"
// 								  << 1 << endl;					
// 						}		
// 					}
// 				Sclust.close();
// 			}
// 		}	

// 		//
// 		// Use linear alignment to ligand the gap between local SDP chain;
// 		//
// 		GenomePos btc_curReadEnd, btc_nextReadStart, btc_curGenomeEnd, btc_nextGenomeStart;
// 		btc_curReadEnd = curReadEnd; btc_curGenomeEnd = curGenomeEnd;
// 		for (int btc = EndChain.size() - 1; btc >= 0; btc--) {
// 			btc_nextGenomeStart = ExtendEndPairs[EndChain[btc]].second.pos;
// 			btc_nextReadStart = ExtendEndPairs[EndChain[btc]].first.pos;
// 			RefineByLinearAlignment(btc_curReadEnd, btc_curGenomeEnd, btc_nextReadStart, btc_nextGenomeStart, str, chromIndex, alignment, read, genome, strands, scoreMat, pathMat, opts, buff);					
// 			alignment->blocks.push_back(Block(ExtendEndPairs[EndChain[btc]].first.pos, ExtendEndPairs[EndChain[btc]].second.pos, ExtendEndPairsMatchesLength[EndChain[btc]]));
// 			debug_alignment(alignment, read.length, genome.lengths[chromIndex]);
// 			btc_curReadEnd = btc_nextReadStart + ExtendEndPairsMatchesLength[EndChain[btc]];
// 			btc_curGenomeEnd = btc_nextGenomeStart + ExtendEndPairsMatchesLength[EndChain[btc]];			
// 		}
// 		//
// 		// Add the linear alignment after the last anchor on BtwnChain;
// 		//
// 		if (block == 0) {
// 			btc_nextGenomeStart = nextGenomeStart;
// 			btc_nextReadStart = nextReadStart;			
// 			if (btc_nextGenomeStart > btc_curGenomeEnd and btc_nextReadStart > btc_curReadEnd) {
// 				RefineByLinearAlignment(btc_curReadEnd, btc_curGenomeEnd, btc_nextReadStart, btc_nextGenomeStart, str, chromIndex, 
// 										alignment, read, genome, strands, scoreMat, pathMat, opts, buff);					
// 			}				
// 		}	
// 	}
// 	else {return;}
// 	return;
// }

//
// For chaining of pure matches
//
void
LocalRefineAlignment( vector<UltimateChain> &ultimatechains, vector<Cluster> &ext_clusters, vector<SegAlignmentGroup> &alignments, const Options &smallOpts, 
	const vector<float> & LookUpTable, Read &read, char *strands[2], int &h, Genome &genome, int &LSC, const Options &tinyOpts, AffineAlignBuffers &buff, ostream *svsigstrm){

	for (int st = 0; st < ultimatechains.size(); st++) {
		if (ultimatechains[st].size() <= 1) continue;
		// Refine and store the alignment; NOTICE: when filling in the space between two adjacent anchors, 
		// the process works in forward direction, so we need to flip the small matches
		// found in the spaces before insert them into the alignment if necessary;
		// INPUT: vector<int> finalchainSeperateStrand; OUTPUT: Alignment *alignment;
		//
		int ad = alignments.back().SegAlignment.size();
		int start = 0, end = ultimatechains[st].size() - 1;
		// assert(start < end); 
		bool str = ultimatechains[st].strand(start);
		int cln = ultimatechains[st].ClusterNum(start);
		int chromIndex = ext_clusters[cln].chromIndex;
		Alignment *alignment = new Alignment(ultimatechains[st].FirstSDPValue, strands[str], read.seq, read.length, 
											 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
											 genome.header.names[chromIndex], chromIndex); 
		// alignment->NumOfAnchors0 = ultimatechains[st].NumOfAnchors0;
		alignment->NumOfAnchors1 = ultimatechains[st].NumOfAnchors1;
		alignments.back().SegAlignment.push_back(alignment);
		vector<int> scoreMat;
		vector<Arrow> pathMat;
		//
		// Set the secondary or supplymentary flag for the alignment; 
		//
		if (h > 0) alignment->ISsecondary = 1;
		if (st != LSC) alignment->Supplymentary = 1;
		GenomePos qPos = 0;
		bool inv_str = 0, inversion = 0, breakalignment = 0, prev_inv = 0;
		if (str == 0) {
			int last = end;
			int fl = end; inv_str = 1; prev_inv = end;
			while (fl > start) {
				int cur = fl;
				int next = fl - 1;
				//
				// Check the distance between two anchors. If the distance is too large, refine anchors and apply 3rd SDP;
				// Otherwise apply linear alignment. 
				//
				RefinedAlignmentbtwnAnchors(cur, next, str, inv_str, chromIndex, ultimatechains[st], alignments, alignments.back().SegAlignment.back(), 
											read, genome, strands, scoreMat, pathMat, tinyOpts, buff, 
											LookUpTable, inversion, breakalignment, svsigstrm);
				if (inversion == 1 or breakalignment == 1) { // an inversion found in between two anchors
					if (inversion == 1) alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
					else alignments.back().SegAlignment.back()->UpdateParameters(str, smallOpts, LookUpTable, svsigstrm, strands);
					alignments.back().SegAlignment.back()->NumOfAnchors0 = ultimatechains[st].NumOfAnchors0;
 					alignments.back().SegAlignment.back()->NumOfAnchors1 = last - fl;
					Alignment *next_alignment = new Alignment(ultimatechains[st].FirstSDPValue, strands[str], read.seq, read.length, 
											 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
											 genome.header.names[chromIndex], chromIndex); 
					last = fl;
					next_alignment->Supplymentary = 1;
					alignments.back().SegAlignment.push_back(next_alignment);
					prev_inv = fl;
					inversion = 0;
					breakalignment = 0;
				} 
				fl--;
			}
			alignments.back().SegAlignment.back()->NumOfAnchors0 = ultimatechains[st].NumOfAnchors0;
			alignments.back().SegAlignment.back()->NumOfAnchors1 = last - fl;
			alignments.back().SegAlignment.back()->blocks.push_back(Block(ultimatechains[st][start].first.pos, ultimatechains[st][start].second.pos, ultimatechains[st].length(start)));
		}
		else {
			int last = start;
			int fl = start; 
			inv_str = 0; prev_inv = start;
			while (fl < end) {
				int cur = fl;
				int next = fl + 1;
				RefinedAlignmentbtwnAnchors(cur, next, str, inv_str, chromIndex, ultimatechains[st], alignments, alignments.back().SegAlignment.back(), 
											read, genome, strands, scoreMat, pathMat, tinyOpts, buff, 
											LookUpTable, inversion, breakalignment, svsigstrm);
				if (inversion == 1 or breakalignment == 1) {
					if (inversion == 1) alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
					else alignments.back().SegAlignment.back()->UpdateParameters(str, smallOpts, LookUpTable, svsigstrm, strands);
					alignments.back().SegAlignment.back()->NumOfAnchors0 = ultimatechains[st].NumOfAnchors0;
 					alignments.back().SegAlignment.back()->NumOfAnchors1 = fl - last;
					Alignment *next_alignment = new Alignment(ultimatechains[st].FirstSDPValue, strands[str], read.seq, read.length, 
											 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
											 genome.header.names[chromIndex], chromIndex); 
					last = fl;
					next_alignment->Supplymentary = 1;
					alignments.back().SegAlignment.push_back(next_alignment);
					prev_inv = fl;
					inversion = 0;
					breakalignment = 0;
				} 
				fl++;
			}	
			alignments.back().SegAlignment.back()->NumOfAnchors0 = ultimatechains[st].NumOfAnchors0;
			alignments.back().SegAlignment.back()->NumOfAnchors1 = fl - last;
			alignments.back().SegAlignment.back()->blocks.push_back(Block(read.length - ultimatechains[st][end].first.pos - ultimatechains[st].length(end), ultimatechains[st][end].second.pos, ultimatechains[st].length(end)));			
		}
		// // debug checking
		// for (int bb = 1; bb < alignments.back().SegAlignment.back()->blocks.size(); bb++) {
		// 	assert(alignments.back().SegAlignment.back()->blocks[bb-1].qPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= read.length);
		// 	assert(alignments.back().SegAlignment.back()->blocks[bb-1].tPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= genome.lengths[chromIndex]);
		// 	assert(alignments.back().SegAlignment.back()->blocks[bb-1].qPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= alignments.back().SegAlignment.back()->blocks[bb].qPos);
		// 	assert(alignments.back().SegAlignment.back()->blocks[bb-1].tPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= alignments.back().SegAlignment.back()->blocks[bb].tPos);	
		// }
		//
		// Set some parameters in Alignment alignment
		//
		alignments.back().SegAlignment.back()->UpdateParameters(str, smallOpts, LookUpTable, svsigstrm, strands);
		//
		// Decide the inversion flag for each segment alignment (PAF format); seek +,-,+ or -,+,-
		//
		int js = ad + 2;
		while (js < alignments.back().SegAlignment.size()) {
			if ((alignments.back().SegAlignment[js-2]->strand == 0 and alignments.back().SegAlignment[js-1]->strand == 1
				and alignments.back().SegAlignment[js]->strand == 0) or ((alignments.back().SegAlignment[js-2]->strand == 1 
				and alignments.back().SegAlignment[js-1]->strand == 0 and alignments.back().SegAlignment[js]->strand == 1))) {
				//
				// inversion cannot be too far (<10,000) from alignments at both sides;
				//
				if (alignments.back().SegAlignment[js-1]->tStart > alignments.back().SegAlignment[js-2]->tEnd + 10000
					or alignments.back().SegAlignment[js]->tStart > alignments.back().SegAlignment[js-1]->tEnd + 10000) {
					js++;
					continue;
				}
				//
				// length requirements; 
				//
				else if (alignments.back().SegAlignment[js]->nm < 500 or alignments.back().SegAlignment[js-1]->nm < 500
						or alignments.back().SegAlignment[js-2]->nm < 40 or alignments.back().SegAlignment[js-2]->nm > 15000) {
					js++;
					continue;
				}
				else if (alignments.back().SegAlignment[js-2]->typeofaln != 3) {
					alignments.back().SegAlignment[js-1]->typeofaln=3;
				}
			}
			js++;
		}
	}
	ultimatechains.clear();
	
}

#endif
