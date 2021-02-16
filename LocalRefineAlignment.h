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
void 
debug_alignment (Alignment *alignment) {
	// This is for debugging
	if (alignment->blocks.size() > 1) {
		int last=alignment->blocks.size();
		assert(alignment->blocks[last-2].qPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].qPos);
		assert(alignment->blocks[last-2].tPos + alignment->blocks[last-2].length <= alignment->blocks[last-1].tPos);
	}	
}

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
PassgenomeThres(int cur, GenomePos & genomeThres, FinalChain & finalchain) {
	if (genomeThres != 0 and finalchain[cur].second.pos < genomeThres) return 1;
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

//
// This function removes spurious anchors after SDP;
//
void 
RemoveSpuriousAnchors(FinalChain & chain, Options & opts) {

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
void 
RemovePairedIndels (FinalChain &chain) {

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
				int Gap = (int)(((long)chain[c].second.pos - (long)chain[c].first.pos) - 
							    ((long)chain[c-1].second.pos - (long)chain[c-1].first.pos));

				if (abs(Gap) > 30) { //30
					SV.push_back(Gap);	
					SVgenome.push_back(chain[c].second.pos);
					SVpos.push_back(c);
				}
			}
			else {
				int Gap = (int)((long)(chain[c].first.pos + chain[c].second.pos) - 
								(long)(chain[c-1].first.pos + chain[c-1].second.pos));
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
void 
RemovePairedIndels (GenomePairs &matches, vector<unsigned int> &chain, vector<int> &lengths) {

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

int 
AlignSubstrings(char *qSeq, GenomePos &qStart, GenomePos &qEnd, char *tSeq, GenomePos &tStart, GenomePos &tEnd,
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
						GenomePos genomeSubEnd, vector<int> &scoreMat, vector<Arrow> &pathMat, Alignment &aln, Options &opts,
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
						bool & str, int & chromIndex, Alignment * alignment, Read & read, Genome & genome, char *strands[2], 
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
			if (opts.dotPlot and !opts.readname.empty() and read.name == opts.readname) {
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
RefinedAlignmentbtwnAnchors(int &cur, int &next, bool &str, bool &inv_str, int &chromIndex, FinalChain &finalchain, vector<SegAlignmentGroup> &alignments,
							Alignment *alignment, Read &read, Genome &genome, char *strands[2], vector<int> &scoreMat, 
							vector<Arrow> &pathMat, Options &tinyOpts, GenomePos &genomeThres, AffineAlignBuffers &buff, 
							const vector<float> & LookUpTable, bool &inversion, ostream *svsigstrm) {

	if (str == 0) alignment->blocks.push_back(Block(finalchain[cur].first.pos, finalchain[cur].second.pos, finalchain.length(cur)));
	else alignment->blocks.push_back(Block(read.length - finalchain[cur].first.pos - finalchain.length(cur), finalchain[cur].second.pos, finalchain.length(cur)));
	debug_alignment(alignment);
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
			GenomePairs *BtwnPairs;
			GenomePairs for_BtwnPairs;
			GenomePairs rev_BtwnPairs;
			//cerr << "abs(read_dist-genome_dist): "  << abs(read_dist-genome_dist)<< endl;
			// cerr << "min(read_dist, genome_dist): " << min(read_dist, genome_dist) << endl;
			// cerr << "curReadEnd: " << curReadEnd << "  curGenomeEnd: " << curGenomeEnd << "  nextReadStart: " << nextReadStart << "  nextGenomeStart: " << nextGenomeStart << endl;
			if (abs(read_dist-genome_dist) < 50) tinyOpts.refineSpaceDiag = 30;
			else tinyOpts.refineSpaceDiag = 80;
			Options nanoOpts = tinyOpts;
			if (min(read_dist, genome_dist) >= 3000 and abs(read_dist-genome_dist) >= 3000) { // big tandem repeat
				nanoOpts.globalK += 3;
			}
			RefineSpace(0, for_BtwnPairs, nanoOpts, genome, read, strands, chromIndex, nextReadStart, curReadEnd, nextGenomeStart, 
							curGenomeEnd, str);

			//
			// check larger refineSpaceDiag;
			//
			if (min(read_dist, genome_dist) > 30000 and (for_BtwnPairs.size() / (float) min(read_dist, genome_dist)) < tinyOpts.anchorstoosparse) { // this case likely means there is a large SV events. 
				//cerr << "1: " << (BtwnPairs.size() / (float) min(read_dist, genome_dist)) << endl;
				for_BtwnPairs.clear();
				nanoOpts.refineSpaceDiag = 8000;
				RefineSpace(0, for_BtwnPairs, nanoOpts, genome, read, strands, chromIndex, nextReadStart, curReadEnd, nextGenomeStart, 
							curGenomeEnd, str);		
				//cerr << "2: " << (BtwnPairs.size() / (float) min(read_dist, genome_dist)) << endl;
			} 
			//
			// try finding seeds in the inversed direction
			//
			if ((for_BtwnPairs.size() / (float) min(read_dist, genome_dist)) < tinyOpts.anchorstoosparse 
					and alignments.back().SegAlignment.back()->blocks.size() >=5 ) { // try inversion
				//BtwnPairs.clear();
				GenomePos temp = curReadEnd;
				curReadEnd = read.length - nextReadStart;
				nextReadStart = read.length - temp;
				nanoOpts.refineSpaceDiag = 100; //80
				RefineSpace(0, rev_BtwnPairs, nanoOpts, genome, read, strands, chromIndex, nextReadStart, curReadEnd, nextGenomeStart, 
							curGenomeEnd, inv_str);	
				inversion = 1;	

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
				LinearExtend(BtwnPairs, ExtendBtwnPairs, ExtendBtwnPairsMatchesLength, nanoOpts, genome, read, chromIndex);

				if (inversion == 0) {
					//
					// insert the next anchor;
					//
					ExtendBtwnPairs.push_back(GenomePair(GenomeTuple(0, nextReadStart), GenomeTuple(0, nextGenomeStart)));
					ExtendBtwnPairsMatchesLength.push_back(finalchain.length(next));
					// 
					// insert the previous anchor;
					//
					ExtendBtwnPairs.push_back(GenomePair(GenomeTuple(0, curReadEnd-finalchain.length(cur)), GenomeTuple(0, curGenomeEnd-finalchain.length(cur))));
					ExtendBtwnPairsMatchesLength.push_back(finalchain.length(cur));			
				}

				TrimOverlappedAnchors(ExtendBtwnPairs, ExtendBtwnPairsMatchesLength);
				// nanoOpts.coefficient=12;//9 this is calibrately set to 12

				float inv_value = 0; int inv_NumofAnchors = 0;
				SparseDP_ForwardOnly(ExtendBtwnPairs, ExtendBtwnPairsMatchesLength, BtwnChain, nanoOpts, LookUpTable, inv_value, inv_NumofAnchors, 1); //6
				RemovePairedIndels(ExtendBtwnPairs, BtwnChain, ExtendBtwnPairsMatchesLength);
				//
				// Use linear alignment to ligand the gap between local SDP chain;
				//
				GenomePos btc_curReadEnd, btc_nextReadStart, btc_curGenomeEnd, btc_nextGenomeStart;
				btc_curReadEnd = curReadEnd;
				btc_curGenomeEnd = curGenomeEnd;

				if (nanoOpts.dotPlot and !nanoOpts.readname.empty() and read.name == nanoOpts.readname) {
					ofstream pSclust("BtwnPairs.tab", std::ofstream::app);
					for (int bp = BtwnPairs->size()-1; bp >= 0; bp--) {
						if (inversion == 0) {
							if (str == 0) {
								pSclust << (*BtwnPairs)[bp].first.pos << "\t"
									  << (*BtwnPairs)[bp].second.pos << "\t"
									  << (*BtwnPairs)[bp].first.pos + nanoOpts.globalK << "\t"
									  << (*BtwnPairs)[bp].second.pos + nanoOpts.globalK << "\t"
									  << 0 << endl;
							}
							else if (str == 1) {
								pSclust << read.length - (*BtwnPairs)[bp].first.pos - nanoOpts.globalK << "\t"
									  << (*BtwnPairs)[bp].second.pos + nanoOpts.globalK << "\t"
									  << read.length - (*BtwnPairs)[bp].first.pos << "\t"
									  << (*BtwnPairs)[bp].second.pos << "\t"
									  << 1 << endl;					
							}						
						} 
						else {
							if (inv_str == 0) {
								pSclust << (*BtwnPairs)[bp].first.pos << "\t"
									  << (*BtwnPairs)[bp].second.pos << "\t"
									  << (*BtwnPairs)[bp].first.pos + nanoOpts.globalK << "\t"
									  << (*BtwnPairs)[bp].second.pos + nanoOpts.globalK << "\t"
									  << 0 << endl;
							}
							else if (inv_str == 1) {
								pSclust << read.length - (*BtwnPairs)[bp].first.pos - nanoOpts.globalK << "\t"
									  << (*BtwnPairs)[bp].second.pos + nanoOpts.globalK << "\t"
									  << read.length - (*BtwnPairs)[bp].first.pos << "\t"
									  << (*BtwnPairs)[bp].second.pos << "\t"
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
									  << ExtendBtwnPairs[bp].second.pos << "\t"
									  << ExtendBtwnPairs[bp].first.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << 0 << endl;
							}
							else if (str == 1) {
								eSclust << read.length - ExtendBtwnPairs[bp].first.pos - ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << read.length - ExtendBtwnPairs[bp].first.pos << "\t"
									  << ExtendBtwnPairs[bp].second.pos << "\t"
									  << 1 << endl;					
							}							
						}
						else {
							if (inv_str == 0) {
								eSclust << ExtendBtwnPairs[bp].first.pos << "\t"
									  << ExtendBtwnPairs[bp].second.pos << "\t"
									  << ExtendBtwnPairs[bp].first.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << 0 << endl;
							}
							else if (inv_str == 1) {
								eSclust << read.length - ExtendBtwnPairs[bp].first.pos - ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << ExtendBtwnPairs[bp].second.pos + ExtendBtwnPairsMatchesLength[bp] << "\t"
									  << read.length - ExtendBtwnPairs[bp].first.pos << "\t"
									  << ExtendBtwnPairs[bp].second.pos << "\t"
									  << 1 << endl;					
							}							
						}
					}
					eSclust.close();
				}	

				if (nanoOpts.dotPlot and !nanoOpts.readname.empty() and read.name == nanoOpts.readname) {
					ofstream Sclust("SparseDP_Forward.tab", std::ofstream::app);
						for (int btc = BtwnChain.size()-1; btc >= 0; btc--) {
							if (inversion == 0) {
								if (str == 0) {
									Sclust << ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].first.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << 0 << endl;
								}
								else if (str == 1) {
									Sclust << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos - ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos << "\t"
										  << 1 << endl;					
								}		
							}
							else {
								if (inv_str == 0) {
									Sclust << ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].first.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << 0 << endl;
								}
								else if (inv_str == 1) {
									Sclust << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos - ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos + ExtendBtwnPairsMatchesLength[BtwnChain[btc]] << "\t"
										  << read.length - ExtendBtwnPairs[BtwnChain[btc]].first.pos << "\t"
										  << ExtendBtwnPairs[BtwnChain[btc]].second.pos << "\t"
										  << 1 << endl;					
								}							
							}

						}
					Sclust.close();
				}
				int btc_end = BtwnChain.size()-1; int btc_start = 0;
				if (BtwnChain.back() == ExtendBtwnPairs.size()-1) btc_end = BtwnChain.size()-2;
				if (BtwnChain[0] == ExtendBtwnPairs.size()-2) btc_start = 1;

				if (inversion == 1) {
					alignment->UpdateParameters(str, tinyOpts, LookUpTable, svsigstrm, strands);
					Alignment *inv_alignment = new Alignment(inv_value, strands[str], read.seq, read.length, read.name, inv_str, read.qual, genome.seqs[chromIndex],  
																		 genome.lengths[chromIndex], genome.header.names[chromIndex], chromIndex, 0); 
					alignment->NumOfAnchors1 = BtwnChain.size();
					alignments.back().SegAlignment.push_back(inv_alignment);	
				}
				
				for (int btc = btc_end; btc >= btc_start; btc--) {
					btc_nextGenomeStart = ExtendBtwnPairs[BtwnChain[btc]].second.pos;
					btc_nextReadStart = ExtendBtwnPairs[BtwnChain[btc]].first.pos;
					RefineByLinearAlignment(btc_curReadEnd, btc_curGenomeEnd, btc_nextReadStart, btc_nextGenomeStart, str, chromIndex, 
											alignments.back().SegAlignment.back(), read, genome, strands, scoreMat, pathMat, nanoOpts, genomeThres, buff);				

					alignments.back().SegAlignment.back()->blocks.push_back(Block(ExtendBtwnPairs[BtwnChain[btc]].first.pos, ExtendBtwnPairs[BtwnChain[btc]].second.pos, 
														ExtendBtwnPairsMatchesLength[BtwnChain[btc]]));
					debug_alignment(alignments.back().SegAlignment.back());
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
											alignments.back().SegAlignment.back(), read, genome, strands, scoreMat, pathMat, nanoOpts, genomeThres, buff);					
				}				
			}
			else {
				// inversion = 0;
				// GenomePos temp = curReadEnd;
				// curReadEnd = read.length - nextReadStart;
				// nextReadStart = read.length - temp;
				RefineByLinearAlignment(curReadEnd, curGenomeEnd, nextReadStart, nextGenomeStart, str, chromIndex, 
									alignment, read, genome, strands, scoreMat, pathMat, tinyOpts, genomeThres, buff);
			}
		}
		else {
			RefineByLinearAlignment(curReadEnd, curGenomeEnd, nextReadStart, nextGenomeStart, str, chromIndex, 
									alignment, read, genome, strands, scoreMat, pathMat, tinyOpts, genomeThres, buff);				
		}
	}
	else genomeThres = curGenomeEnd;
}

void
SparseDP_and_RefineAlignment_btwn_anchors(vector<Primary_chain> &Primary_chains, vector<SplitChain> &splitchains, vector<Cluster> &ExtendClusters, 
		vector<SegAlignmentGroup> &alignments, Options &smallOpts, const vector<float> & LookUpTable, Read &read, char *strands[2], int &p, int &h, 
		Genome &genome, int &LSC, Options &tinyOpts, AffineAlignBuffers &buff, ostream *svsigstrm){

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
		// RemovePairedIndels(finalchain); (This is deleting lots of good anchors)
		// RemoveSpuriousAnchors(finalchain, smallOpts);

		//cerr << "2nd SDP done!" << endl;
		if (finalchain.size() == 0) continue; // cannot be mapped to the genome!
		if (smallOpts.dotPlot and !smallOpts.readname.empty() and read.name == smallOpts.readname) {
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
		int ad = alignments.back().SegAlignment.size();
		for (int fsc = 0; fsc < finalSeperateChain.size(); fsc++) {

			assert(finalSeperateChain[fsc].size() == 2);

			int start = finalSeperateChain[fsc][0];
			int end = finalSeperateChain[fsc][1];
			assert(start < end); 
			bool str = finalchain.strand(start);
			int cln = finalchain.ClusterNum(start);
			int chromIndex = ExtendClusters[cln].chromIndex;	
			Alignment *alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, 
												 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
												 genome.header.names[chromIndex], chromIndex, 0); 
			alignment->NumOfAnchors0 = Primary_chains[p].chains[h].NumOfAnchors0;
			alignment->NumOfAnchors1 = finalchain.size();
			alignment->SecondSDPValue = finalchain.SecondSDPValue;
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
			bool inv_str = 0, inversion = 0, prev_inv = 0;
			if (str == 0) {
				int fl = end - 1;
				inv_str = 1; prev_inv = end - 1;
				while (fl > start) {
					int cur = fl;
					int next = fl - 1;
					if (!PassgenomeThres(cur, genomeThres, finalchain)) {
						//
						// Check the distance between two anchors. If the distance is too large, refine anchors and apply 3rd SDP;
						// Otherwise apply linear alignment. 
						//
						RefinedAlignmentbtwnAnchors(cur, next, str, inv_str, chromIndex, finalchain, alignments, alignments.back().SegAlignment.back(), 
													read, genome, strands, scoreMat, pathMat, tinyOpts, genomeThres, buff, 
													LookUpTable, inversion, svsigstrm);

						if (inversion == 1) {
							//alignment.FirstSDPValue /= (float)(prev_inv-fl) / (end-start)
							//alignment->NumOfAnchors = prev_inv-fl;
							alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
							Alignment *next_alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, 
													 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
													 genome.header.names[chromIndex], chromIndex, 0); 
							next_alignment->NumOfAnchors1 = fl-start;
							alignments.back().SegAlignment.push_back(next_alignment);
							prev_inv = fl;
							inversion = 0;
						} 
					}
					fl--;
				}
				if (!PassgenomeThres(start, genomeThres, finalchain)) {
					alignments.back().SegAlignment.back()->blocks.push_back(Block(finalchain[start].first.pos, finalchain[start].second.pos, 
														finalchain.length(start)));
				}
			}
			else {
				int fl = start; 
				inv_str = 0; prev_inv = start;
				while (fl < end - 1) {
					int cur = fl;
					int next = fl + 1;
					if (!PassgenomeThres(cur, genomeThres, finalchain)) {
						RefinedAlignmentbtwnAnchors(cur, next, str, inv_str, chromIndex, finalchain, alignments, alignments.back().SegAlignment.back(), 
													read, genome, strands, scoreMat, pathMat, tinyOpts, genomeThres, buff, 
													LookUpTable, inversion, svsigstrm);
						if (inversion == 1) {
							//alignment.FirstSDPValue /= (float)(prev_inv-fl) / (end-start)
							alignments.back().SegAlignment.back()->UpdateParameters(inv_str, smallOpts, LookUpTable, svsigstrm, strands);
							//alignment->NumOfAnchors = fl-start;
							Alignment *next_alignment = new Alignment(Primary_chains[p].chains[h].value, strands[str], read.seq, read.length, 
													 read.name, str, read.qual, genome.seqs[chromIndex], genome.lengths[chromIndex], 
													 genome.header.names[chromIndex], chromIndex, 0); 
							next_alignment->NumOfAnchors1 = end-fl;
							alignments.back().SegAlignment.push_back(next_alignment);
							prev_inv = fl;
							inversion = 0;
						} 
					}
					fl++;
				}	
				if (!PassgenomeThres(end - 1, genomeThres, finalchain)) {
					alignments.back().SegAlignment.back()->blocks.push_back(Block(read.length-finalchain[end - 1].first.pos-finalchain.length(end - 1), 
														finalchain[end - 1].second.pos, finalchain.length(end - 1)));
				}
			}
		
			for (int bb = 1; bb < alignments.back().SegAlignment.back()->blocks.size(); bb++) {
				assert(alignments.back().SegAlignment.back()->blocks[bb-1].qPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= alignments.back().SegAlignment.back()->blocks[bb].qPos);
				assert(alignments.back().SegAlignment.back()->blocks[bb-1].tPos + alignments.back().SegAlignment.back()->blocks[bb-1].length <= alignments.back().SegAlignment.back()->blocks[bb].tPos);	
			}
			//
			// Set some parameters in Alignment alignment
			//
			alignments.back().SegAlignment.back()->UpdateParameters(str, smallOpts, LookUpTable, svsigstrm, strands);
		}
		//
		// Decide the inversion flag for each segment alignment; seek +,-,+ or -,+,-
		//
		int js = ad+2;
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


#endif