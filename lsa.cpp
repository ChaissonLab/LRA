#include "htslib/hts.h"
#include "htslib/kseq.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <vector>
#include "htslib/kseq.h"
#include <algorithm>
#include <queue>
#include "SeqUtils.h"
#include "TupleOps.h"
#include "MinCount.h"
#include "MMIndex.h"
#include "CompareLists.h"
#include "Sorting.h"
#include <algorithm>
#include "Options.h"
#include "Clustering.h"
#include "seqan/seeds.h"
#include "seqan/align.h"
#include "Path.h"
#include "KBandAlign.h"
#include "GlobalChain.h"

// Print results to stdout.
typedef seqan::String<seqan::Dna> TSequence;                 // sequence type
//typedef seqan::Infix<seqan::String<seqan::Dna> >::Type  Substring;
typedef seqan::Infix<char* >::Type  Substring;

typedef seqan::Align<Substring, seqan::ArrayGaps> TAlign;     // align type
typedef seqan::Row<TAlign>::Type TRow;                 // gapped sequence typ

typedef seqan::Align<TSequence, seqan::ArrayGaps> TSeqAlign;   
typedef seqan::Row<TSeqAlign>::Type TSeqRow;                 // gapped sequence typ


bool ArgIs(const char* a, const char* b) {
	return strcmp(a,b) == 0;
}

void HelpMap() {
	cout << "Usage: lsc align genome.fa reads.fa [options]" << endl;
	cout << "Options:" << endl
			 << "  -p  (flag)  View pairwise alignment." << endl;
}

void SwapStrand(kseq_t* read, Options &opts, GenomePairs &matches) {
	for (int m=0; m < matches.size(); m++) {
		matches[m].first.pos = read->seq.l - (matches[m].first.pos + opts.k - 1);
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
		
int SetStrand(kseq_t *read, Genome &genome, Options &opts, GenomePairs &matches) { 
	int nSame=0;
	int nDifferent=0;
	for (int m=0; m< matches.size(); m++) {
		int chromIndex = genome.header.Find(matches[m].second.pos);
		char *chrom=genome.seqs[chromIndex];
		int chromPos = matches[m].second.pos - genome.header.pos[chromIndex];
		GenomeTuple readTup, genomeTup;
		StoreTuple(read->seq.s, matches[m].first.pos, opts.k, readTup);
		StoreTuple(chrom, chromPos, opts.k, genomeTup);
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

void ReverseClusterStrand(kseq_t* read, Genome &genome, Options &opts, 
											vector<Cluster> &clusters) {
	for (int c = 0; c < clusters.size(); c++) {
			SwapStrand(read, opts, clusters[c].matches);
			clusters[c].strand = 1;
	}
}


void SetClusterStrand(kseq_t* read, Genome &genome, Options &opts, 
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
		if (clusters[c].start == clusters[c].end or clusters[c].matches.size() < minSize ) {
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
	RemoveEmptyClusters(order.clusters);
}

typedef seqan::Iterator<TSeqRow>::Type TRowIterator;
typedef char TChar;                             // character type

int AlignSubstrings(char *qSeq, GenomePos &qStart, GenomePos &qEnd, 
										char *tSeq, GenomePos &tStart, GenomePos &tEnd, 
										vector<int> &scoreMat, vector<Arrow> &pathMat, Alignment &aln) {
	
	int qLen = qEnd-qStart;
	int tLen = tEnd-tStart;
	int drift = abs(qLen - tLen);
	int k = 7;
	k = min(max(k,drift+4),30);
	/*
	int score = KBandAlign(&qSeq[qStart], qEnd-qStart, &tSeq[tStart], tEnd-tStart, 
												 -5,3,2,2, k, // make these smart later.
												 scoreMat, pathMat, aln);*/

	string readSeq(&qSeq[qStart], qEnd-qStart);
	string chromSeq(&tSeq[tStart],tEnd-tStart);
	TSeqAlign align;
	seqan::resize(seqan::rows(align), 2);
	seqan::assignSource(seqan::row(align, 0), readSeq.c_str());
	seqan::assignSource(seqan::row(align, 1), chromSeq.c_str());
	int score = seqan::globalAlignment(align, seqan::Score<int, seqan::Simple>(4, -4, -20, -0), -k, k);

	TSeqRow & row1 = seqan::row(align, 0);
	TSeqRow & row2 = seqan::row(align, 1);
	int q=0, t=0;
	int rowLen = seqan::length(align);
	int i=0;


	TRowIterator qit = begin(row1),
		qitEnd = end(row1),
		tit = begin(row2);
	while ( qit != qitEnd ) {
		while(qit  < qitEnd and (seqan::isGap(qit) == true or seqan::isGap(tit) == true)) {
			/*
			TChar qc = isGap(qit) ? '-' : readSeq[q];
			TChar tc = isGap(tit) ? '-' : chromSeq[t];
			cout << qc << " " << tc << endl;
			*/
			if (seqan::isGap(qit)) { t+=1;}
			if (seqan::isGap(tit)) { q+=1;}
			qit++;
			tit++;
		}
		if (qit != qitEnd and seqan::isGap(qit) == false and seqan::isGap(tit) == false) {
			int bqs = q;
			int bts = t;
			int is=i;
			while(qit != qitEnd and seqan::isGap(qit) == false and seqan::isGap(tit) == false) {
					qit++;
					tit++;
					q++;
					t++;
			}
			aln.blocks.push_back(Block(bqs, bts, q-bqs));
		}
	}
	return score;
}
/*
void AddGaps(int p, int g, TRow &r){ 
	if (g > 0) {
		seqan::insertGaps(r, p, g);
	}
}
*/


void RefineSubstrings(char *read,   GenomePos readSubStart, GenomePos readSubEnd, 
											char *genome, GenomePos genomeSubStart, GenomePos genomeSubEnd, 
											vector<int> &scoreMat, vector<Arrow> &pathMat, 
											Alignment &aln) {

	aln.blocks.clear();
	AlignSubstrings( read, readSubStart, readSubEnd, genome, genomeSubStart, genomeSubEnd, scoreMat, pathMat, aln);
	/*
	string q,a,t;
	aln.CreateAlignmentStrings(&read[readSubStart], &genome[genomeSubStart], q, a, t);
	cout << "refined " << endl;
	cout << q << endl << a << endl << t << endl;
	*/
	for (int b = 0; b < aln.blocks.size(); b++) {
		aln.blocks[b].qPos += readSubStart;
		aln.blocks[b].tPos += genomeSubStart;
	}
	
}


int Matched(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te) {
	return min(qe-qs, te-ts);
}

void SetMatchAndGaps(GenomePos qs, GenomePos qe, GenomePos ts, GenomePos te, int &m, int &qg, int &tg) {
	m=Matched(qs, qe, ts, te);
	qg=qe-qs-m;
	tg=te-ts-m;
}

struct IndexedSeed_;
typedef seqan::Tag<IndexedSeed_> IndexedSeed;;

template<typename TConfig>
class seqan::Seed<IndexedSeed, TConfig> : public seqan::Seed<Simple, TConfig>{
public:
	int index;
	Seed(int i, int j, int k, int l) : seqan::Seed<Simple,TConfig>(i,j,k,l){index=-1;}
	Seed(int i, int j, int k) : seqan::Seed<Simple,TConfig>(i,j,k){index=-1;}
	Seed() : seqan::Seed<Simple,TConfig>(){index=-1;}	
	Seed(int i, int j, int k, int l, int idx) : seqan::Seed<Simple,TConfig>(i,j,k,l), index(idx) {}
};

void RunAlign(int argc, char* argv[], Options &opts ) {
	// open query file for reading; you may use your favorite FASTA/Q parser
	int argi = 0;
	string genomeFile = "", reads = "";
	string indexFile="";
	int w=10;
	bool storeAll = false;
	for (argi = 0; argi < argc; ) {
		if (ArgIs(argv[argi], "-i")) {
			++argi;
			indexFile=argv[argi];
		}		
		else if (ArgIs(argv[argi], "-a")) {
			++argi;
			storeAll=true;
		}		
		else if (ArgIs(argv[argi], "-w")) {
			++argi;
			opts.w=atoi(argv[argi]);
		}		
		else if (ArgIs(argv[argi], "-M")) {
			++argi;
			opts.minClusterSize= atoi(argv[argi]);
		}		
		else if (ArgIs(argv[argi], "-m")) {
			++argi;
			opts.minRefinedClusterSize= atoi(argv[argi]);
		}		
		else if (ArgIs(argv[argi], "-f")) {
			++argi;
			opts.maxFreq = atoi(argv[argi]);
		}		
		else if (ArgIs(argv[argi], "-k")) {
			++argi;
			opts.k=atoi(argv[argi]);
		}		
		else if (ArgIs(argv[argi], "-p")) {
			opts.viewPairwise=true;
		}
		else {
			if (genomeFile == "") {
				genomeFile = argv[argi];
			}
			else if (reads == "") {
				reads = argv[argi];
			}
		}
		++argi;
	}

	if (genomeFile == "" || reads == "") {
		HelpMap();
		exit(1);
	}
	if (indexFile == "") {
		indexFile = genomeFile + ".mmi";
	}
	Header header;
	vector<GenomeTuple> genomemm;
	LocalIndex glIndex;
	if (ReadIndex(indexFile, genomemm, header, opts) == 0) {
		StoreIndex(genomeFile, genomemm, header, opts);
	}
	if (glIndex.Read(genomeFile+".gli") == 0) {
		glIndex.IndexFile(genomeFile);
	}
	GenomePos mm=0;
	for(GenomePos mi =0; mi < genomemm.size(); mi++) {
		if (genomemm[mi].pos > mm) {
			mm = genomemm[mi].pos;
		}
	}
	Genome genome;
	genome.Read(genomeFile);

	
	gzFile f = gzopen(reads.c_str(), "r");
	kseq_t *ks = kseq_init(f);
	int offset=0;
	vector<GenomeTuple> readmm;
	vector<pair<GenomeTuple, GenomeTuple> > matches;
	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		readmm.clear();
		matches.clear();
		if (storeAll) {
			Options allOpts = opts;
			allOpts.w=1;
			StoreMinimizers<GenomeTuple, Tuple>(ks->seq.s, ks->seq.l,
																					allOpts.k, allOpts.w, readmm);			
		}
		else {
			StoreMinimizers<GenomeTuple, Tuple>(ks->seq.s, ks->seq.l, opts.k, opts.w, readmm);
		}
		sort(readmm.begin(), readmm.end());
		CompareLists(readmm, genomemm, matches, opts);
		/*
			stringstream allm;
			allm << "allm.tab";
			ofstream callm(allm.str().c_str());
			for (int m=0; m < matches.size(); m++) {
				callm << matches[m].second.pos << "\t" << matches[m].first.pos << endl;
			}
			callm.close();
		*/
#ifdef _TESTING_				
				PrintPairs(matches, opts.k);
#endif
		
		DiagonalSort<GenomeTuple>(matches);

		CleanOffDiagonal(matches, opts);

		vector<Cluster> clusters;

		StoreDiagonalClusters(matches, clusters, opts);
		//cerr << "Warning, need to add back diagonal clusters." << endl;
		AntiDiagonalSort<GenomeTuple>(matches, genome.GetSize());

		SwapStrand(ks, opts, matches);
	
		vector<Cluster> revClusters;
		StoreDiagonalClusters(matches, revClusters, opts, 1);
 
		SwapStrand(ks, opts, matches);

		//
		// Add pointers to seq that make code more readable.
		//
		char *readRC;
		char *read=ks->seq.s;

		int readLen = ks->seq.l;
		CreateRC(read, readLen, readRC);
		char *strands[2] = { read, readRC };

		clusters.insert(clusters.end(), revClusters.begin(), revClusters.end());

		//
		// Build local index for refining alignments.
		//
		LocalIndex forwardIndex(glIndex);
		LocalIndex reverseIndex(glIndex);

		LocalIndex *localIndexes[2] = {&forwardIndex, &reverseIndex};
		forwardIndex.IndexSeq(read, readLen);
		reverseIndex.IndexSeq(readRC, readLen);

		Options smallOpts = opts;
		smallOpts.k=glIndex.k;
		smallOpts.w=glIndex.w;
		smallOpts.maxFreq=100;
		smallOpts.maxDiag=25;
		smallOpts.minDiagCluster=3;

		Options tinyOpts = opts;
		tinyOpts.maxFreq=3;
		tinyOpts.maxDiag=5;
		tinyOpts.minDiagCluster=2;

		//
		// Merge overlapping clusters
		//
		RemoveEmptyClusters(clusters, opts.minClusterSize);
		if (opts.mergeGapped) {
			ClusterOrder clusterOrder(clusters);
			clusterOrder.Sort();
			MergeOverlappingClusters(clusterOrder);
		}

		vector<Cluster> refinedClusters(clusters.size());
		vector<Fragment> fragments;
		vector<int>      optFragmentIndices;
		vector<Endpoint> endpoints;
		vector<int> clusterStart, clusterEnd;

		for (int c = 0; c < clusters.size(); c++) {


			if (clusters[c].start == clusters[c].end) {
				continue;
			}			
			//
			// Get the boundaries of the cluster in both sequences.
			//
			CartesianTargetSort<GenomeTuple>(clusters[c].matches.begin(), 
																			 clusters[c].matches.end());
			int nMatch = clusters[c].matches.size();
			GenomePos tPos=clusters[c].matches[0].second.pos;
			int firstChromIndex   = genome.header.Find(tPos);
			int lastChromIndex;
			if (nMatch > 1 ) {
				tPos = clusters[c].matches[nMatch-1].second.pos;
				lastChromIndex = genome.header.Find(tPos);
			} else { lastChromIndex = firstChromIndex; }
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
			/*
			stringstream name;
			cout << "cluster " << c << " strand " << (int) clusters[c].strand << endl;
			name << "cluster." << c << ".tab";
			ofstream clust(name.str().c_str());
			for (int m=0; m < clusters[c].matches.size(); m++) {
				clust << clusters[c].matches[m].second.pos << "\t" << clusters[c].matches[m].first.pos << endl;
			}
			clust.close();
			*/
			//
			// Get shorthand access to alignment boundaries.
			//

		  GenomePos readClusterStart, readClusterEnd, genomeClusterStart, genomeClusterEnd;
			readClusterStart = clusters[c].matches[0].first.pos;
			genomeClusterStart = clusters[c].matches[0].second.pos + chromOffset;

			int cl = clusters[c].size();
			readClusterEnd = clusters[c].matches[cl-1].first.pos + opts.k;
			genomeClusterEnd = clusters[c].matches[cl-1].second.pos + opts.k + chromOffset;

			/*
			vector<GenomeTuple> readClustMM, genomeClustMM;

			StoreMinimizers<GenomeTuple, Tuple>(&ks->seq.s[readClusterStart], readClusterEnd-readClusterStart,
																					glIndex.k, glIndex.w, readClustMM); 
			sort(readClustMM.begin(), readClustMM.end());

			StoreMinimizers<GenomeTuple, Tuple>(&genome.seqs[firstChromIndex][genomeClusterStart-chromOffset],
																					genomeClusterEnd - genomeClusterStart,
																					glIndex.k, glIndex.w, genomeClustMM); 
			sort(genomeClustMM.begin(), genomeClustMM.end());			
			vector<pair<GenomeTuple, GenomeTuple> > testMatches;

			CompareLists<GenomeTuple>(readClustMM.begin(), readClustMM.end(),
															 genomeClustMM.begin(), genomeClustMM.end(),
															 testMatches, smallOpts);
			*/

			int ls, le;

			GenomePos chromEndOffset   = header.GetNextOffset(genomeClusterEnd);
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
			/*
			stringstream locname;
			locname << "loc." << c << ".tab";
			ofstream loc(locname.str().c_str());
			*/

			for (int lsi=ls; lsi <= le; lsi++) {
				//
				// Find the coordinates in the cluster that start in this local index.
				//
				GenomePos genomeLocalIndexStart = 
					glIndex.seqOffsets[lsi]  - 
					chromOffset;
				GenomePos genomeLocalIndexEnd   = 
					glIndex.seqOffsets[lsi+1] - 1 - chromOffset;

				int matchStart = 
					CartesianTargetLowerBound<GenomeTuple>(clusters[c].matches.begin(),
																								 clusters[c].matches.end(),
																								 genomeLocalIndexStart);

				int matchEnd   = 
					CartesianTargetUpperBound<GenomeTuple>(clusters[c].matches.begin(),
																								 clusters[c].matches.end(),
																								 genomeLocalIndexEnd);

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
					readEnd = clusters[c].matches[matchStart].first.pos + opts.k;
				}
				//
				// Expand boundaries of read to match.
				if (lsi == le) {
					if (readEnd + opts.window > ks->seq.l) {
						readEnd = ks->seq.l; 
					}
					else { 
						readEnd += opts.window;	
					}
				}			
				
				//
				// Find the boundaries where in the query the matches should be added.
				//
				
				int queryIndexStart = readIndex->LookupIndex(readStart);
				int queryIndexEnd   = 
					readIndex->LookupIndex(min(readEnd, 
																		 (GenomePos) ks->seq.l-1));
				assert(queryIndexEnd < readIndex->seqOffsets.size()+1);

				for (int qi = queryIndexStart; qi <= queryIndexEnd; ++qi){ 
					LocalPairs smallMatches;

					GenomePos qStartBoundary = readIndex->tupleBoundaries[qi];
					GenomePos qEndBoundary   = readIndex->tupleBoundaries[qi+1];
					GenomePos readSegmentStart= readIndex->seqOffsets[qi];
					GenomePos readSegmentEnd  = readIndex->seqOffsets[qi+1];
					//					cout << lsi << " " << qi << " Comparing read mm " << qStartBoundary << " " << qEndBoundary 
					//							 << " and ref " << glIndex.tupleBoundaries[lsi] << " " << glIndex.tupleBoundaries[lsi+1] << endl;

					
					//					loc << readSegmentStart << "\t" << genomeLocalIndexStart << "\t" << readIndex->seqOffsets[qi+1] << "\t" << genomeLocalIndexEnd << endl;
	
					CompareLists<LocalTuple>(readIndex->minimizers.begin() + 
																	 qStartBoundary,

																	 readIndex->minimizers.begin() +
																	 qEndBoundary,

																	 glIndex.minimizers.begin()+
																	 glIndex.tupleBoundaries[lsi], 

																	 glIndex.minimizers.begin()+
																	 glIndex.tupleBoundaries[lsi+1], 

																	 smallMatches, smallOpts);

					lmIndex+=smallMatches.size();


					//
					// Do local processing of matches to ensure the region that is searched returns reasonable anchors.
					//

					DiagonalSort<LocalTuple>(smallMatches);
					CleanOffDiagonal(smallMatches, smallOpts);					

					AppendValues<LocalPairs>(refinedClusters[c].matches, 
																	 smallMatches.begin(), smallMatches.end(), 
																	 readSegmentStart, genomeLocalIndexStart );
				}
				
			}
			//			loc.close();
			refinedClusters[c].SetClusterBoundariesFromMatches();
			refinedClusters[c].strand = clusters[c].strand;
			refinedClusters[c].chromIndex = clusters[c].chromIndex;
		}
		RemoveEmptyClusters(refinedClusters);

		for (int r =0; r < refinedClusters.size(); r++) {

			if (refinedClusters[r].matches.size() < opts.minRefinedClusterSize) {
				continue;
			}



			/*
			stringstream ref;
			ref << "ref_cluster." << r << ".tab";
			ofstream rclust(ref.str().c_str());
			for (int m=0; m < refinedClusters[r].matches.size(); m++) {
				rclust << refinedClusters[r].matches[m].second.pos << "\t" << refinedClusters[r].matches[m].first.pos << endl;
			}
			rclust.close();
			*/
			//
			// Clean local matches to reduce chaining burden.
			//
			DiagonalSort<GenomeTuple>(refinedClusters[r].matches);
			CleanOffDiagonal(refinedClusters[r].matches, smallOpts);

			if (refinedClusters[r].matches.size() == 0) {
				continue;
			}

			
			// Build SeedSet.
			seqan::SeedSet<seqan::Seed<seqan::Simple>, seqan::Unordered> seedSet;

			for (int m=0; m< refinedClusters[r].matches.size(); m++) {
			  seqan::addSeed(seedSet, 
					 seqan::Seed<seqan::Simple>(refinedClusters[r].matches[m].second.pos, 
								    refinedClusters[r].matches[m].first.pos, 
								    glIndex.k),
					 seqan::Single());
			}

			// Perform sparse chaining, uses time O(n log n).
			seqan::String<seqan::Seed<seqan::Simple> > chain;
			seqan::chainSeedsGlobally(chain, seedSet, seqan::SparseChaining());
			
			vector<GenomePair> tupChain;
			int qPrev=0, tPrev=0;
			for (int ch=0; ch < seqan::length(chain); ch++) {

				tupChain.push_back(GenomePair(GenomeTuple(0, beginPositionV(chain[ch])),
																			GenomeTuple(0, beginPositionH(chain[ch]))));
			}
			vector<Cluster> chainClust;
			Options diagOpts;
			diagOpts = smallOpts;
			diagOpts.maxDiag=5;
			diagOpts.minClusterSize=1;
			StoreDiagonalClustersLite(tupChain, chainClust, diagOpts, refinedClusters[r].strand);
			//			sort(smallClusters.begin(), smallClusters.end());
			/*
			if (seqan::length(chain)> 0) {
				cout << "Removing paired indels for refined cluster " << r << " " << beginPositionV(chain[0]) << " " << beginPositionH(chain[0]) << endl;
				}
			*/
			RemovePairedIndels(tupChain, chainClust, smallOpts);
			
			seqan::clear(chain);
			
			for (int m=0; m< tupChain.size(); m++) {
			  seqan::append(chain, 
											 seqan::Seed<seqan::Simple>(tupChain[m].second.pos,
																									tupChain[m].first.pos, glIndex.k));

				/*				cout << m << "\t" << tupChain[m].first.pos << "\t" 
						 << tupChain[m].first.pos - qPrev << "\t" 
						 << tupChain[m].second.pos - tPrev << endl;
				qPrev =  tupChain[m].first.pos ;
				tPrev =  tupChain[m].second.pos;*/
			}
			
			/*
			stringstream refd;
			refd << "diag_cluster." << r << ".tab";
			ofstream rdclust(refd.str().c_str());
			for (int m=0; m < tupChain.size(); m++) {
				rdclust << tupChain[m].second.pos << "\t" << tupChain[m].first.pos << "\t" << glIndex.k << endl;
			}
			rdclust.close();
			*/
			GenomePos chainGenomeStart = seqan::beginPositionH(chain[0]);
			GenomePos chainGenomeEnd   = seqan::endPositionH(chain[seqan::length(chain)-1]);

			GenomePos chainReadStart = seqan::beginPositionV(chain[0]);
			GenomePos chainReadEnd   = seqan::endPositionV(chain[seqan::length(chain)-1]);


			GenomePos globalStart = refinedClusters[r].tStart;
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
			gapOpts.maxFreq=5;
			gapOpts.k=7;
			vector< seqan::String<seqan::Seed<seqan::Simple> > > refinedChains(seqan::length(chain)-1);
			// Build SeedSet.
			seqan::SeedSet<seqan::Seed<seqan::Simple>, seqan::Unordered> gapSeedSet;


			vector<int> scoreMat;
			vector<Arrow> pathMat;
			int gapK=gapOpts.k-2;
			int chainLength = seqan::length(chain);
			for (int c = 0; chainLength > 0 and c < seqan::length(chain)-1; c++) {
				GenomePos curGenomeEnd = seqan::endPositionH(chain[c]);
				GenomePos nextGenomeStart = seqan::beginPositionH(chain[c+1]);

				GenomePos curReadEnd = seqan::endPositionV(chain[c]);
				GenomePos nextReadStart = seqan::beginPositionV(chain[c+1]);
				
				int rg=nextReadStart-curReadEnd, gg=nextGenomeStart-curGenomeEnd;
				int mg=min(rg, gg);
				int rm=rg-mg;
				int gm=gg-mg;
				rg-=mg;
				gg-=mg;

				GenomePos subreadLength = nextReadStart-curReadEnd;
				GenomePos subgenomeLength = nextGenomeStart-curGenomeEnd;

				if (nextReadStart > curReadEnd and nextGenomeStart > curGenomeEnd) {

					if ((subreadLength > 30 or subgenomeLength > 30) and subreadLength < 4000 and subgenomeLength < 4000) {

						gapGenomeTup.clear();
						gapReadTup.clear();
						gapPairs.clear();

						//
						// Find matches between read and reference in the coordinate space of read and chromosome
						//
						StoreMinimizers<GenomeTuple, Tuple>( genome.seqs[chromIndex] + curGenomeEnd,
																								 subgenomeLength, gapK, 1, gapGenomeTup, false);

						sort(gapGenomeTup.begin(), gapGenomeTup.end());
						StoreMinimizers<GenomeTuple, Tuple>( strands[refinedClusters[r].strand] + curReadEnd,
																								 subreadLength, gapK, 1, gapReadTup, false);
						sort(gapReadTup.begin(), gapReadTup.end());
						CompareLists(gapReadTup.begin(),
												 gapReadTup.end(),
												 gapGenomeTup.begin(),
												 gapGenomeTup.end(),
												 gapPairs, gapOpts);
						//
						// Remove egregious off diagonal seeds
						//
						DiagonalSort<GenomeTuple>(gapPairs);
						CleanOffDiagonal(gapPairs, tinyOpts);


						for(int rm=0; rm < gapPairs.size(); rm++) {
							gapPairs[rm].first.pos  += curReadEnd;
							gapPairs[rm].second.pos += curGenomeEnd;
						}
						int gp=0;
						int nSaved =0;

						seqan::clear(gapSeedSet);
						int gpStart;
						int clusterIndex=0;


						fragments.clear();
						optFragmentIndices.clear();
						endpoints.clear();
						seqan::SeedSet<seqan::Seed<seqan::Simple>, seqan::Unordered> gapSeedSet;
						for (int m=0; m< gapPairs.size(); m++) {
							seqan::addSeed(gapSeedSet, 
														 seqan::Seed<seqan::Simple>(gapPairs[m].second.pos, 
																												gapPairs[m].first.pos, gapK),
														 seqan::Single());
						}
						seqan::String<seqan::Seed<seqan::Simple> > gapChain;
						if (seqan::length(gapSeedSet) > 0) {
							seqan::chainSeedsGlobally(gapChain, gapSeedSet, seqan::SparseChaining());
						}
						
						/*
						ofstream region;
						if (curReadEnd == 5187) {
							region.open("5187.tab");
						}
						for (gp=0; gp < gapPairs.size(); gp++) {

							fragments.push_back(Fragment(gapPairs[gp].second.pos,gapPairs[gp].first.pos,
																					 gapPairs[gp].second.pos+gapK,gapPairs[gp].first.pos+gapK, 1));
							if (curReadEnd == 5187) {
								region << gapPairs[gp].second.pos << "\t" << gapPairs[gp].first.pos << "\t" << gapPairs[gp].second.pos+gapK <<"\t" << gapPairs[gp].first.pos+gapK << "\t0" << endl;
							
							}
						}
						// Perform sparse chaining, uses time O(n log n).
			
						seqan::chainSeedsGlobally(chain, seedSet, seqan::SparseChaining());

						if (fragments.size() > 0) {
							GlobalChain(fragments, optFragmentIndices, endpoints);
						}
						StoreSubset(gapPairs, optFragmentIndices);


							if (curReadEnd == 5187) {
								for (int gpi=0; gpi < gapPairs.size(); gpi++) {
									region << gapPairs[gpi].second.pos << "\t" << gapPairs[gpi].first.pos << "\t" << gapPairs[gpi].second.pos+gapK <<"\t" << gapPairs[gpi].first.pos+gapK << "\t1" << endl;
								}
								region.close();
							}

						vector<Cluster> smallClusters;
						StoreDiagonalClustersLite(gapPairs, smallClusters, smallOpts, clusters[c].strand);
						sort(smallClusters.begin(), smallClusters.end());
						cout << "Removing paired indels with tiny opts" << endl;
						RemovePairedIndels(gapPairs, smallClusters, tinyOpts);
						
						//
						// Turn the indexed seeds into a gap chain.
						//
						seqan::String<seqan::Seed<seqan::Simple> > gapChain;

						for (int clusterIndex = 0; clusterIndex < optFragmentIndices.size(); clusterIndex++) {
							seqan::appendValue(gapChain, seqan::Seed<seqan::Simple>(fragments[optFragmentIndices[clusterIndex]].xl, 
																																			fragments[optFragmentIndices[clusterIndex]].yl, gapK));
						}
						*/														 
						refinedChains[c]=gapChain;						
					}
				}
			}
			
			//
			// Refine and store the alignment
			//
			//
			// The alignment is on a substring that starts at the beginning of the first chain.
			//
			GenomePos alnReadPos = seqan::beginPositionV(chain[0]);
			GenomePos alnRefPos  = seqan::beginPositionH(chain[0]);
			Alignment alignment;
			/*
			stringstream frc;
			frc << "full_ref_chaindiag_cluster." << r << ".tab";
			ofstream rdc(frc.str().c_str());
			*/
			for (int c = 0; chainLength> 0 and  c < chainLength-1; c++) {
				
				//
				// Chain is with respect to full sequence
				//
				GenomePos curGenomeEnd     = seqan::endPositionH(chain[c]);
				GenomePos nextGenomeStart  = seqan::beginPositionH(chain[c+1]);

				GenomePos curReadEnd       = seqan::endPositionV(chain[c]);
				GenomePos nextReadStart    = seqan::beginPositionV(chain[c+1]);
				int curRefinedReadEnd      = curReadEnd;
				int curRefinedGenomeEnd    = curGenomeEnd;
				int nextRefinedReadStart   = nextReadStart;
				int nextRefinedGenomeStart = nextGenomeStart;

				alignment.blocks.push_back(Block(seqan::beginPositionV(chain[c]),
																				 seqan::beginPositionH(chain[c]), glIndex.k));
				//				rdc << seqan::beginPositionH(chain[c]) << "\t" << seqan::beginPositionV(chain[c]) << "\t" << glIndex.k << "\t" << 0 << endl;
				for (int cs = 0; cs < seqan::length(refinedChains[c]); cs++) {
					//
					// Refined anchors are with respect to the chained sequence
					nextRefinedReadStart   = seqan::beginPositionV(refinedChains[c][cs]);
					nextRefinedGenomeStart = seqan::beginPositionH(refinedChains[c][cs]);


					//					rdc << seqan::beginPositionH(refinedChains[c][cs]) << "\t" << seqan::beginPositionV(refinedChains[c][cs]) << "\t" << glIndex.k <<  "\t" << c << endl;

					int m, rg, gg;
					SetMatchAndGaps(curRefinedReadEnd, nextRefinedReadStart,
													curRefinedGenomeEnd, nextRefinedGenomeStart, m, rg, gg);

					if (m > 0) {
						Alignment aln;
						RefineSubstrings(strands[refinedClusters[r].strand], curRefinedReadEnd, nextRefinedReadStart,
														 genome.seqs[chromIndex], curRefinedGenomeEnd, nextRefinedGenomeStart,
														 scoreMat, pathMat, aln);
						alignment.blocks.insert(alignment.blocks.end(), aln.blocks.begin(), aln.blocks.end());
						aln.blocks.clear();
					}



					curRefinedReadEnd   = seqan::endPositionV(refinedChains[c][cs]);
					curRefinedGenomeEnd = seqan::endPositionH(refinedChains[c][cs]);
					alignment.blocks.push_back(Block(nextRefinedReadStart, nextRefinedGenomeStart, 
																					 curRefinedReadEnd - nextRefinedReadStart));
				}

				// Add the last gap, or the only one if no refinements happened here.
							 
				int match, readGap, genomeGap;
				SetMatchAndGaps(curRefinedReadEnd, nextReadStart,
												curRefinedGenomeEnd, nextGenomeStart, match, readGap, genomeGap);
				
				if (match > 0) {
					Alignment aln;
					RefineSubstrings(strands[refinedClusters[r].strand], curRefinedReadEnd, nextReadStart,
													 genome.seqs[chromIndex], curRefinedGenomeEnd, nextGenomeStart,
													 scoreMat, pathMat, aln);
					alignment.blocks.insert(alignment.blocks.end(), aln.blocks.begin(), aln.blocks.end());
					aln.blocks.clear();					
				}
			}
			//			rdc.close();
			
			alignment.blocks.push_back(Block(seqan::beginPositionV(chain[chainLength-1]),
																			 seqan::beginPositionH(chain[chainLength-1]), 
																			 glIndex.k));

			string qStr, matStr, gStr;
			if (opts.viewPairwise) {
				alignment.CreateAlignmentStrings(strands[refinedClusters[r].strand],
																				 genome.seqs[chromIndex],
																				 qStr, matStr, gStr);
				
				
				int i=0;
				int q=0;
				int t=0;
				int nBlocks = alignment.blocks.size();
				
				while (i < qStr.size()) {
					int end = min((int) qStr.size(), i+50);
					string qsub = qStr.substr(i,end-i);
					std::cout.width(10);
					std::cout << q + alignment.GetQStart() << " q: " << qsub << endl;
					q+= qsub.size() - count(qsub.begin(),qsub.end(),'-');
					std::cout << "              " << matStr.substr(i,end-i) << endl;
					string tsub = gStr.substr(i,end-i);
					std::cout.width(10);
					std::cout << t + alignment.GetTStart() << " t: " << tsub << endl;
					t+= tsub.size() - count(tsub.begin(), tsub.end(),'-');
					cout <<endl;
					i=end;
				}
			}

			/*
			int pv=0, ph=0;
			int start = seqan::beginPositionV(refinedChains[0][0]);

			for (int c=0; c < refinedChains.size(); c++ ) {
				for (int cs = 0; cs < seqan::length(refinedChains[c]); cs++) {
					cout << c << "\t" << cs << "\t" 
							 << seqan::beginPositionV(refinedChains[c][cs]) << "\t" 
							 << seqan::beginPositionH(refinedChains[c][cs]) << "\t"
							 << seqan::beginPositionV(refinedChains[c][cs]) - pv << "\t"
							 << seqan::beginPositionH(refinedChains[c][cs]) - ph << "\t"
							 << seqan::beginPositionV(refinedChains[c][cs]) - start << endl;
					pv = seqan::beginPositionV(refinedChains[c][cs]);
					ph = seqan::beginPositionH(refinedChains[c][cs]);
				}
			}
			*/

			int nm=0;
			for(int b=0; b< alignment.blocks.size(); b++) {
				nm+= alignment.blocks[b].length;
			}
			cout << header.names[chromIndex] << "\t" 
					 << refinedClusters[r].tStart << "\t" 
					 << refinedClusters[r].tEnd << "\t"
					 << ks->name.s << "\t" << ks->seq.l << "\t" << r << "\t" 
				   << clusters[r].matches.size() << "\t"
					 << refinedClusters[r].matches.size() << "\t" << nm << endl;

				}

		//
		// Done with one read. Clean memory.
		//
		delete[] readRC;
	}
}

void HelpStore() {
	cout << "Usage: lsa index file.fa [options]" << endl
			 << "   -w (int) Minimizer window size (10)." << endl
			 << "   -f (int) Maximum minimizer frequency (200)." << endl
			 << "   -k (int) Word size" << endl
			 << "   -h Print help." << endl;
	
}
void HelpStoreLocal() {
	cout << "Usage: lsa local file.fa [options]" << endl
			 << "   -w (int) Minimizer window size (10)." << endl
			 << "   -f (int) Maximum minimizer frequency (5)." << endl
			 << "   -k (int) Word size (10)" << endl
			 << "   -h Print help." << endl;
}

void RunStoreLocal(int argc, char* argv[], LocalIndex &glIndex, Options &opts) {
	int argi = 0;
	string genome;
	string indexFile="";
	bool printIndex = false;
	opts.w=10;
	for (argi = 0; argi < argc; ) {
		if (ArgIs(argv[argi], "-h")) {
			HelpStoreLocal();
			exit(1);
		}
		else if (ArgIs(argv[argi], "-k")) {
			glIndex.k=atoi(argv[++argi]);
		}
		else if (ArgIs(argv[argi], "-w")) {
			glIndex.w=atoi(argv[++argi]);
		}
		else if (ArgIs(argv[argi], "-f")) {
			glIndex.maxFreq=atoi(argv[++argi]);
		}
		else if (strlen(argv[argi]) > 0 and argv[argi][0] == '-') {
			HelpStoreLocal();
			cout << "Invalid option " << argv[argi] << endl;
			exit(1);
		}
		else {
			genome = argv[argi];
			cerr << "genome " << genome << endl;
		}
		++argi;
	}
	if (genome == "") {
		HelpStore();
		exit(1);
	}


	glIndex.IndexFile(genome);
	glIndex.Write(genome + ".gli");
}

void RunStore(int argc, char* argv[], vector<GenomeTuple> &minimizers, Header &header, Options &opts) {
	// open query file for reading; you may use your favorite FASTA/Q parser
	int argi = 0;
	string genome;
	string indexFile="";
	bool printIndex = false;
	bool compress=false;
	opts.w=10;
	for (argi = 0; argi < argc; ) {
		if (ArgIs(argv[argi], "-w")) {
			++argi;
			opts.w = atoi(argv[argi]);
		}
		else if (ArgIs(argv[argi], "-f")) {
			++argi;
			opts.maxFreq = atoi(argv[argi]);
		}		
		else if (ArgIs(argv[argi], "-i")) {
			++argi;
			indexFile=argv[argi];
		}
		else if (ArgIs(argv[argi], "-c")) {
			cerr << "WARNING: Compressing index" << endl;
			compress = true;
		}
		else if (ArgIs(argv[argi], "-p")) {
			++argi;
			printIndex = true;
		}
		else if (ArgIs(argv[argi], "-k")) {
			++argi;
			opts.k=atoi(argv[argi]);
		}		
		else if (ArgIs(argv[argi], "-h")) {
			++argi;
			HelpStore();
			exit(0);
		}		

		else if (strlen(argv[argi]) > 0 && argv[argi][0] == '-') {
			HelpStore();
			cout << "Invalid option " << argv[argi] << endl;
			exit(1);
		}
		else {
			genome = argv[argi];
			cerr << "genome " << genome << endl;
		}
		++argi;
	}
	if (genome == "") {
		HelpStore();
		exit(1);
	}
	if (indexFile == "") {
		indexFile = genome + ".mmi";
	}

	if (printIndex and ReadIndex(indexFile, minimizers, header, opts)) {
		PrintIndex(minimizers, opts.k);
		exit(0);
	}

	StoreIndex(genome, minimizers, header, opts);
	WriteIndex(indexFile, minimizers, header, opts);
}


void Usage() {
	cout << "Program: lsa (long sequence alignment)" << endl;
	cout << "Version: beta" << endl;
	cout << "Contact: Mark Chaisson (mchaisso@usc.edu)" << endl << endl;
	cout << "Usage:   lsa <command> [options]"<< endl << endl;
	cout << "Command: index   - Build a mm index on sequences." << endl;
	cout << "         align   - Map reads using the index" << endl;
	cout << "         local   - Build local index" << endl;
}

int main(int argc, char *argv[]) {
	if (argc < 2) {
		Usage();
		return 1;
	}

	Options opts;
	opts.k=21;

  int argi;
	vector<GenomeTuple>  minimizers;
	LocalIndex lIndex;
	Header header;
	for (argi = 1; argi < argc; ){
		if (ArgIs(argv[argi], "index")) {
			argc -=2;
      RunStore(argc,  &argv[2], minimizers, header, opts);		
			exit(0);
		}
		else if (ArgIs(argv[argi], "align")) {
			argc -=2;
			RunAlign(argc, &argv[2], opts);
			exit(0);
		}
		else if (ArgIs(argv[argi], "local")) {
			argc -=2;
			RunStoreLocal(argc, &argv[2], lIndex, opts);
			exit(0);
		}

		else {
			Usage();
			exit(1);
		}
	}

}
