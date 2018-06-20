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

bool Is(const char* a, const char* b) {
	return strcmp(a,b) == 0;
}

void HelpMap() {
	cout << "Usage: lsc map genome.fa reads.fa" << endl;
}

void SwapStrand(kseq_t* read, Options &opts, vector<pair<GenomeTuple, GenomeTuple> > &matches, int start, int end) {
	for (int m=start; m < end; m++) {
		matches[m].first.pos = read->seq.l - (matches[m].first.pos + opts.k - 1);
	}
}



int SetStrand(kseq_t *read, Genome &genome, Options &opts, vector<pair<GenomeTuple, GenomeTuple> > &matches, int start, int end) {
	int nSame=0;
	int nDifferent=0;
	for (int m=start; m< end; m++) {
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

void SetClusterStrand(kseq_t* read, Genome &genome, Options &opts, 
											vector<pair<GenomeTuple, GenomeTuple> > &matches, 			
											vector<Cluster> &clusters) {
	for (int c = 0; c < clusters.size(); c++) {
		clusters[c].strand = SetStrand(read, genome, opts, matches, clusters[c].start, clusters[c].end);
		if (clusters[c].strand == 1) {
			SwapStrand(read, opts, matches, clusters[c].start, clusters[c].end);
		}
	}
}

void RunMap(int argc, char* argv[], Options &opts ) {
	// open query file for reading; you may use your favorite FASTA/Q parser
	int argi = 0;
	string genomeFile = "", reads = "";
	string indexFile="";
	int w=10;
	bool storeAll = false;
	for (argi = 0; argi < argc; ) {
		if (Is(argv[argi], "-i")) {
			++argi;
			indexFile=argv[argi];
		}		
		else if (Is(argv[argi], "-a")) {
			++argi;
			storeAll=true;
		}		
		else if (Is(argv[argi], "-w")) {
			++argi;
			opts.w=atoi(argv[argi]);
		}		
		else if (Is(argv[argi], "-m")) {
			++argi;
			opts.minClusterSize= atoi(argv[argi]);
		}		
		else if (Is(argv[argi], "-f")) {
			++argi;
			opts.maxFreq = atoi(argv[argi]);
		}		
		else if (Is(argv[argi], "-k")) {
			++argi;
			opts.k=atoi(argv[argi]);
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

		DiagonalSort<GenomeTuple>(matches);

		CleanOffDiagonal(matches, opts);
		vector<Cluster> clusters;
		PrintDiagonal(matches);
		StoreDiagonalClusters(matches, clusters, opts);
		
		//
		// Add pointers to seq that make code more readable.
		//
		char *readRC;
		char *read=ks->seq.s;
		int readLen = ks->seq.l;
		CreateRC(read, readLen, readRC);

		SetClusterStrand(ks, genome, opts,
										 matches, clusters);

		LocalIndex forwardIndex;
		LocalIndex reverseIndex;
		forwardIndex.IndexSeq(read, readLen);
		reverseIndex.IndexSeq(readRC, readLen);

		for (int c = 0; c < clusters.size(); c++) {

			if (clusters[c].start == clusters[c].end) {
				continue;
			}
			
			int cs = clusters[c].start;
			int ce = clusters[c].end;

			int chromIndex = genome.header.Find(matches[cs].second.pos);
			GenomePos chromOffset = genome.header.pos[chromIndex];
			cout << c << "\t" << ce - cs <<"\t" << genome.header.names[chromIndex] 
					 << "\t" << matches[cs].second.pos - chromOffset 
					 << "\t" << matches[ce].second.pos - chromOffset 
					 << "\t" << chromIndex 
					 << "\t" << matches[cs].second.pos << endl;
			vector<pair<GenomeTuple, GenomeTuple> > tmpm;
			vector<pair<GenomeTuple, GenomeTuple> > localMatches;
			for (int m = cs; m < ce; m++) { tmpm.push_back(matches[m]);}
				
			PrintPairs(tmpm, opts.k, c);

			//
			// Get the boundaries of the cluster in both sequences.
			//
			CartesianTargetSort<GenomeTuple>(matches.begin()+cs, matches.begin()+ce);
			
			//
			// Get shorthand access to alignment boundaries.
			//

		  GenomePos qs, qe, ts, te;
			qs = matches[cs].first.pos;
			ts = matches[cs].second.pos;

			qe = matches[ce-1].first.pos + opts.k;
			te = matches[ce-1].second.pos + opts.k;

			int ls, le;
			GenomePos chromStartOffset = header.GetOffset(ts);
			GenomePos chromEndOffset = header.GetNextOffset(te);
			// Search region starts in window, or beginning of chromosome
			GenomePos wts, wte;
			if ( chromStartOffset + opts.window > ts ) {wts = chromStartOffset;	}
			else {wts = ts - opts.window; }
				
			if (te + opts.window > chromEndOffset) {wte = chromEndOffset-1;	}
			else { wte = te + opts.window;	}
			
			ls = glIndex.LookupIndex(wts);
			le = glIndex.LookupIndex(wte);
				

			LocalIndex *readIndex;

			if (clusters[c].strand == 0) {
				readIndex = &forwardIndex;
			}
			else {
				readIndex = &reverseIndex;
			}
			int lmIndex=0;
			for (int lsi=ls; lsi <= le; lsi++) {
				//
				// Find the coordinates in the cluster that start in this local index.
				//
				GenomePos genomeIndexStart = glIndex.offsets[lsi];
				GenomePos genomeIndexEnd   = glIndex.offsets[lsi+1]-1;
				int matchStart = CartesianTargetLowerBound<GenomeTuple>(matches.begin()+cs, 
																																matches.begin()+ce,
																																genomeIndexStart);

				int matchEnd   = CartesianTargetUpperBound<GenomeTuple>(matches.begin()+cs, 
																																matches.begin()+ce, 
																																genomeIndexEnd);

				GenomePos readStart = matches[matchStart+cs].first.pos;
				if (lsi == ls) {
					if (readStart < opts.window) {
						readStart = 0;
					}
					else {
						readStart -= opts.window;
					}
				}
				GenomePos readEnd;
				if (matchEnd > matchStart) {readEnd = matches[matchEnd-1+cs].first.pos;	}
				else { readEnd = matches[matchStart+cs].first.pos + opts.k;}
				//
				// Expand boundaries of read to match.
				if (lsi == le) {
					if (readEnd + opts.window > ks->seq.l) { readEnd = ks->seq.l; }
					else { readEnd += opts.window;	}
				}			
				
				//
				// Find the boundaries of where in the query the matches should be added.
				//
				GenomePos minQuery=-1;
				GenomePos maxQuery=0;
				
				for (int m=cs; m < ce; m++) {
					if (minQuery == -1) {
						if ( lsi == ls ){
							minQuery = readStart;
						}
						else {
							minQuery = matches[cs].first.pos;
						}
					}
					else {
						minQuery = min(minQuery, matches[m].first.pos);
					}

					if (maxQuery == 0) {
						if (lsi == le) {
							maxQuery = readEnd-1;
						}
						else {
							maxQuery=max(maxQuery, matches[m].first.pos);
						}
					}
				}
				
				int queryIndexStart = readIndex->LookupIndex(minQuery);
				int queryIndexEnd   = readIndex->LookupIndex(maxQuery);
				Options smallOpts = opts;
				smallOpts.maxFreq=3;
				smallOpts.maxDiag=25;

				for (int qi = queryIndexStart; qi <= queryIndexEnd; qi++) {
					vector<pair<LocalTuple, LocalTuple> > smallMatches;

					GenomePos qStartBoundary=readIndex->boundaries[qi];
					GenomePos qEndBoundary=readIndex->boundaries[qi+1];
					

					CompareLists<LocalTuple>(readIndex->minimizers.begin() + qStartBoundary,
																	 readIndex->minimizers.begin() + qEndBoundary,
																	 glIndex.minimizers.begin()+glIndex.boundaries[lsi], 
																	 glIndex.minimizers.begin()+glIndex.boundaries[lsi+1], smallMatches, smallOpts);

					DiagonalSort<LocalTuple>(smallMatches);
					CleanOffDiagonal(smallMatches, smallOpts);

					localMatches.resize(localMatches.size()+smallMatches.size());
#ifdef _TESTING_
					CartesianSort<LocalTuple>(smallMatches);
#endif
					for(int i=0; i < smallMatches.size(); i++) {
						localMatches[lmIndex+i].first.pos  = smallMatches[i].first.pos+readStart;
						localMatches[lmIndex+i].second.pos = smallMatches[i].second.pos+genomeIndexStart;
						localMatches[lmIndex+i].first.t = localMatches[lmIndex+i].second.t = smallMatches[i].second.t;
					}
					lmIndex+=smallMatches.size();
				}

#ifdef _TESTING_				
				PrintPairs(localMatches, glIndex.k);
#endif
			}
			cout << ks->name.s << "\t" << ks->seq.l << "\t" << c << "\t" << clusters[c].size() << "\t" << localMatches.size() << endl;
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
			 << "   -i (string) Index file for alternative store." << endl
			 << "   -k (int) Word size" << endl
			 << "   -h Print help." << endl;
	
}
void HelpStoreLocal() {
	cout << "Usage lsa local file.fa";
}

void RunStoreLocal(int argc, char* argv[], LocalIndex &glIndex, Options &opts) {
	int argi = 0;
	string genome;
	string indexFile="";
	bool printIndex = false;
	opts.w=10;
	for (argi = 0; argi < argc; ) {
		if (strlen(argv[argi]) > 0 and argv[argi][0] == '-') {
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
	opts.w=10;
	for (argi = 0; argi < argc; ) {
		if (Is(argv[argi], "-w")) {
			++argi;
			opts.w = atoi(argv[argi]);
		}
		else if (Is(argv[argi], "-f")) {
			++argi;
			opts.maxFreq = atoi(argv[argi]);
		}		
		else if (Is(argv[argi], "-i")) {
			++argi;
			indexFile=argv[argi];
		}
		else if (Is(argv[argi], "-p")) {
			++argi;
			printIndex = true;
		}
		else if (Is(argv[argi], "-k")) {
			++argi;
			opts.k=atoi(argv[argi]);
		}		
		else if (Is(argv[argi], "-h")) {
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
	cout << "         map     - Map reads using the index" << endl;
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
		if (Is(argv[argi], "index")) {
			argc -=2;
      RunStore(argc,  &argv[2], minimizers, header, opts);		
			exit(0);
		}
		else if (Is(argv[argi], "map")) {
			argc -=2;
			RunMap(argc, &argv[2], opts);
			exit(0);
		}
		else if (Is(argv[argi], "local")) {
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
