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

bool Is(const char* a, const char* b) {
	return strcmp(a,b) == 0;
}


void HelpMap() {
	cout << "Usage: lsc map genome.fa reads.fa" << endl;
}
void RunMap(int argc, char* argv[], Options &opts ) {
	// open query file for reading; you may use your favorite FASTA/Q parser
	int argi = 0;
	string genome = "", reads = "";
	opts.maxFreq = 200;
	string indexFile="";
	int w=10;
	bool storeAll = false;
	for (argi = 0; argi < argc; ) {
		if (Is(argv[argi], "-i")) {
			++argi;
			indexFile=argv[argi];
		}		
		if (Is(argv[argi], "-a")) {
			++argi;
			storeAll=true;
		}		
		if (Is(argv[argi], "-w")) {
			++argi;
			opts.w=atoi(argv[argi]);
		}		
		if (Is(argv[argi], "-k")) {
			++argi;
			opts.k=atoi(argv[argi]);
		}		
		
		else {
			if (genome == "") {
				genome = argv[argi];
			}
			else if (reads == "") {
				reads = argv[argi];
			}
		}
		++argi;
	}

	if (genome == "" || reads == "") {
		HelpMap();
		exit(1);
	}
	if (indexFile == "") {
		indexFile = genome + ".mmi";
	}

	vector<GenomeTuple> genomemm;
	if (ReadIndex(indexFile, genomemm, opts) == 0) {
		StoreIndex(genome, genomemm, opts);
	}
	
	gzFile f = gzopen(reads.c_str(), "r");
	kseq_t *ks = kseq_init(f);
	int offset=0;
	vector<GenomeTuple> readmm;
	vector<pair<GenomeTuple, GenomeTuple> > matches;
	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		readmm.clear();
		matches.clear();
		if (storeAll) {
			StoreAll(ks, opts.k, readmm);
		}
		else {
			StoreMinimizers(ks, opts.k, opts.w, readmm);
		}
		CompareLists(readmm, genomemm, matches);
		DiagonalSort(matches);
		for (int m=0; m < matches.size(); m++) {
			cout << matches[m].first.pos << "\t" << matches[m].second.pos << "\t" << matches[m].first.tuple;
			if (m > 0) {
				cout << "\t" << matches[m].first.pos - matches[m].second.pos - (matches[m-1].first.pos - matches[m-1].second.pos);
			}
			cout << endl;
		}
		cout << ks->name.s << "\t" << ks->seq.l << "\t" << matches.size() << endl;
	}
}



void RunStore(int argc, char* argv[], vector<GenomeTuple> &minimizers, Options &opts) {
	// open query file for reading; you may use your favorite FASTA/Q parser
	int argi = 0;
	string genome;
	int maxFreq = 200;
	string indexFile="";
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
		else if (Is(argv[argi], "-k")) {
			++argi;
			opts.k=atoi(argv[argi]);
			InitMask(opts.k);
		}		
		else {
			genome = argv[argi];
			cerr << genome << endl;
		}
		++argi;
	}
	
	if (indexFile == "") {
		indexFile = genome + ".mmi";
	}
	StoreIndex(genome, minimizers, opts);
	WriteIndex(indexFile, minimizers, opts);
}


void Usage() {
	cout << "Program: lsa (long sequence alignment)" << endl;
	cout << "Version: beta" << endl;
	cout << "Contact: Mark Chaisson (mchaisso@usc.edu)" << endl << endl;
	cout << "Usage:   lsa <command> [options]"<< endl << endl;
	cout << "Command: index   - Build a mm index on sequences." << endl;
	cout << "         map     - Map reads using the index" << endl;
}
int main(int argc, char *argv[]) {

	if (argc < 2) {
		Usage();
		return 1;
	}
	InitSeqMap();
	Options opts;
	opts.k=21;
	InitMask(opts.k);

  int argi;
	vector<GenomeTuple>  minimizers;
	for (argi = 1; argi < argc; ){
		if (Is(argv[argi], "index")) {
			argc -=2;
      RunStore(argc,  &argv[2], minimizers, opts);		
			exit(0);
		}
		else if (Is(argv[argi], "map")) {
			argc -=2;
			RunMap(argc, &argv[2], opts);
			exit(0);
		}
		else {
			Usage();
			exit(1);
		}
	}

}
