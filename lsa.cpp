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

void RunMap(int argc, char* argv[], Options &opts ) {
	// open query file for reading; you may use your favorite FASTA/Q parser
	int argi = 0;
	string genome = "", reads = "";
	string indexFile="";
	int w=10;
	bool storeAll = false;
	bool cartesianSort =false;
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
		else if (Is(argv[argi], "-c")) {
			++argi;
			cartesianSort= true;
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
	Header header;
	vector<GenomeTuple> genomemm;
	if (ReadIndex(indexFile, genomemm, header, opts) == 0) {
		StoreIndex(genome, genomemm, header, opts);
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
			StoreMinimizers<GenomeTuple, Tuple>(ks->seq.s, ks->seq.l, opts.k, opts.w, readmm);
		}
		sort(readmm.begin(), readmm.end());
		CompareLists(readmm, genomemm, matches, opts);
		if (cartesianSort) {
			CartesianTargetSort(matches);
		}
		else {
			DiagonalSort(matches);

			CleanOffDiagonal(matches, opts);
			vector< vector<pair<GenomeTuple, GenomeTuple> > > clusters;
			OffDiagonalClusters(matches, clusters, opts);

			for (int c =0; c < clusters.size(); c++) {
				if (clusters[c].size() == 0) {
					continue;
				}
				int l=clusters[c].size()-1;
				
				cout << ks->name.s << "\t" << ks->seq.l << "\t" << c << "\t" << clusters[c].size() << endl;
				/*
					
				for (int m=0; m < clusters[c].size(); m++) {
					string s;			
					
					clusters[c][m].first.ToString(opts.k,s);
					cout << clusters[c][m].first.pos << "\t" << clusters[c][m].second.pos << "\t" << clusters[c][m].first.tuple << "\t" << s;
					if (m > 0) {
						cout << "\t" << clusters[c][m].first.pos - clusters[c][m].second.pos - (clusters[c][m-1].first.pos - clusters[c][m-1].second.pos);
					}
					cout << endl;
				}
				cout << ks->name.s << "\t" << ks->seq.l << "\t" << clusters[c].size() << endl;
				*/
			}
		}
	}
}

void HelpStore() {
	cout << "Usage: lsa index file.fa [options]" << endl
			 << "   -w (int) Minimizer window size (10)." << endl
			 << "   -f (int) Maximum minimizer frequency (200)." << endl
			 << "   -i (string) Index file for alternative store." << endl
			 << "   -k (int) Word size" << endl;
}
void HelpStoreLocal() {
	cout << "Usage lsa local file.fa";
}

void RunStoreLocal(int argc, char* argv[], GenomeLocalIndex &glIndex, Options &opts) {
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


	StoreLocalIndex(genome, glIndex, opts);
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
	InitSeqMap();
	Options opts;
	opts.k=21;

  int argi;
	vector<GenomeTuple>  minimizers;
	GenomeLocalIndex glIndex;
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
			RunStoreLocal(argc, &argv[2], glIndex, opts);
			exit(0);
		}

		else {
			Usage();
			exit(1);
		}
	}

}
