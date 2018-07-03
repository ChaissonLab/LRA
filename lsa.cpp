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
#include "SeqUtils.h"
#include "TupleOps.h"
#include "MinCount.h"
#include "MMIndex.h"
#include "Options.h"
#include "Alignment.h"
#include "MapRead.h"
#include "Input.h"
#include <math.h>

const char* version="0.1-alpha";

bool ArgIs(const char* a, const char* b) {
	return strcmp(a,b) == 0;
}

void HelpMap() {
	cout << "Usage: lsc align genome.fa reads.fa [options]" << endl;
	cout << "Options:" << endl
			 << "  -p  (flag)  View pairwise alignment." << endl;
}
		


void RunAlign(int argc, char* argv[], Options &opts ) {
// open query file for reading; you may use your favorite FASTA/Q parser
	int argi = 0;
	string genomeFile = "", reads = "";
	string indexFile="";
	int w=10;

	for (argi = 0; argi < argc; ) {
		if (ArgIs(argv[argi], "-i")) {
			++argi;
			indexFile=argv[argi];
		}		
		else if (ArgIs(argv[argi], "-a")) {
			++argi;
			opts.storeAll=true;
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
		else if (ArgIs(argv[argi], "-H")) {
			opts.hardClip=true;
		}
		else if (ArgIs(argv[argi], "-p")) {
			opts.printFormat = argv[++argi][0];
		}
		else if (ArgIs(argv[argi], "-n")) {
			opts.bestn = atoi(argv[++argi]);
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

	Input reader;
	reader.Initialize(reads);
	int offset=0;
	Read read;

	if (opts.printFormat == 's') {
		stringstream cl;
		cl << "lsa align";
		for (int i=0; i < argc; i++) {
			cl << " " << argv[i];
		}
		cout << "@PG\tID:lsa\tPN:lsa\tVN:"<<version<<"\tCL:"<<cl.str() << endl;
		genome.header.WriteSAMHeader(cout);
	}
	while (reader.GetNext(read)) {
		MapRead(read, genome, genomemm, glIndex, opts, cout);
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
