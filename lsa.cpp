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
#include "htslib/sam.h"
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
const char* GetArgv(const char* argv[], int argc, int argi) {
	if (argi +1 >= argc) {
		cout << "ERROR, argument " << argv[argi] << " requires a value." << endl;
		exit(1);
	}
	return argv[argi+1];
}


void HelpMap() {
	cout << "Usage: lsc align genome.fa reads.fa [options]" << endl;
	cout << "Options:" << endl
			 << "   -p  [FMT]   Print alignment format FMT='b' bed, 's' sam 'p' pair ." << endl
			 << "   -H          Use hard-clipping for SAM output format" << endl
			 << "   -M  M(int)  Do not refine clusters with fewer than M global matches (20)." << endl
			 << "   -m  m(int)  Do not align clusters with fewer than m refined"<< endl
			 << "               matches (40). Typically m > 3*M" << endl
			 << "   -a          Query all positions in a read, not just minimizers. " << endl
			 << "               This is ~20% slower, with an increase in specificity " << endl
			 << "               to be determined." << endl;
}
		
class MapInfo {
public:
	Header*header;
	Genome *genome;
	vector<GenomeTuple> *genomemm;
	LocalIndex *glIndex;
	Input *reader;
	Options *opts;
	ostream *out;
	sem_t *semaphore;
};

void MapReads(MapInfo *mapInfo) {
	Read read;
	
	while (mapInfo->reader->GetNext(read)) {
		MapRead(read, 
						*mapInfo->genome, 
						*mapInfo->genomemm, 
						*mapInfo->glIndex, 
						*mapInfo->opts, 
						*mapInfo->out,
						mapInfo->semaphore);
	}
	pthread_exit(NULL); 
}



void RunAlign(int argc, const char* argv[], Options &opts ) {
// open query file for reading; you may use your favorite FASTA/Q parser
	int argi = 0;
	string genomeFile = "", reads = "";
	string indexFile="";
	int w=10;


	for (argi = 0; argi < argc; ) {
		if (ArgIs(argv[argi], "-a")) {
			++argi;
			opts.storeAll=true;
		}		
		else if (ArgIs(argv[argi], "-W")) {
			opts.globalW=atoi(GetArgv(argv, argc, argi));
			++argi;
		}		
		else if (ArgIs(argv[argi], "-M")) {
			opts.minClusterSize=atoi(GetArgv(argv, argc, argi));
			++argi;
		}		
		else if (ArgIs(argv[argi], "-m")) {
			opts.minRefinedClusterSize=atoi(GetArgv(argv, argc, argi));
			++argi;
		}		
		else if (ArgIs(argv[argi], "-f")) {
			opts.globalMaxFreq=atoi(GetArgv(argv, argc, argi));
			++argi;
		}		
		else if (ArgIs(argv[argi], "-K")) {
			opts.globalK=atoi(GetArgv(argv, argc, argi));
			++argi;
		}		
		else if (ArgIs(argv[argi], "-H")) {
			opts.hardClip=true;
		}
		else if (ArgIs(argv[argi], "-p")) {			
			opts.printFormat = GetArgv(argv, argc, argi)[0];
			++argi;
		}
		else if (ArgIs(argv[argi], "-n")) {
			opts.bestn=atoi(GetArgv(argv, argc, argi));
		}
		else if (ArgIs(argv[argi], "-t")) {
			opts.nproc=atoi(GetArgv(argv, argc, argi));
			++argi;
		}		

		else if (ArgIs(argv[argi], "-o")) {
			opts.outfile = argv[++argi];
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
	ostream *outPtr;
	ofstream outfile;
	if (opts.outfile == "") {
		outPtr = &cout;
	}
	else {
		outfile.open(opts.outfile.c_str());
		outPtr = &outfile;
	}

	if (opts.printFormat == 's') {
		stringstream cl;
		cl << "lsa align";
		for (int i=0; i < argc; i++) {
			cl << " " << argv[i];
		}
		cout << "@PG\tID:lsa\tPN:lsa\tVN:"<<version<<"\tCL:"<<cl.str() << endl;
		genome.header.WriteSAMHeader(*outPtr);
	}

	pthread_attr_t *threadAttr = new pthread_attr_t[opts.nproc];
	for (int procIndex = 0; procIndex < opts.nproc; procIndex++ ){
		pthread_attr_init(&threadAttr[procIndex]);
	}
	if (opts.nproc > 1) {
		pthread_t *threads = new pthread_t[opts.nproc];
		MapInfo mapInfo;
		mapInfo.genome = &genome;
		mapInfo.genomemm = &genomemm;
		mapInfo.glIndex = &glIndex;
		mapInfo.reader = &reader;
		mapInfo.opts= &opts;
		mapInfo.out = outPtr;

		mapInfo.semaphore = sem_open("/writer",     O_CREAT, 0644, 1);


		for (int procIndex = 0; procIndex < opts.nproc; procIndex++ ){ 
			pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*))MapReads, &mapInfo);
		}

		for (int procIndex = 0; procIndex < opts.nproc; procIndex++) {
			pthread_join(threads[procIndex], NULL);
		}
		
	}
	else {
		while (reader.GetNext(read)) {
			MapRead(read, genome, genomemm, glIndex, opts, *outPtr);
		}
	}
}
void HelpStoreIndex() {
	cout << "Usage: lsa index file.fa [options]" << endl
			 << "  Global index options " << endl
			 << "   -W (int) Minimizer window size (10)." << endl
			 << "   -F (int) Maximum minimizer frequency (200)." << endl
			 << "   -K (int) Word size" << endl
			 << "  Local index options: "<< endl
			 << "   -w (int) Local minimizer window size (10)." << endl
			 << "   -f (int) Local maximum minimizer frequency (5)." << endl
			 << "   -k (int) Local word size (10)" << endl
			 << "   -h Print help." << endl;
}

void HelpStoreGlobal() {
	cout << "Usage: lsa index file.fa [options]" << endl
			 << "   -W (int) Minimizer window size (10)." << endl
			 << "   -F (int) Maximum minimizer frequency (200)." << endl
			 << "   -K (int) Word size" << endl
			 << "   -h Print help." << endl;
	
}
void HelpStoreLocal() {
	cout << "Usage: lsa local file.fa [options]" << endl
			 << "   -w (int) Local minimizer window size (10)." << endl
			 << "   -f (int) Local maximum minimizer frequency (5)." << endl
			 << "   -k (int) Local word size (10)" << endl
			 << "   -h Print help." << endl;
}

void RunStoreLocal(int argc, const char* argv[], 
									 LocalIndex &glIndex, Options &opts) {
	int argi = 0;
	string genome;
	string indexFile="";
	bool printIndex = false;
	opts.localK=glIndex.k;
	for (argi = 0; argi < argc; ) {
		if (ArgIs(argv[argi], "-h")) {
			HelpStoreLocal();
			exit(1);
		}
		else if (ArgIs(argv[argi], "-k")) {
			opts.localK=atoi(argv[++argi]);
			glIndex.k=opts.localK;
		}
		else if (ArgIs(argv[argi], "-w")) {
			opts.localW=atoi(argv[++argi]);
			glIndex.w=opts.localW;
		}
		else if (ArgIs(argv[argi], "-f")) {
			opts.localMaxFreq=atoi(argv[++argi]);
			glIndex.maxFreq=opts.localMaxFreq;
		}
		else if (ArgIs(argv[argi], "-K") or ArgIs(argv[argi], "-W") or ArgIs(argv[argi], "-F")) {
			argi+=2;
			continue;
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
		HelpStoreGlobal();
		exit(1);
	}

	glIndex.IndexFile(genome);
	glIndex.Write(genome + ".gli");
}

void RunStoreGlobal(int argc, const char* argv[], 
										vector<GenomeTuple> &minimizers, Header &header, Options &opts) {
	// open query file for reading; you may use your favorite FASTA/Q parser
	int argi = 0;
	string genome;
	string indexFile="";
	bool printIndex = false;
	bool compress=false;
	opts.globalW=10;
	for (argi = 0; argi < argc; ) {
		if (ArgIs(argv[argi], "-W")) {
			++argi;
			opts.globalW = atoi(argv[argi]);
		}
		else if (ArgIs(argv[argi], "-F")) {
			++argi;
			opts.globalMaxFreq = atoi(argv[argi]);
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
		else if (ArgIs(argv[argi], "-K")) {
			++argi;
			opts.globalK=atoi(argv[argi]);
		}		
		else if (ArgIs(argv[argi], "-k") or ArgIs(argv[argi], "-w") or ArgIs(argv[argi], "-f")) {
			argi+=2;
			continue;
		}
		else if (ArgIs(argv[argi], "-h")) {
			++argi;
			HelpStoreGlobal();
			exit(0);
		}		

		else if (strlen(argv[argi]) > 0 && argv[argi][0] == '-') {
			HelpStoreGlobal();
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
		HelpStoreGlobal();
		exit(1);
	}
	if (indexFile == "") {
		indexFile = genome + ".mmi";
	}

	if (printIndex and ReadIndex(indexFile, minimizers, header, opts)) {
		PrintIndex(minimizers, opts.globalK);
		exit(0);
	}

	StoreIndex(genome, minimizers, header, opts);
	WriteIndex(indexFile, minimizers, header, opts);
}

void RunStoreIndex(int argc, const char* argv[]) {
	LocalIndex glIndex;
	vector<GenomeTuple> minimizers;
	Header header;
	Options opts;

	RunStoreGlobal(argc, argv, minimizers, header, opts);
  RunStoreLocal(argc, argv, glIndex, opts);
}



void Usage() {
	cout << "Program: lsa (long sequence alignment)" << endl;
	cout << "Version: beta" << endl;
	cout << "Contact: Mark Chaisson (mchaisso@usc.edu)" << endl << endl;
	cout << "Usage:   lsa <command> [options]"<< endl << endl;
	cout << "Command: index   - Build global and local indexes on a genome." << endl;
	cout << "         align   - Map reads using the index" << endl;
	cout << "         global  - Build a global index." << endl;
	cout << "         local   - Build local index" << endl;
}

int main(int argc, const char *argv[]) {
	if (argc < 2) {
		Usage();
		return 1;
	}

	Options opts;

  int argi;
	vector<GenomeTuple>  minimizers;
	LocalIndex lIndex;
	Header header;
	for (argi = 1; argi < argc; ){
		if (ArgIs(argv[argi], "index")) {
			argc -=2;
      RunStoreIndex(argc,  &argv[2]);
			exit(0);
		}
		else if (ArgIs(argv[argi], "global")) {
			argc -=2;
      RunStoreGlobal(argc,  &argv[2], minimizers, header, opts);		
			exit(0);
		}
		else if (ArgIs(argv[argi], "local")) {
			argc -=2;
			RunStoreLocal(argc, &argv[2], lIndex, opts);
			exit(0);
		}
		else if (ArgIs(argv[argi], "align")) {
			argc -=2;
			RunAlign(argc, (const char**) &argv[2], opts);
			exit(0);
		}

		else {
			Usage();
			exit(1);
		}
	}

}
