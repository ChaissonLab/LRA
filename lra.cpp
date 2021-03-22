#include "htslib/hts.h"
#include "htslib/kseq.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <thread>
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
#include "Input.h"
#include "MMIndex.h"
#include "TupleOps.h"
#include "MinCount.h"
#include "MapRead.h"
#include <algorithm>
#include "SeqUtils.h"
#include "Options.h"
#include "Alignment.h"
#include "LogLookUpTable.h"


using namespace std;

int IO_BUFFER_SIZE=10000000;
const char* lraVersion="V1.1.2";

bool ArgIs(const char* a, const char* b) {
	return strcmp(a,b) == 0;
}
const char* GetArgv(const char* argv[], int argc, int argi) {
	if (argi + 1 >= argc) {
		cout << "ERROR, argument " << argv[argi] << " requires a value." << endl;
		exit(1);
	}
	return argv[argi + 1];
}


void HelpMap() {
	cout << "Usage: lra align [options] genome.fa reads [reads2 ...]" << endl << endl;
	cout << "   The genome should be indexed using the 'lra index' program." << endl
			 << "   'reads' may be either fasta, sam, or bam, and multiple input files may be given." << endl << endl;
	cout << "Options:" << endl
			 << "   -CCS (flag) Align CCS reads. " << endl
			 << "   -CLR (flag) Align CLR reads. " << endl
			 << "   -ONT (flag) Align Nanopore reads. " << endl
			 << "   -CONTIG (flag) Align large contigs." << endl
			 << "   -p  [FMT]   Print alignment format FMT='b' bed, 's' sam, 'p' PAF, 'pc' PAF with cigar, 'a' pairwise alignment." << endl
			 << "   -H          Use hard-clipping for SAM output format" << endl
     		 << "   -Flag  F(int)  Skip reads with any flags in F set (bam input only)." << endl
     		 << "   -t  n(int)   Use n threads (1)" << endl
			 //<< "   -M  M(int)  Do not refine clusters with fewer than M global matches (20)." << endl
			//<< "   -m  m(int)  Do not align clusters with fewer than m refined"<< endl
			// << "               matches (40). Typically m > 3*M" << endl
			 << "   -a  (flag)  Query all positions in a read, not just minimizers. " << endl
			 //<< "               This is 10-20% slower, with an increase in specificity. " << endl
			 // << "   -b  (flag)  Skip banded alignment. This is about a 15% speedup." << endl
			 << "   -SV  (int) (path to svsig file)  Print sv signatures for each alignment with length above the given threshold (DEFAULT:25). And the path of output svsig file" << endl
			 //<< "   -R  (flag)  MeRge clusters before sparse dynamic programming." << endl
			 //<< "   -N  (flag)  Use Naive dynamic programming to find the global chain." << endl
			// << "	-S 	(flag)  Use Sparse dynamic programming to find the global chain." << endl
			// << "	-T 	(flag)  Use log LookUpTable when gap length is larger than 501." << endl
			 << "   -at  (float) a float in (0, 1), Threshold to decide secondary alignments based on chaining value (DEFAULT:0.7)." << endl
			 << "   --start  (int)   Start aligning at this read." << endl
			 << "   --stride (int)   Read stride (for multi-job alignment of the same file)." << endl
			 << "   -d 	(flag)  Enable dotPlot" << endl
			 << "   -PAl (int) Print at most how many alignments for one read" << endl
			 << "   -Al (int) Compute at most how many alignments for one read" << endl
			 << "   --passthrough Pass auxilary tags from the input unaligned bam to the output" << endl;
	cout << "Examples: " << endl
			 << "Aligning CCS reads:  lra align -CCS -t 16 ref.fa input.fasta/input.bam/input.sam -p s > output.sam" << endl
			 << "Aligning CLR reads:  lra align -CLR -t 16 ref.fa input.fasta/input.bam/input.sam -p s > output.sam" << endl
			 << "Aligning Nanopore reads:  lra align -ONT -t 16 ref.fa input.fasta/input.bam/input.sam -p s > output.sam" << endl;

}

class MapInfo {
public:
	vector<float> *LookUpTable;
	Header*header;
	Genome *genome;
	vector<GenomeTuple> *genomemm;
	LocalIndex *glIndex;
	Input *reader;
	Options *opts;
	ostream *out;
	ostream *svsigOut;
	int thread;
	int numThreads;
	pthread_mutex_t *semaphore;
	int *numAligned;
	int *numRead;
	Timing timing;
};

void MapReads(MapInfo *mapInfo) {
	Read read;
	stringstream strm;
	stringstream svsigstrm;
	vector<Read> reads;
	IndelRefineBuffers indelRefineBuffers;
	while (mapInfo->reader->BufferedRead(reads,IO_BUFFER_SIZE,(*mapInfo->opts))) {
		if (mapInfo->opts->readStride != 1 and
				mapInfo->reader->nReads % mapInfo->opts->readStride != mapInfo->opts->readStart ) {
			continue;
		}
		else {
			for (int i = 0; i< reads.size(); i++) {
				*mapInfo->numAligned+=MapRead(*mapInfo->LookUpTable, reads[i], *mapInfo->genome, *mapInfo->genomemm,
																			*mapInfo->glIndex, *mapInfo->opts, &strm, &svsigstrm, mapInfo->timing, indelRefineBuffers, 
																			mapInfo->semaphore);
				reads[i].Clear();
			}
			reads.clear();
			//
			// Print progress to the screen.
			//
			if (mapInfo->reader->basesRead > 100000000) {
				if (mapInfo->semaphore != NULL) {
					pthread_mutex_lock(mapInfo->semaphore);
				}
				// Check to see if another thread hit this spot at the same time
				if (mapInfo->reader->basesRead > 100000000) {
					
					clock_t cur = clock();		
					cerr << "lra aligned " << *mapInfo->numAligned << " from " << mapInfo->reader->nReads << ", " 
						 << mapInfo->reader->totalRead /1000000<< "M bases (" << std::setprecision(4) <<  
						 ((float)(cur - mapInfo->reader->timestamp))/CLOCKS_PER_SEC  << "s)." << endl;
					mapInfo->reader->timestamp=cur;
					mapInfo->reader->basesRead = 0;
				}
			
				if (mapInfo->semaphore != NULL) {
					pthread_mutex_unlock(mapInfo->semaphore);
				}
			}
			if (strm.str().size() > IO_BUFFER_SIZE) {
				if (mapInfo->semaphore != NULL) {
					pthread_mutex_lock(mapInfo->semaphore);
				}
				*mapInfo->out << strm.str();
				*mapInfo->svsigOut << svsigstrm.str(); 
				if (mapInfo->semaphore != NULL) {
					pthread_mutex_unlock(mapInfo->semaphore);
				}
				strm.str("");
				strm.clear();
				svsigstrm.str("");
				svsigstrm.clear();
			}
				
		}
	}
	if (mapInfo->semaphore != NULL) {
		pthread_mutex_lock(mapInfo->semaphore);
	}
	*mapInfo->out << strm.str();
	*mapInfo->svsigOut << svsigstrm.str();
	if (mapInfo->semaphore != NULL) {
		pthread_mutex_unlock(mapInfo->semaphore);
	}

	pthread_exit(NULL);
}

void RunAlign(int argc, const char* argv[], Options &opts ) {
// open query file for reading; you may use your favorite FASTA/Q parser
	string genomeFile = "", reads = "";
	string indexFile="";
	int w=10;
	vector<string> allreads;
	int argi = 0;
	for (argi = 0; argi < argc; ) {
		if (ArgIs(argv[argi], "-a")) {
			opts.storeAll=true;
		}				
		if (ArgIs(argv[argi], "-m")) {
			opts.minRefinedClusterSize = atoi(GetArgv(argv, argc, argi));
			++argi;
		}	
		else if (ArgIs(argv[argi], "-M")) {
			opts.minClusterSize = atoi(GetArgv(argv, argc, argi));
			++argi;
		}	
		else if (ArgIs(argv[argi], "-Flag")) {
			opts.flagRemove = atoi(GetArgv(argv, argc, argi));
			++argi;
		}	
		else if (ArgIs(argv[argi], "-H")) {
			opts.hardClip = true;
		}
		else if (ArgIs(argv[argi], "-p")) {			
			opts.printFormat = GetArgv(argv, argc, argi);
			++argi;
		}
		else if (ArgIs(argv[argi], "--maxDiag")) {
			opts.maxDiag = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "-t")) {
			opts.nproc = atoi(GetArgv(argv, argc, argi));
			++argi;
		}		
		else if (ArgIs(argv[argi], "-r")) {
			opts.refineLevel = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "-o")) {
			opts.outfile = GetArgv(argv, argc, argi);
			++argi;
		}
		else if (ArgIs(argv[argi], "--stride")) {
			opts.readStride = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--start")) {
			opts.readStart = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--PrintNumAln")) {
			opts.PrintNumAln = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--Al")) {
			opts.NumAln = atoi(GetArgv(argv, argc, argi));
			++argi;
		}		
		else if (ArgIs(argv[argi], "-S")) {
			opts.SparseDP = true;
		}
		else if (ArgIs(argv[argi], "--timing")) {
			opts.timing = GetArgv(argv, argc, argi);
			++argi;
		}			
		else if (ArgIs(argv[argi], "--timeRead")) {
			opts.storeTiming = true;
		}
		else if (ArgIs(argv[argi], "--passthrough")) {
			opts.passthroughtag = true;
		}
		else if (ArgIs(argv[argi], "--CalculateMinimizerStats")) {
			opts.CalculateMinimizerStats = true;
		}
    	else if (ArgIs(argv[argi], "--refineBand")) {
			opts.refineBand=atoi(GetArgv(argv, argc, argi));
			++argi;
		}
    	else if (ArgIs(argv[argi], "--anchorstoosparse")) {
			opts.anchorstoosparse=atof(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--skipBandedRefine")) {
			opts.skipBandedRefine = true;
		}
		else if (ArgIs(argv[argi], "-CONTIG")) {
			opts.readType=Options::contig;
			opts.HighlyAccurate = true;
			opts.anchor_rate=1.0; // boost the match score for 1st SDP
			opts.cleanMaxDiag=100; // For CleanOffDiagonal
			opts.maxDiag=100; // For StoreFineCluster
			opts.maxGap=500; // for StoreFineCluster; cannot be too large, otherwise lose lots of INV
			opts.RoughClustermaxGap=500;  // for SplitRoughCluster; cannot be too large, otherwise cannot find INV
			opts.NumAln=2; 
   			opts.PrintNumAln = 1;
   			opts.anchorstoosparse=0.01; // For Cluster Refinement
   			opts.minDiagCluster=10;
   			opts.minClusterLength=100;
   			opts.minClusterSize=30;
   			opts.firstcoefficient=24;
     		opts.merge_dist = 100;
     		opts.SecondCleanMinDiagCluster=40;
     		opts.minDiagCluster=30;
     		opts.refineSpaceDist=100000;
   			// opts.secondcoefficient=15;
   			// opts.rate_FirstSDPValue=0;
			// opts.rate_value=1;	
		}
		else if (ArgIs(argv[argi], "-CCS")) {
			opts.readType=Options::ccs;
			opts.HighlyAccurate = false;
			opts.anchor_rate=18.0;
			opts.NumAln=3;
   			opts.PrintNumAln = 1;
      		opts.merge_dist = 100;
      		opts.RoughClustermaxGap=500; 
     		opts.maxGap = 400;
     		opts.SecondCleanMinDiagCluster=40;
     		opts.minDiagCluster=10;
     		opts.cleanClustersize=100;
     		opts.punish_anchorfreq=10;
     		opts.anchorPerlength=10;
     		opts.refineSpaceDist=10000;

			// opts.rate_FirstSDPValue=0;
			// opts.rate_value=1;

			// opts.maxDiag=500;
			// opts.maxGap=5000;

			// opts.globalMaxFreq = 30;
			// opts.NumOfminimizersPerWindow = 1;	
			//opts.minClusterSize=5; 
			//opts.minClusterLength=50;  
			// opts.maxGapBtwnAnchors=1500;
		}
		else if (ArgIs(argv[argi], "-CLR")) {
			opts.HighlyAccurate = false;
			opts.anchor_rate=8.0;
			opts.NumAln=3;
   			opts.PrintNumAln = 1;
     		opts.merge_dist = 100;
     		opts.RoughClustermaxGap = 500; 
     		opts.maxGap = 500;
     		opts.SecondCleanMinDiagCluster=20;
     		opts.minDiagCluster=10;
      		opts.refineSpaceDist=10000;
    	

    		opts.minDiagCluster = 10;
    		opts.cleanMaxDiag = 30;
    		opt.RemovePairedIndels = false;
    		opts.RemoveSpuriousAnchors = false;
    		opts.bypassClustering = true;
    		opts.anchor_rate = 6;
    		opts.SecondCleanMinDiagCluster = 5;
    		opts.punish_anchorfreq = 5;
    		opts.anchorPerlength = 5;
    		opts.cleanClustersize = 100;
    		opts.SecondCleanMaxDiag = 30;
    		opts.maxGap = 1000;
    		opts.minClusterSize = 2;
    		opts.anchorstoosparse = 0.02;
			// opts.rate_FirstSDPValue=0;
			// opts.rate_value=1;		
			// opts.NumAln=3;
   			// opts.PrintNumAln = 1;
			// opts.rate_FirstSDPValue=0;
			// opts.rate_value=1;	
			// opts.maxDiag=800;
			// opts.maxGap=5000;
			// opts.globalMaxFreq = 50;
			// opts.NumOfminimizersPerWindow = 1;
			//opts.minClusterSize=5; 
			//opts.minClusterLength=50; 
			// opts.maxGapBtwnAnchors=1800;
		}		
		else if (ArgIs(argv[argi], "-ONT")) {
			opts.HighlyAccurate = false;
			opts.anchor_rate=8.0;
			opts.NumAln=3;
   			opts.PrintNumAln = 1;
   			opts.merge_dist = 30;
   			opts.RoughClustermaxGap=5000; 
     		opts.maxGap = 500;
     		opts.SecondCleanMinDiagCluster=30;
     		opts.minDiagCluster=10;
      		opts.refineSpaceDist=10000;
    		
			// opts.rate_FirstSDPValue=0;
			// opts.rate_value=1;	
			// opts.HighlyAccurate = false;
			// opts.anchor_rate=8.0;
			// opts.maxDiag=800;
			// opts.maxGap=5000;
			// opts.globalMaxFreq = 50;
			// opts.NumOfminimizersPerWindow = 1;
			// opts.maxGapBtwnAnchors=1800;
		}
		else if (ArgIs(argv[argi], "--bypassClustering")) {
			opts.bypassClustering = true;
			opts.NumAln=2;
		}
		else if (ArgIs(argv[argi], "--debug")) {
			opts.debug = true;
		}
		else if (ArgIs(argv[argi], "--read")) {
			opts.readname = string(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--firstcoefficient")) {
			opts.firstcoefficient = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--secondcoefficient")) {
			opts.secondcoefficient = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--merge_dist")) {
			opts.merge_dist = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--minDiagCluster")) {
			opts.minDiagCluster = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--SecondCleanMinDiagCluster")) {
			opts.SecondCleanMinDiagCluster = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--anchorPerlength")) {
			opts.anchorPerlength = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--punish_anchorfreq")) {
			opts.punish_anchorfreq = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--cleanClustersize")) {
			opts.cleanClustersize = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--RoughClustermaxGap")) {
			opts.RoughClustermaxGap = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--CheckTrueIntervalInFineCluster")) {
			opts.CheckTrueIntervalInFineCluster = true;
		}
		else if (ArgIs(argv[argi], "--anchor_rate")) {
			opts.anchor_rate = atof(GetArgv(argv,argc,argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--second_anchor_rate")) {
			opts.second_anchor_rate = atof(GetArgv(argv,argc,argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--maxGap")) {
			opts.maxGap= atoi(GetArgv(argv,argc,argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "-o")) {
			opts.outfile = argv[++argi];
		}
		else if (ArgIs(argv[argi], "--refineEnd")) {
			opts.refineEnd = true;
		}	
		else if (ArgIs(argv[argi], "-d")) {
			opts.dotPlot = true;
		}		
		else if (ArgIs(argv[argi], "--SkipRemovePairedIndels")) {
			opts.RemovePairedIndels = false;
		}
		else if (ArgIs(argv[argi], "--SkipRemoveSpuriousAnchors")) {
			opts.RemoveSpuriousAnchors = false;
		}
		else if (ArgIs(argv[argi], "--locBand")) {
			opts.localBand = atoi(GetArgv(argv,argc,argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--sseBand")) {
			opts.sseBand = atoi(GetArgv(argv,argc,argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--locIndel")) {
			opts.localIndel = atoi(GetArgv(argv,argc,argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--SecondCleanMaxDiag")) {
			opts.SecondCleanMaxDiag = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--cleanMaxDiag")) {
			opts.cleanMaxDiag = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--minClusterSize")) {
			opts.minClusterSize = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--minClusterLength")) {
			opts.minClusterLength = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--minDiagCluster")) {
			opts.minDiagCluster = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else if (ArgIs(argv[argi], "--maxCandidates")) {
			opts.maxCandidates = atoi(GetArgv(argv, argc, argi));
			++argi;
		}
		else {
			if (genomeFile == "") {
				genomeFile = argv[argi];
			}
			else {
				allreads.push_back(string(argv[argi]));
			}
		}
		++argi;
	}

	if (genomeFile == "" || allreads.size() == 0) {
		HelpMap();
		exit(1);
	}
	if (indexFile == "") {
		indexFile = genomeFile + ".mmi";
	}
	Header header;
	vector<GenomeTuple> genomemm;
	LocalIndex glIndex(opts.localIndexWindow);

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
	reader.Initialize(allreads);
	reader.flagRemove = opts.flagRemove;
	int offset=0;
	Read read;
	ostream *outPtr;
	ostream *outSVsig;
	ofstream outfile;
	ofstream outsvfile;
	if (opts.outfile == "" or opts.outfile=="-") {
		outPtr = &cout;
	}
	else {
		outfile.open(opts.outfile.c_str());
		outPtr = &outfile;
	}

	if (opts.outsvfile== "" or opts.outsvfile =="-") {
		outSVsig = &cout;
	}
	else {
		outsvfile.open(opts.outsvfile.c_str());
		outSVsig = &outsvfile;
	}

	if (opts.printFormat == "s") {
		stringstream cl;
		cl << "lra align";
		for (int i=0; i < argc; i++) {
			cl << " " << argv[i];
		}
		*outPtr << "@PG\tID:lra\tPN:lra\tVN:"<<lraVersion<<"\tCL:"<<cl.str() << endl;
		genome.header.WriteSAMHeader(*outPtr);
	}

	vector<float> LookUpTable;
	CreateLookUpTable(LookUpTable);

	if (opts.nproc > 1) {
		pthread_attr_t *threadAttr = new pthread_attr_t[opts.nproc];
		for (int procIndex = 0; procIndex < opts.nproc; procIndex++ ){
			pthread_attr_init(&threadAttr[procIndex]);
		}

		pthread_t *threads = new pthread_t[opts.nproc];
		vector<MapInfo> mapInfo(opts.nproc);		
		pthread_mutex_t semaphore;		
		pthread_mutex_init(&semaphore, NULL);
		int numAligned=0;
		for (int procIndex = 0; procIndex < opts.nproc; procIndex++){ 
			mapInfo[procIndex].LookUpTable = &LookUpTable;
			mapInfo[procIndex].genome = &genome;
			mapInfo[procIndex].genomemm = &genomemm;
			mapInfo[procIndex].glIndex = &glIndex;
			mapInfo[procIndex].reader = &reader;
			mapInfo[procIndex].opts= &opts;
			mapInfo[procIndex].out = outPtr;
			mapInfo[procIndex].svsigOut = outSVsig;
			mapInfo[procIndex].thread=procIndex;
			mapInfo[procIndex].semaphore=&semaphore;
			mapInfo[procIndex].numAligned=&numAligned;
			pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*))MapReads, &mapInfo[procIndex]);
		}

		for (int procIndex = 0; procIndex < opts.nproc; procIndex++) {
			pthread_join(threads[procIndex], NULL);
		}
		
		for (int procIndex = 1; procIndex < opts.nproc; procIndex++) {
			mapInfo[0].timing.Add(mapInfo[procIndex].timing);
		}
		if (opts.timing != "") {
			mapInfo[0].timing.Summarize(opts.timing);
		}
	}
	else {
		Timing timing;
		IndelRefineBuffers indelRefineBuffers;
		while (reader.GetNext(read, opts)) {
			int rstmm = 0;			

			MapRead(LookUpTable, read, genome, genomemm, glIndex, opts, outPtr, outSVsig, timing, indelRefineBuffers);
			if (opts.timing != "") {
				timing.Summarize(opts.timing);
			}
		}
	}
	outfile.close();
	outsvfile.close();
}

void HelpStoreIndex() {
	cout << "Usage: lra index file.fa [options]" << endl
			 << "  Global index options " << endl
			 << "	-CCS (flag) Index for aligning CCS reads" << endl
			 << "	-CLR (flag) Index for aligning CLR reads" << endl
			 << "	-ONT (flag) Index for aligning Nanopore reads" << endl
			 << "   -CONTIG (flag) Index for aligning large contigs" << endl
			 << "   -W (int) Minimizer window size (10)." << endl
			 << "   -F (int) Maximum minimizer frequency (DEFAULT: 80)." << endl
			 << "   -K (int) Word size" << endl
			 << "  Local index options: "<< endl
			 << "   -w (int) Local minimizer window size (10)." << endl
			 << "   -f (int) Local maximum minimizer frequency (5)." << endl
			 << "   -k (int) Local word size (10)" << endl
			 << "   -h Print help." << endl;
}

void HelpStoreGlobal() {
	cout << "Usage: lra index file.fa [options]" << endl;
	cout << "Options: " << endl
			 << "   -CCS (flag) Index for aligning CCS reads" << endl
			 << "   -CLR (flag) Index for aligning CLR reads" << endl
			 << "   -ONT (flag) Index for aligning Nanopore reads" << endl
			 << "   -CONTIG (flag) Index for aligning large contigs" << endl
			 << "   -W (int) Minimizer window size (10)." << endl
			 << "   -F (int) Maximum minimizer frequency. (default: 60 for CLR and NANO reads; 50 for CCS reads)" << endl
			 << "   -K (int) Word size" << endl
			 << "   -h Print help." << endl;	
	cout << "Examples: " << endl
			 << "Index reference for aligning CCS reads: lra index -CCS ref.fa" << endl
			 << "Index reference for aligning CLR reads: lra index -CLR ref.fa" << endl
			 << "Index reference for aligning Nanopore reads: lra index -ONT ref.fa" << endl
			 << "Index reference for aligning contig: lra index -CONTIG ref.fa" << endl;
}

void HelpStoreLocal() {
	cout << "Usage: lra local file.fa [options]" << endl
			 << "   -w (int) Local minimizer window size (10)." << endl
			 << "   -f (int) Local maximum minimizer frequency (5)." << endl
			 << "   -k (int) Local word size (10)" << endl
			 << "   -h Print help." << endl;
}

void RunStoreLocal(int argc, const char* argv[], LocalIndex &glIndex, Options &opts) {
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
		else if (ArgIs(argv[argi], "-CCS") or ArgIs(argv[argi], "-CONTIG") or ArgIs(argv[argi], "-CLR") or ArgIs(argv[argi], "-ONT") 
			or ArgIs(argv[argi], "--CalculateMinimizerStats")) {
			argi+=1;
			continue;
		}	
		else if (ArgIs(argv[argi], "-k")) {
			argi++;
			opts.localK = atoi(argv[argi]);
			glIndex.k = opts.localK;
		}
		else if (ArgIs(argv[argi], "-w")) {
			argi++;
			opts.localW=atoi(argv[argi]);
			glIndex.w=opts.localW;
		}
		else if (ArgIs(argv[argi], "-f")) {
			argi++;
			opts.localMaxFreq=atoi(argv[argi]);
			glIndex.maxFreq=opts.localMaxFreq;
		}
		else if (ArgIs(argv[argi], "-K") or ArgIs(argv[argi], "-W") or ArgIs(argv[argi], "-F")
				 or ArgIs(argv[argi], "-N") or ArgIs(argv[argi], "--globalWinsize")) {
			argi+=2;
			continue;
		}
		else if (ArgIs(argv[argi], "-d")) {
			argi+=1;
			continue;
		}
		else if (ArgIs(argv[argi], "--localIndexWindow")) {
			opts.localIndexWindow=atoi(GetArgv(argv, argc, argi));
			glIndex.StoreLocalIndexWindow(opts.localIndexWindow);
			++argi;		
		}
		else if (strlen(argv[argi]) > 0 and argv[argi][0] == '-') {
			HelpStoreLocal();
			cout << "Invalid option " << argv[argi] << endl;
			exit(1);
		}
		else {
			genome = argv[argi];
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

void RunStoreGlobal(int argc, const char* argv[], vector<GenomeTuple> &minimizers, Header &header, Options &opts) {
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
		else if (ArgIs(argv[argi], "-K")) {
			++argi;
			opts.globalK = atoi(argv[argi]);
		}
		else if (ArgIs(argv[argi], "-CONTIG")) {
			opts.globalK = 17;
			opts.globalW = 10;
			opts.globalMaxFreq = 60;
			opts.globalWinsize = 12;
			opts.NumOfminimizersPerWindow = 1;	
		}
		else if (ArgIs(argv[argi], "-CCS")) {
			opts.globalK = 15;
			opts.globalW = 10;
			opts.globalMaxFreq = 200;
			opts.globalWinsize = 9;
			opts.NumOfminimizersPerWindow = 1;			
		}	
		else if (ArgIs(argv[argi], "-CLR")) {
			opts.globalK = 15;
			opts.globalW = 10;
			opts.globalMaxFreq = 200;
			opts.globalWinsize = 9;
			opts.NumOfminimizersPerWindow = 1;	

			opts.globalK = 12;
			opts.globalW = 10;
			opts.globalWinsize = 7;

		}
		else if (ArgIs(argv[argi], "-ONT")) {
			opts.globalK = 15;
			opts.globalW = 10;
			opts.globalMaxFreq = 200;
			opts.globalWinsize = 9;
			opts.NumOfminimizersPerWindow = 1;			
		}
		else if (ArgIs(argv[argi], "-F")) {
			opts.globalMaxFreq=atoi(GetArgv(argv, argc, argi));
			++argi;
		}	
		else if (ArgIs(argv[argi], "-N")) {
			opts.NumOfminimizersPerWindow=atoi(GetArgv(argv, argc, argi));
			++argi;
		}			
		else if (ArgIs(argv[argi], "--globalWinsize")) {
			opts.globalWinsize=atoi(GetArgv(argv, argc, argi));
			++argi;
		}					
		else if (ArgIs(argv[argi], "-d")) {
			opts.dotPlot = true;
		}
		else if (ArgIs(argv[argi], "--SkipLocalMinimizer")) {
			opts.SkipLocalMinimizer = true;
		}
		else if (ArgIs(argv[argi], "--SkipClusering")) {
			opts.SkipClusering = true;
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
			printIndex = true;
		}	
		else if (ArgIs(argv[argi], "--CalculateMinimizerStats")) {
			opts.CalculateMinimizerStats = true;
		}
		else if (ArgIs(argv[argi], "-k") or ArgIs(argv[argi], "-w") or ArgIs(argv[argi], "-f")) {
			argi+=2;
			continue;
		}
		else if (ArgIs(argv[argi], "--localIndexWindow")) {
			opts.localIndexWindow=atoi(GetArgv(argv, argc, argi));
			++argi;		
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
	cout << "Program: lra (long sequence alignment)" << endl;
	cout << "Version: " << lraVersion << endl;
	cout << "Contact: Mark Chaisson (mchaisso@usc.edu) and Jingwen Ren (jingwenr@usc.edu)" << endl << endl;
	cout << "Usage:   lra <command> [options]"<< endl << endl;
	cout << "Command: index   - Build global and local indexes on a genome." << endl;
	cout << "         align   - Map reads using the index." << endl;
	cout << "         global  - Build a global index." << endl;
	cout << "         local   - Build local index." << endl;
}

void InitStatic() {
	Tuple mask = 1;
	GenomeTuple::for_mask_s = ~(mask << (sizeof(mask)*8-1)); // for_mask_s = 0111...11
	Tuple rmask = 1;
	GenomeTuple::rev_mask_s = (rmask << (sizeof(rmask)*8-1)); // rev_mask_s = 1000...00
	LocalTuple::for_mask_s = 1;
	for (int i = 1; i < 32 - LOCAL_POS_BITS; i++) {
		LocalTuple::for_mask_s = LocalTuple::for_mask_s << 1;
		LocalTuple::for_mask_s += 1; // for_mask_s = 111...11 (20 bit)
	}
	LocalTuple::rev_mask_s = 0;
	// CreateLookUpTable(LookUpTable);
	// cerr << "GenomeTuple::for_mask_s: " << GenomeTuple::for_mask_s << endl;
	// cerr << "LocalTuple::for_mask_s: " << LocalTuple::for_mask_s << endl;
	// u_int32_t local_mask = 1;
	// for (int i = 0; i < LOCAL_POS_BITS; i++) {
	// 	LocalTuple::for_mask_s = LocalTuple::for_mask_s << 2;
	// 	LocalTuple::for_mask_s += 3; // for_mask_s = 111...11 (20 bit)
	// }
}

int main(int argc, const char *argv[]) {
	if (argc < 2) {
		Usage();
		return 1;
	}
	Options opts;
	InitStatic();
  	int argi;
	vector<GenomeTuple> minimizers;
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
