#ifndef OPTIONS_H_
#define OPTIONS_H_

const unsigned int REF_LOC=1;
const unsigned int REF_DYN=2;
const unsigned int REF_DP=4;

class Options {
public:
	int globalK;
	int localK;
	int globalW;
	int localW;
	int globalMaxFreq;
	int localMaxFreq;
	int maxDiag;
	int cleanMaxDiag;
	int minClusterSize;
	int minClusterLength;
	int window;
	bool dotPlot;
	bool mergeClusters;
	bool mergeGapped;
	int minDiagCluster;
	int minRefinedClusterSize;
	bool viewPairwise;
	bool hardClip;
	char printFormat;
	int bestn;
	bool storeAll;
	int nproc;
	string outfile;
	int maxCandidates;
	int refineLevel;
	bool doBandedAlignment;
	int maxGap;
	int maxGapBtwnAnchors;
	bool NaiveDP;
	bool SparseDP;
	bool LookUpTable;
	int readStart;
	int readStride;
	bool seqan;
	int localMatch;
 	int localMismatch;
	int localIndel;
	int localBand;
	int MergeSplit;
	int mintupChainClustersize;
	int flagRemove;
	int minRemovePairedIndelsLength; // if an anchor's length is larger than this parameter, 
									// then even if it has paired indels before and after it, we do not delete this anchor.
	int maxRemoveSpuriousAnchorsDist; 
	int minRemoveSpuriousAnchorsNum;
	int minRemoveSpuriousAnchorsLength;

	Options() {
		localMatch=4;
		localMismatch=-3;
		localIndel=-3;
		localBand=15;
		readStart=0;
		readStride=1;
		dotPlot=false;
		globalK=17;
		globalW=10; 
		localK=7;
		localW=5;
		bestn=1;
		globalMaxFreq=20;
		localMaxFreq=30;
		maxDiag=500; // We want maxDiag to be a small number // used to be 500
		cleanMaxDiag=100; // used to be 50
		minDiagCluster=20;
		minClusterSize=10;
		minClusterLength=100;
		minRefinedClusterSize=40;
		window=100;
		mergeGapped=false;
		viewPairwise=false;
		hardClip=false;
		printFormat='b';
		storeAll=false;
		nproc=1;
		outfile="";
		maxCandidates=3;
		doBandedAlignment=true;
		refineLevel= REF_LOC | REF_DYN | REF_DP;
		maxGap=10000; 
		maxGapBtwnAnchors=1500; // no larger than 2000 // used to be 2000
		mergeClusters=false;
		NaiveDP=false;
		seqan=false;
		SparseDP=true;
		LookUpTable=true;
		MergeSplit=true;
		mintupChainClustersize=10;
   		flagRemove=0;
   		minRemovePairedIndelsLength=400;
   		maxRemoveSpuriousAnchorsDist=500;
   		minRemoveSpuriousAnchorsNum=20;
   		minRemoveSpuriousAnchorsLength=100;
	}
};
#endif
