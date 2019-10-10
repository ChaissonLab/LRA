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
	int flagRemove;
	int minRemovePairedIndelsLength; // if an anchor's length is larger than this parameter, 
									// then even if it has paired indels before and after it, we do not delete this anchor.
	int maxRemoveSpuriousAnchorsDist; 
	int minRemoveSpuriousAnchorsNum;
	int minRemoveSpuriousAnchorsLength;
	int SecondaryAln;
	int PrimaryAln;
	int BtnSubClusterswindow;

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
		maxDiag=500; // We want maxDiag to be a small number 
		cleanMaxDiag=50; // used to be 50  
		minDiagCluster=10; 	// used to be 20
							// This parameter is used in CleanOffDiagonal function; It's better not to set it to a single value. 
							// This parameter is used in another CleanOFFDiagonal function
							// This parameter can be deleted here

		minClusterSize=5; // used to be 10 
		minClusterLength=50;  // used to be 100
		minRefinedClusterSize=40;
		window=0;
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
		maxGapBtwnAnchors=20000; // no larger than 2000 // used to be 1500
		mergeClusters=true;
		NaiveDP=false;
		seqan=false;
		SparseDP=true;
		LookUpTable=true;
		MergeSplit=true;
   		flagRemove=0;
   		minRemovePairedIndelsLength=50; // used to be 50
   		maxRemoveSpuriousAnchorsDist=200;
   		minRemoveSpuriousAnchorsNum=15;
   		minRemoveSpuriousAnchorsLength=100;
   		SecondaryAln = 0;
   		PrimaryAln = 1;
   		BtnSubClusterswindow = 800;

	}
};
#endif
