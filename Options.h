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
	int minClusterSize;
	int window;
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
	bool NaiveDP;

	
	Options() {
		globalK=17;
		globalW=10; 
		localK=7;
		localW=5;
		bestn=1;
		globalMaxFreq=20;
		localMaxFreq=30;
		maxDiag=3000;
		minDiagCluster=2;
		minClusterSize =20;
		minRefinedClusterSize = 40;
		window=1000;
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
		mergeClusters=false;
		NaiveDP = false;
	}
};
#endif
