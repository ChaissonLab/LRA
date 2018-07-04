#ifndef OPTIONS_H_
#define OPTIONS_H_
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
	Options() {
		globalK=15;
		globalW=10; 
		localK=7;
		localW=5;
		bestn=1;
		globalMaxFreq=200;
		localMaxFreq=30;
		maxDiag=3000;
		minDiagCluster=2;
		minClusterSize =5;
		minRefinedClusterSize = 30;
		window=1000;
		mergeGapped=true;
		viewPairwise=false;
		hardClip=false;
		printFormat='b';
		storeAll=false;
		nproc=1;
		outfile="";
		maxCandidates=4;
	}
};
#endif
