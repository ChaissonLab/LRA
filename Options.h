#ifndef OPTIONS_H_
#define OPTIONS_H_
class Options {
public:
	int k;
	int w;
	int maxFreq;
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
	Options() {
		k=15;
		w=10; 
		bestn=1;
		maxFreq=200;
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
	}
};
#endif
