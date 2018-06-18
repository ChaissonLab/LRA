#ifndef OPTIONS_H_
#define OPTIONS_H_
class Options {
public:
	int k;
	int w;
	int maxFreq;
	int maxDiag;
	int minClusterSize;
	Options() {
		k=15; w=10; maxFreq=200;
		maxDiag=10000;
		minClusterSize =5;
	}
};
#endif
