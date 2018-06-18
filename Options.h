#ifndef OPTIONS_H_
#define OPTIONS_H_
class Options {
public:
	int k;
	int w;
	int maxFreq;
	Options() {
		k=15; w=10; maxFreq=200;
	}
};
#endif
