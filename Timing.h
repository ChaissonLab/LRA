#ifndef _TIMING_H_
#define _TIMING_H_

#include <vector>
#include <ostream>
#include <fstream>
#include <string>

class Timing {
public:
	vector<long> ticks;
	vector<string> labels;
	int index;
	clock_t curTime;
	clock_t startTime;
	void Start() {
		index=0;
		curTime=clock()/1000;
		startTime=curTime;
	}
	int Elapsed() {
		if (ticks.size() == 0) {
			return 0;
		}
		else {
			return clock()/1000 - startTime;
		}
	}
		
	void Tick(string label) {
		if (index +1 > ticks.size()) {
			ticks.push_back(0);
		}
		long prev=ticks[index];
		ticks[index] += clock()/1000 - curTime;
		if (prev > ticks[index]) {
			cerr << "Warning, tick has wrapped over long" << endl;
		}
		curTime=clock()/1000;
		if (index + 1 > labels.size()) {
			labels.push_back(label);
		}
		index+=1;
	}
	void Add(Timing &t) {
		//
		// Only add ticks the same size -- sometimes a thread 
		// will be created without adding ticks.
		//
		if (t.ticks.size() == ticks.size()) {
			for (int i =0; i < t.ticks.size(); i++) {
				ticks[i] += t.ticks[i];
			}
		}
	}
	void Summarize(const string &outFileName) {
		long total=0;
		ofstream outFile(outFileName.c_str());
		for (int i=0; i<ticks.size(); i++) {
			total+= ticks[i];
		}
		for (int i=0; i < ticks.size(); i++) {
			outFile << labels[i] << "\t" << ticks[i] << "\t" << ((double)ticks[i])/total << endl;
		}
		outFile.close();
	}
		
};

#endif
