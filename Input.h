#ifndef INPUT_H_
#define INPUT_H_

#include <iostream>
#include <stdlib.h>

class Input {
 public:
	int inputType;
	istream *strmPtr;
	ifstream strm;
	gzFile fastaFile;
	kseq_t *ks;
	bool Initialize(string &filename) {
		if (filename == "-") {
			strmPtr = &std::cin;
		}
		else {
			strm.open(filename.c_str());
			strmPtr = &strm;
		}
		if (SetInputType(strmPtr) == false) {
			cout << "ERROR. Could not determine type of input." << endl;
			exit(1);
		}
		if (inputType == 0) {
			//
			// Initialize for fasta
			//
			fastaFile = gzopen(filename.c_str(), "r");
			ks = kseq_init(fastaFile);
		}
		else if (inputType == 1) {
			//
			// Initialize for sam
			//
			cerr << "SAM input not finished." << endl;
			exit(1);
		}
	}

	bool GetNext(Read &read) {
		if (inputType == 0) {
			if (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
				read.seq=ks->seq.s;
				read.length=ks->seq.l;
				read.name=string(ks->name.s);
				read.qual = ks->qual.s;
				read.passthrough=NULL;
				return true;
			}
			else {
				return false;
			}
		}
	}
	
	bool SetInputType(istream *strmPtr) {
		if (strmPtr->peek() != EOF && strmPtr->peek() == '>') {
			inputType=0;
			return true;
		}
		else if (strmPtr->peek()!= EOF && strmPtr->peek() == '@') {
			inputType=1;
			return true;
		}
		else {
			return false;
		}
	}
		
			



};


#endif
