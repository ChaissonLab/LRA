#ifndef INPUT_H_
#define INPUT_H_

#include <iostream>
#include <stdlib.h>

#include <pthread.h>
#include <semaphore.h>

class Input {
 public:
	int inputType;
	istream *strmPtr;
	ifstream strm;
	gzFile fastaFile;
	htsFile *htsfp;
	kseq_t *ks;
	bam_hdr_t *samHeader;			
	sem_t *semaphore;
	bool doInit;
	int curFile;
	int basesRead;
	long totalRead;
	vector<string> allReads;
	Input() {
		doInit = true;
		curFile = 0;
		basesRead = 0;
		totalRead = 0;
	}
	bool StreamIsFasta(istream &s) {
		if (s.eof() or s.good() == false) {
			return false;
		}
		if (s.peek() == '>') {
			return true;
		}
		return false;
	}

	bool StreamIsFastq(istream &s) {

		if (s.eof() or s.good() == false) {
			return false;
		}
		vector<string> lines(2);
		string line;
		getline(s,lines[0]);
		getline(s,lines[1]);
		if (s.peek() == '+') {
			s.putback('\n');
			for (int j=2; j > 0; j--) {
				for (int i=lines[j-1].size(); i > 0; i--) {
					s.putback(lines[j-1][i]);
				}
			}
			return true;
		}
		return false;
	}
	bool Initialize(string &filename) {
		istream *strmPtr;
		if (filename == "-") {
			strmPtr = &std::cin;
		}
		else {
			strm.close();
			strm.open(filename.c_str());
			strmPtr = &strm;
		}
		if (StreamIsFasta(*strmPtr) or StreamIsFastq(*strmPtr) ) {

			inputType=0;
			fastaFile = gzopen(filename.c_str(), "r");
			ks = kseq_init(fastaFile);
			return true;
		}
		else {
			//
			// possibly sam 
			//
			htsfp = hts_open(filename.c_str(),"r");
			const htsFormat *fmt = hts_get_format(htsfp);
			if (fmt == NULL or (fmt->format != sam and fmt->format != bam)) {
				cout << "Cannot determine format of input reads." << endl;
				exit(1);
			}

      samHeader = sam_hdr_read(htsfp);
			inputType=1;
			return true;
		}
		return false;
	}

	bool Initialize(vector<string> &_allReads) {
		//
		// Check to see if the input is fasta
		//
		allReads = _allReads;
		if (allReads.size() == 0) {
			exit(0);
		}

		semaphore = sem_open("/reader",     O_CREAT, 0644, 1);
		sem_init(semaphore, 0, 1);

		if (Initialize(allReads[curFile]) == false) {
			return 0;
		}
		doInit = false;
	}

	bool GetNext(Read &read, bool top=true) {
		if (top == true) {
			sem_wait(semaphore);
		}

		if (doInit) {
			if (curFile >= allReads.size() or 
					Initialize(allReads[curFile]) == false) {
				if (top) {
					doInit = false;
					sem_post(semaphore);
				}
				return 0;
			}
			doInit = false;
		}

		bool readOne=false;
		if (inputType == 0) {
			if (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
				read.seq = new char[ks->seq.l];
				read.length=ks->seq.l;
				memcpy(read.seq, ks->seq.s,read.length);
				read.name=string(ks->name.s);
				read.qual = ks->qual.s;
				read.passthrough=NULL;
				readOne=true;

			}
		}
		else if (inputType == 1) {
			int res;
			bam1_t *b = bam_init1();

			res= sam_read1(htsfp, samHeader, b);
			if (res > 0) {
				
				read.length = b->core.l_qseq;			
				read.seq = new char[read.length];
				read.name = string(bam_get_qname(b));
				uint8_t *q = bam_get_seq(b);
				for (int i=0; i < read.length; i++) {
					read.seq[i]=seq_nt16_str[bam_seqi(q,i)];
				}
				read.qual = NULL;
				
				// 
				// Eventually this will store the passthrough data
				//
				bam_destroy1(b);
				readOne=true;
			}
		}

		if (readOne == false) {
			++curFile;
			doInit=true;
			readOne=GetNext(read, false);
		}
		basesRead += read.length;
		totalRead += read.length;
		if (basesRead > 10000000) {
			cerr << "lra processed " << totalRead << endl;
			basesRead = 0;
		}
		sem_post(semaphore);
		return readOne;
	}

};


#endif
