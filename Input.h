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


	bool Initialize(string &filename) {
		//
		// Check to see if the input is fasta
		//
	

		semaphore = sem_open("/reader",     O_CREAT, 0644, 1);
		sem_init(semaphore, 0, 1);

		istream *strmPtr;
		if (filename == "-") {
			strmPtr = &std::cin;
		}
		else {
			strm.open(filename.c_str());
			strmPtr = &strm;
		}

		if (strmPtr->peek() != EOF && strmPtr->peek() == '>') {
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
		}
	}

	bool GetNext(Read &read) {
		sem_wait(semaphore);
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
			
#define bam_get_seq(b)   ((b)->data + ((b)->core.n_cigar<<2) + (b)->core.l_qname)
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
		sem_post(semaphore);
		return readOne;
	}

};


#endif
