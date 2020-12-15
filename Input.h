#ifndef INPUT_H_
#define INPUT_H_

#include <iostream>
#include <stdlib.h>

#include <pthread.h>
#include <semaphore.h>
#include <time.h>
#include "Read.h"
#include <fcntl.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <assert.h>
#include "Options.h"

KSEQ_INIT(gzFile, gzread)

class Input {
 public:
	int inputType;
	istream *strmPtr;
	ifstream strm;
	gzFile fastaFile;
	htsFile *htsfp;
	kseq_t *ks;
	bam_hdr_t *samHeader;			
	pthread_mutex_t semaphore;
	bool doInit;
	int curFile;
	int basesRead;
	long totalRead;
	bool done;
	int nReads;
	int flagRemove;
	clock_t timestamp;

	vector<string> allReads;
	Input() {
		flagRemove=0;
		doInit = true;
		curFile = 0;
		basesRead = 0;
		totalRead = 0;
		done=false;
		ks=NULL;
		htsfp=NULL;
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
		nReads=0;
		istream *strmPtr;
		doInit = false;
		if (filename == "-" or filename == "stdin") {
			strmPtr = &std::cin;
		}
		else {
			strm.close();
			strm.open(filename.c_str());
			strmPtr = &strm;
		}
		if (StreamIsFasta(*strmPtr) or StreamIsFastq(*strmPtr) ) {
			inputType=0;
			if (ks != NULL) {
				kseq_destroy(ks);
				gzclose(fastaFile);
			}
			if (filename == "-" or filename == "/dev/stdin" or filename == "stdin") {
				gzFile fp = gzdopen(fileno(stdin), "r");
				ks = kseq_init(fp);
			}
			else {
				fastaFile = gzopen(filename.c_str(), "r");
				ks = kseq_init(fastaFile);
			}
			return true;
		}
		else {
			//
			// possibly sam 
			//

			if (htsfp != NULL) {
				hts_close(htsfp);
				bam_hdr_destroy(samHeader);

			}	
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


		pthread_mutex_init(&semaphore, NULL);


		if (Initialize(allReads[curFile]) == false) {
			return 0;
		}
		timestamp = clock();
		doInit = false;
		return true;
	}

	bool GetNext(Read &read, Options &opt, bool overrideSemaphore=false, bool top=true) {
		read.Clear();

		if (overrideSemaphore == false and top == true) {
			pthread_mutex_lock(&semaphore);
		}
		++nReads;
		if (doInit) {
			assert(top == false);
			doInit = false;

			if (curFile >= allReads.size() or 
					 Initialize(allReads[curFile]) == false) {
				//
				// Do not continue to read.
				done=true;
				return false;
			}
			doInit = false;
		}

		bool readOne=false;
		if (done == false) {
			if (inputType == 0) {
				if (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
					read.seq = new char[ks->seq.l];
					read.length=ks->seq.l;
					memcpy(read.seq, ks->seq.s,read.length);
					for (int i=0;i<ks->seq.l;i++) { read.seq[i] = toupper(ks->seq.s[i]);}
					read.name=string(ks->name.s);
					if (ks->qual.s != NULL) {
						read.qual = new char[ks->seq.l];
						memcpy(read.qual, ks->qual.s, ks->seq.l);
					}
					else {
						read.qual = new char[2];
						read.qual[0] = '*';
						read.qual[1] = '\0';
					}					
					read.passthrough=NULL;
					readOne=true;
				}
			}
			else if (inputType == 1) { // sam input
				int res;
				bam1_t *b = bam_init1();
				res= sam_read1(htsfp, samHeader, b);
				while (res > 0 and readOne == false) {
					if (res > 0) {	
					if ((b->core.flag & flagRemove) == 0) {		
							// // get auxilary tags
							// if (opt.passthroughtag and bam_get_aux(b))	{
							// 	unsigned char *pq = bam_get_aux(b);
							// 	int pq_len = strlen((char*)pq);
							// 	read.passthrough = new unsigned char[pq_len + 1];
							// 	for (int p=0; p<pq_len; p++) {
							// 		read.passthrough[p] = pq[p];
							// 	}
							// 	read.passthrough[pq_len] = '\0';

							// }				
							read.length = b->core.l_qseq;			
							read.seq = new char[read.length];
							read.name = string(bam_get_qname(b));
							read.flags = b->core.flag;
							uint8_t *q = bam_get_seq(b);
							for (int i=0; i < read.length; i++) {read.seq[i]=seq_nt16_str[bam_seqi(q,i)];	}
							char* qual=(char*) bam_get_qual(b);
							if (qual[0] == char(0xff)) {
								read.qual = new char[2];
								read.qual[1] = '\0';
								read.qual[0] = '*';
							}
							else {
								read.qual=new char[read.length+1];
								for (int q=0; q < read.length; q++) {
									read.qual[q] = qual[q]+33;
								}
								read.qual[read.length]='\0';
							}
				
							// 
							// Eventually this will store the passthrough data
							//
							readOne=true;
							if (opt.passthroughtag) {
								int ksLen;
								kstring_t fullKs;
								int fullLen;
								fullKs = { 0, 0, NULL };
								fullLen = sam_format1(samHeader, b, &fullKs);
								int t=0;
								int numTab=0;							
								while (t < fullKs.l and numTab < 11) 
									{
										if (fullKs.s[t] == '\t') 
											{
												numTab++;
											}
										t+=1;									
									}
								if (t < fullKs.l) 
									{
										int lenPassthrough=fullKs.l-t;									
										if (lenPassthrough > 0) 
											{											
												read.passthrough=new char[lenPassthrough+1];
												read.passthrough[lenPassthrough]='\0';
												memcpy(read.passthrough, fullKs.s + t, lenPassthrough);											
											}
										else
											{
												read.passthrough=NULL;
											}
									}
								free(fullKs.s);								
							}
							bam_destroy1(b);
							//bam1_t *b = bam_init1();
						}
						else {
							bam_destroy1(b);
							b = bam_init1();
							res= sam_read1(htsfp, samHeader, b);
						}
					}
				}
				if (res == 0) {	
					if (b != NULL) {
						bam_destroy1(b);
						readOne = false;
					}
				}
			}

			if (readOne == false and top == true ) {
				++curFile;
				doInit=true;
				readOne=GetNext(read, opt, overrideSemaphore, false);
			}
			basesRead += read.length;
			totalRead += read.length;
		}

		if (overrideSemaphore== false and top == true) {
			pthread_mutex_unlock(&semaphore);
		}
		return readOne;
	}
	bool BufferedRead(vector<Read> &reads, int maxBufferSize, Options &opt) {
		int totalSize=0;

		pthread_mutex_lock(&semaphore);
		Read read;

		while(totalSize < maxBufferSize and GetNext(read, opt, true, true)) {
			reads.resize(reads.size()+1);
			reads[reads.size()-1]=read;
			totalSize += read.length;
			read.Clear();
		}

		pthread_mutex_unlock(&semaphore);

		return reads.size();
	}


};


#endif

