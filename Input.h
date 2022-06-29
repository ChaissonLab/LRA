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
#include <algorithm>
#include <cctype>
#include "Options.h"

KSEQ_INIT(gzFile, gzread)

class Input {
public:
  enum InputType { FASTA, FASTQ, HTS };
  int inputType;
  istream *strmPtr;
  ifstream strm;
  gzFile fp;
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
  string format;
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
    fp=NULL;
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
    if (s.peek() != '@') { return false;}
    getline(s,lines[0]);
    getline(s,lines[1]);
    bool res=false;
    if (s.peek() == '+') {
      res=true;
    }
    for (int j=2; j > 0; j--) {
      s.putback('\n');
      for (int i=lines[j-1].size(); i > 0; i--) {
	s.putback(lines[j-1][i-1]);
      }
    }
    return res;
  }
  
  bool Initialize(string &filename) {
    nReads=0;

    doInit = false;
		
    /*
    // When HTSLIB 1.13 is released, this may be used to initialize input
    htsfp = hts_open(filename.c_str(),"r");
    const htsFormat *fmt = hts_get_format(htsfp);
    format=hts_format_file_extension(fmt);
		
    if (format == "fq" or format == "fa") {
    inputType=0;
    return true;
    }
    else if (format == "sam" or format == "bam" or format == "cram") {
    samHeader = sam_hdr_read(htsfp);
    return true;
    }
    else {
    cerr << "Cannot determine type of input " << endl;
    exit(1);
    }
    */

    if (filename == "-" or filename == "stdin" or filename == "/dev/stdin") {
      strmPtr = &cin;
    }
    else {
      strm.open(filename.c_str());
      strmPtr=&strm;
    }
		   
    if (StreamIsFasta(*strmPtr)) {
      inputType=FASTA;
      return true;		  
    }
    else if (StreamIsFastq(*strmPtr) ) {
      inputType=FASTQ;
      return true;		  		
    }
    else {
		
      //
      // possibly sam 
      //
      if (filename == "-" or filename == "stdin" or filename == "/dev/stdin") {
	cout << "Streaming of sam/bam/cram input is not supported. You can convert to fasta/fastq, e.g.:" << endl
	     << "samtools fasta input.bam  | lra align ref.fa -" << endl;
	exit(1);
      }
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
      inputType=HTS;
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
    bool readOne=false;
    if (overrideSemaphore == false and top == true) {
      pthread_mutex_lock(&semaphore);
    }
    string name;
    string seq;
    string qual;
	  //----------------------------------------------
             if (inputType == FASTA and strmPtr->eof())  // Any more FASTA files?
             {
                   strm.close(); // at eof so close before checking if another file is in input list
                   ++curFile;

                   if (curFile >= allReads.size()) // any more input file in list
                   { // no more read files?
                         return 0;
                   }
                   if (Initialize(allReads[curFile]) == false) // does next input file initialise?
                   {
                         return 0;
                   }
             }	  
	  //----------------------------------------------
    if (inputType == FASTA or inputType == FASTQ) {
      if (strmPtr->eof()) {
	     return 0;
      }
      if (inputType == FASTA) {
      	string header;
      	char c;
      	getline(*strmPtr, header);
      	stringstream nameStrm(header);
      	nameStrm >> c >> read.name;
      	c=strmPtr->peek();

      	while (c != EOF and c != '>') {
      	  string line;
      	  getline(*strmPtr, line);
      	  int i=0,j=0;
      	  for (i=0; i < line.size(); i++) { if (line[i] != ' ') { line[j] = toupper(line[i]); j++;} }
      	  line.resize(j);		      
      	  seq+=line;
      	  c=strmPtr->peek();
      	}
      	if (c == EOF) {
      	  strmPtr->get();
      	}
      }

      else if (inputType == FASTQ) {
      	string header;
	string sep;
      	char c;
      	getline(*strmPtr, header);
	getline(*strmPtr, seq);
	getline(*strmPtr, sep);
	getline(*strmPtr, qual);
	if (header.size() ==0 or seq.size() == 0 or sep.size() == 0 or qual.size() == 0) {
		// -------------------------------------------
                  strm.close();
                  ++curFile;
                  if (curFile >= allReads.size())    // Exit if no more input files.
                  { // no more read files?
                        readOne = false;
                        return 0;
                  }
                  if (Initialize(allReads[curFile]) == false)
                  {
                       readOne = false;
                       return 0;
                  }
                  // opened next input file - is it fastq? (set in Initialize())
                  if (inputType == FASTQ)
                  {
                       getline(*strmPtr, header);
                       getline(*strmPtr, seq);
                       getline(*strmPtr, sep);
                       getline(*strmPtr, qual);
                  }
        }		
        if (header.size() == 0 or seq.size() == 0 or sep.size() == 0 or qual.size() == 0)
        {
            readOne = false;
            return 0;
        }		
		// -------------------------------------------
	else {
	  stringstream nameStrm(header);
	  nameStrm >> c >> read.name;
	  int i,j;
	  for (i=0,j=0; i < seq.size(); i++) { if (seq[i] != ' ') { seq[j] = toupper(seq[i]); j++;} }
	  seq.resize(j);

	  for (i=0,j=0; i < qual.size(); i++) { if (qual[i] != ' ') { qual[j] = qual[i]; j++;} }
	  qual.resize(j);
	}
      }
      read.seq = new char[seq.size()+1];
      memcpy(read.seq, seq.c_str(), seq.size());
      read.length=seq.size();
      read.seq[read.length] = '\0';
      if (qual.size() > 0) {
      	assert(qual.size() == seq.size());
      	read.qual = new char[qual.size()+1];
      	memcpy(read.qual, qual.c_str(), qual.size());
	read.qual[read.length] = '\0';
      }
      read.passthrough=NULL;
      readOne=true;
      nReads++;
    }
    else if (inputType == HTS) {
      int res;
      bam1_t *b = bam_init1();
      res= sam_read1(htsfp, samHeader, b);

      while (res >= 0 and readOne == false) {
      	 if (res >= 0) {	
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
      		  if (lenPassthrough > 0) {											
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
      	    nReads++;
      	    bam_destroy1(b);
      	    b=NULL;
      	    //bam1_t *b = bam_init1();
      	  }
      	  else {
      	    bam_destroy1(b);
      	    b = bam_init1();
      	    res= sam_read1(htsfp, samHeader, b);
      	  }
      	}
      }
      if (res < 0) {	
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

