#ifndef ALIGNMENT_TYPE_H_
#define ALIGNMENT_TYPE_H_

#include <vector>


#include "AlignmentBlock.h"
#include "Path.h"
#include "SeqUtils.h"
#include "Options.h"
#include <assert.h>
#include <sstream>
#include <algorithm>
using namespace std;
const unsigned int READ_UNMAPPED=0x4;
const unsigned int READ_REVERSE=0x10;
const unsigned int READ_SECONDARY=0x100;
const unsigned int READ_MULTIPLESEGMENTS=0x1;
const unsigned int READ_FIRSTSEGMENT=0x40;
const unsigned int READ_LASTSEGMENT=0x80;
const unsigned int READ_SUPPLEMENTARY=0x800;

class Alignment {
 public:
	unsigned char mapqv;
	unsigned int flag;
	int chromIndex;
	string chrom;
	string name;
	string queryString, alignString, refString;
	int nblocks;
	int nm, nmm, nins, ndel;
	int preClip, sufClip;
	string cigar;
	bool prepared;
	char *read;
	char *forward;
	int readLen;
	int refLen;
	GenomePos genomeLen;
	string readName;
	char *genome;
	int strand;
	void Clear() {
		queryString=alignString=refString="";
		blocks.clear();
	}
	GenomePos qStart, qEnd, tStart, tEnd;

	void SetSecondary() {
		flag = flag | READ_SECONDARY;
	}
	void SetReverse() {
		flag = flag | READ_REVERSE;
	}
	void SetUnmapped() {
		flag = flag | READ_UNMAPPED;
	}

	Alignment() {
		flag=0;
		mapqv=0;
		nm=nmm=nins=ndel=0;
		nblocks=0;
		preClip=0; sufClip=0;
		prepared=false;
		read=NULL;
		forward=NULL;
	}
 	Alignment(char *_read, char *_forward, 
					 int _rl, string _rn, int _str, 
					 char *_genome, GenomePos _gl, string &_chrom, int _ci) : Alignment() { 
		read=_read; 
		forward=_forward;
		readLen = _rl; 
		readName= _rn;
		strand= _str;
		genome=_genome;
		genomeLen = _gl;
		chrom=_chrom;
		chromIndex=_ci;
	}

	int GetQStart() const {
		if (blocks.size() > 0) {
			return blocks[0].qPos;
		}		
		return 0;
	}
	int GetTStart() const {
		if (blocks.size() > 0) {
			return blocks[0].tPos;
		}
		return 0;
	}
	int qPos, tPos;
	// The positions in every block are relative to qPos and tPos;
	vector<Block> blocks;
	int size() const {
		return blocks.size();
	}

	void AppendAlignment(Alignment &next) {
		int qOffset = next.qPos - qPos;
		int tOffset = next.tPos - tPos;
		Block tempBlock(0,0,0);
		int n;
		for (n = 0; n < next.blocks.size(); n++ ) {
			tempBlock = next.blocks[n];
			tempBlock.qPos += qOffset;
			tempBlock.tPos += tOffset;
			blocks.push_back(tempBlock);
		}
	}

	void ArrowPathToAlignment(vector<Arrow> &optPath) {
		int q, t;
		int a = 1;
		q = 0; t = 0;
		Block b;
		a = 0;
		while (a < optPath.size()) {
			if (optPath[a] == Diagonal) {
				// Start of a block;
				b.qPos = q;
				b.tPos = t;
				b.length = 0;
				while(a < optPath.size() and optPath[a] == Diagonal) {
					b.length++;
					a++;
					t++;
					q++;
				}
				blocks.push_back(b);
			}
			while(a < optPath.size() and optPath[a] == Left) {
				t++;
				a++;
			}
			while(a < optPath.size() and optPath[a] == Up) {
				q++;
				a++;
			}
		}
	}
	

	void CreateAlignmentStrings(char *query, char* text,
															string &queryStr, string &alignStr, string &textStr) {
		GenomePos q = qPos;
		GenomePos t = tPos;
		GenomePos qPos, tPos;
		GenomePos  g;
		char mismatchChar = '*';
		char matchChar = '|';
		char gapChar = '-';
		char gapSeparationChar = ' ';
		refLen=0;
		int refStart=0;
		if (blocks.size() == 0) {
			textStr = "";
			alignStr = "";
			queryStr = "";
			return;
		}

		if (blocks.size() > 0) {
			q=blocks[0].qPos;
			t=blocks[0].tPos;
		}

		for (int b = 0; b < blocks.size() ; b++) {

			for (int bl = 0; bl < blocks[b].length; bl++ ) {
				assert(t < genomeLen);
				queryStr.push_back(query[q]);
				textStr.push_back(text[t]);
				if (seqMap[query[q]] != 
						seqMap[text[t]])
					alignStr.push_back(mismatchChar);
				else
					alignStr.push_back(matchChar);
				q++;
				t++;
			}
			//
			//  There are no gaps to count after the last block, so 
			//  don't add the gapped characters for this.
			//
			if (blocks.size() == 0)
				continue;
			if (b == blocks.size() - 1) {
				continue;
			}


			int queryGapLen = (blocks[b+1].qPos - 
												 blocks[b].qPos - blocks[b].length);
			int textGapLen  = (blocks[b+1].tPos - 
												 blocks[b].tPos - blocks[b].length);
			assert(queryGapLen >= 0);
			assert(textGapLen >= 0);
			if (queryGapLen > 0 or textGapLen > 0) {
				// commonGapLen should be the shorter gap.
				int commonGapLen = queryGapLen; 
				if (queryGapLen > textGapLen) {
					commonGapLen = textGapLen;
				}
				textGapLen -= commonGapLen;
				queryGapLen -= commonGapLen;

				for (g = 0; g < queryGapLen; g++, q++ ){
					assert(t < genomeLen);
					textStr.push_back(gapChar);
					alignStr.push_back(gapSeparationChar);
					queryStr.push_back(query[q]);
				}
				for (g = 0; g < textGapLen; g++, t++ ){
					assert(t < genomeLen);
					textStr.push_back(text[t]);
					alignStr.push_back(gapSeparationChar);
					queryStr.push_back(gapChar);

				}

				for (g = 0; g < commonGapLen; g++ ) {
					assert(t < genomeLen);
					textStr.push_back(text[t]);
					if (seqMap[query[q]] != 
							seqMap[text[t]])
						alignStr.push_back(mismatchChar);
					else
						alignStr.push_back(matchChar);

					queryStr.push_back(query[q]);
					t++;
					q++;
				}
			}
		}
		refLen=t-refStart;
	}

	void AlignStringsToCigar(string &query, string &target, string &cigar, int &nm, int &nmm, int &nins, int &ndel) {
		stringstream cigarstrm;
		int i=0;
		int p=0;
		nm=nmm=nins=ndel=0;
		while (i < query.size()) {
			p=i;
			while (i < query.size() and seqMap[query[i]] == seqMap[target[i]] and query[i] != '-' and target[i] != '-') {	i++;}
			
			if (i > p) {
				cigarstrm << i-p << '=';
				nm+=i-p;
				continue;
			}
			while (i < query.size() and seqMap[query[i]] != seqMap[target[i]] and query[i] != '-' and target[i] != '-') {	i++;}
			if (i > p) {
				cigarstrm << i-p << 'X';
				nmm+=i-p;
				continue;
			}
			while (i < query.size() and query[i] == '-' and target[i] != '-') {	i++;}
			if (i > p) {
				cigarstrm << i-p << 'D';
				ndel+=i-p;
				continue;
			}
			while (i < query.size() and query[i] != '-' and target[i] == '-') {	i++;}
			if (i > p) {
				cigarstrm << i-p << 'I';
				nins+=i-p;
				continue;
			}
		}
		cigar=cigarstrm.str();
	}
	void CalculateStatistics(int size, int cur) {

		CreateAlignmentStrings(read, genome, queryString, alignString, refString);
		AlignStringsToCigar(queryString, refString, cigar, nm, nmm, ndel, nins);
		preClip = 0;
		sufClip=0;
		if (blocks.size() > 0) {
			int last=blocks.size();
			preClip = blocks[0].qPos;
			sufClip = readLen - blocks[last-1].qPos - blocks[last-1].length;
			qStart = blocks[0].qPos;
			qEnd   = blocks[last-1].qPos + blocks[last-1].length;
			tStart = blocks[0].tPos;
			tEnd   = blocks[last-1].tPos + blocks[last-1].length;
		}
		//
		// Flag that the stats are calculated for methods that need them.
		// 
		prepared=true;
		if (size == 1) {
			if (strand == 1) {
				flag = flag | READ_REVERSE;
			}
		}
		else if (cur == 0) { // this is a chimeric alignment
			if (strand == 1) {
				flag = flag | READ_MULTIPLESEGMENTS | READ_FIRSTSEGMENT | READ_SUPPLEMENTARY | READ_REVERSE;
			}
			else {
				flag = flag | READ_MULTIPLESEGMENTS | READ_FIRSTSEGMENT | READ_SUPPLEMENTARY;
			}
		}
		else if (cur == size - 1) {
			if (strand == 1) {
				flag = flag | READ_MULTIPLESEGMENTS | READ_LASTSEGMENT | READ_SUPPLEMENTARY | READ_REVERSE;
			}
			else {
				flag = flag | READ_MULTIPLESEGMENTS | READ_LASTSEGMENT | READ_SUPPLEMENTARY;
			}
		}
		else {
			if (strand == 1) {
				flag = flag | READ_MULTIPLESEGMENTS | READ_LASTSEGMENT | READ_FIRSTSEGMENT | READ_SUPPLEMENTARY | READ_REVERSE;
			}
			else {
				flag = flag | READ_MULTIPLESEGMENTS | READ_LASTSEGMENT | READ_FIRSTSEGMENT | READ_SUPPLEMENTARY;
			}
		}
	}
	
	bool Overlaps(const Alignment &b, float frac) const {
		int ovp=0;
		if (b.qStart >= qStart and b.qStart < qEnd) {
			ovp=min(qEnd, b.qEnd)-b.qStart;
		}
		else if (b.qEnd > qStart and b.qEnd < qEnd) {
			ovp=b.qEnd-max(qStart, b.qStart);
		}
		else if (b.qStart <= qStart and b.qEnd > qEnd) {
			ovp=qEnd-qStart;
		}
		float denom=qEnd-qStart;
		if (ovp/denom > frac) { return true; }
		else { return false; }
	}

	void PrintPairwise(ostream &out) const {
		assert(prepared);				
				
		int i=0;
		int q=0;
		int t=0;
		int nBlocks = blocks.size();
		out << readName << endl;
		if (nBlocks > 0) {
			out << "Interval:\t" << chrom<< ":" << blocks[0].tPos << "-" << blocks[0].tPos +refLen << endl;
		}
		while (i < queryString.size()) {
			int end = min((int) queryString.size(), i+50);
			string qsub = queryString.substr(i,end-i);
			out.width(10);
			out << q + GetQStart() << " q: " << qsub << endl;
			q+= qsub.size() - std::count(qsub.begin(),qsub.end(),'-');
			out << "              " << alignString.substr(i,end-i) << endl;
			string tsub = refString.substr(i,end-i);
			out.width(10);
			out << t + GetTStart() << " t: " << tsub << endl;
			t+= tsub.size() - std::count(tsub.begin(), tsub.end(),'-');
			cout <<endl;
			i=end;
		}
	}

	void PrintBed(ostream &out) {
		out  << chrom << "\t" 
				 << tStart << "\t" 
				 << tEnd << "\t"
				 << (int) mapqv << "\t" 
				 << readName << "\t" << readLen << "\t" << qStart << "\t" << qEnd << "\t"
				 << nm << "\t" << nblocks << endl;
	}

	void PrintSAM(ostream &out, Options &opts, char *passthrough=NULL) {
		stringstream samStrm;
		samStrm << readName << "\t";
		assert(prepared);
		if (blocks.size() == 0) {
			//
			// Create a null alignment
			//
			cerr << "will create this later." << endl;
		}
		else {

			int last = blocks.size();
			samStrm << (unsigned int) flag << "\t"
							<< chrom << "\t" 
							<< tStart+1 << "\t"
							<< (unsigned int) mapqv << "\t";
			char clipOp = 'S';
			if (opts.hardClip) {
				clipOp = 'H';
			}
			if (preClip > 0) {
				samStrm << preClip << clipOp;
			}
			samStrm << cigar;
			if (sufClip > 0) {
				samStrm << sufClip << clipOp;
			}
			// Rnext, Pnext
			samStrm << "\t*\t0\t";
			// Template length
			samStrm << tEnd - tStart << "\t";
			if (opts.hardClip) {
				string subStr;
				subStr=string(read, blocks[0].qPos, blocks[last-1].qPos + blocks[last-1].length);
				samStrm << subStr;
			}
			else {
				string readStr(read, 0, readLen);
				samStrm << readStr;
			}
			samStrm << "\t";
			samStrm << "*";
		}
		out << samStrm.str();
		if (passthrough != NULL ) {
			out << "\t" << passthrough;
		}
		out << endl;
		out.flush();
	}
};



class SegAlignmentGroup {
public:
	vector<Alignment*> SegAlignment;
	GenomePos qStart, qEnd, tStart, tEnd;
	unsigned char mapqv;
	int nm;
	SegAlignmentGroup () {};
	~SegAlignmentGroup () {};

	void SetBoundariesFromSegAlignmentAndnm (Read & read) {
		qStart = 0, qEnd = 0, tStart = 0, tEnd = 0, nm = 0;
		for (int s = 0; s < SegAlignment.size(); s++) {
			if (SegAlignment[s]->strand == 0) {
				qStart = min(qStart, SegAlignment[s]->qStart);
				qEnd   = max(qEnd, SegAlignment[s]->qEnd);
				tStart = min(tStart, SegAlignment[s]->tStart);
				tEnd   = max(tEnd, SegAlignment[s]->tEnd);
			}
			else {
				qStart = min(qStart, read.length - SegAlignment[s]->qEnd);
				qEnd   = max(qEnd, read.length - SegAlignment[s]->qStart);
				tStart = min(tStart, SegAlignment[s]->tStart);
				tEnd   = max(tEnd, SegAlignment[s]->tEnd);				
			}
			nm += SegAlignment[s]->nm;
		}

	}

	bool Overlaps(const SegAlignmentGroup &b, float frac) const {
		int ovp=0;
		if (b.qStart >= qStart and b.qStart < qEnd) {
			ovp=min(qEnd, b.qEnd)-b.qStart;
		}
		else if (b.qEnd > qStart and b.qEnd < qEnd) {
			ovp=b.qEnd-max(qStart, b.qStart);
		}
		else if (b.qStart <= qStart and b.qEnd > qEnd) {
			ovp=qEnd-qStart;
		}
		float denom=qEnd-qStart;
		if (ovp/denom > frac) { return true; }
		else { return false; }
	}

	void SetMapqv () {
		for (int s = 0; s < SegAlignment.size(); s++) {
			SegAlignment[s]->mapqv = mapqv;
		}
	}

};


#endif
