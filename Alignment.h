#ifndef ALIGNMENT_TYPE_H_
#define ALIGNMENT_TYPE_H_

#include <vector>
using namespace std;

#include "AlignmentBlock.h"
#include "Path.h"
#include "SeqUtils.h"
#include "Options.h"
#include <sstream>
const unsigned int READ_UNMAPPED=0x4;
const unsigned int READ_REVERSE=0x10;
const unsigned int READ_SECONDARY=0x100;

class Alignment {
 public:
	unsigned char mapqv;
	unsigned int flag;
	string chrom;
	string name;
	string queryString, alignString, refString;
	int nm, nmm, nins, ndel;
	int preClip, sufClip;
	string cigar;
	bool prepared;
	char *read;
	char *forward;
	int readLen;
	string readName;
	char *genome;
	int strand;

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
		preClip=0; sufClip=0;
		prepared=false;
		read=NULL;
		forward=NULL;
	}
 Alignment(char *_read, char *_forward, 
					 int _rl, string _rn, int _str, 
					 char *_genome, string &_chrom) : Alignment() { 
		read=_read; 
		forward=_forward;
		readLen = _rl; 
		readName= _rn;
		strand= _str;
		genome=_genome;
		chrom=_chrom;
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

			if (queryGapLen > 0 or textGapLen > 0) {
				// commonGapLen should be the shorter gap.
				int commonGapLen = queryGapLen; 
				if (queryGapLen > textGapLen) {
					commonGapLen = textGapLen;
				}
				textGapLen -= commonGapLen;
				queryGapLen -= commonGapLen;

				for (g = 0; g < queryGapLen; g++, q++ ){
					textStr.push_back(gapChar);
					alignStr.push_back(gapSeparationChar);
					queryStr.push_back(query[q]);
				}
				for (g = 0; g < textGapLen; g++, t++ ){
					textStr.push_back(text[t]);
					alignStr.push_back(gapSeparationChar);
					queryStr.push_back(gapChar);

				}

				for (g = 0; g < commonGapLen; g++ ) {
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
	}

	void AlignStringsToCigar(string &query, string &target, string &cigar, int &nm, int &nmm, int &nins, int &ndel) {
		stringstream cigarstrm;
		int i=0;
		int p=0;
		nm=nmm=nins=ndel=0;
		while (i < query.size()) {
			p=i;
			while (seqMap[query[i]] == seqMap[target[i]] and query[i] != '-' and target[i] != '-') {	i++;}
			if (i > p) {
				cigarstrm << i-p << '=';
				nm+=i-p;
				continue;
			}
			while (seqMap[query[i]] != seqMap[target[i]] and query[i] != '-' and target[i] != '-') {	i++;}
			if (i > p) {
				cigarstrm << i-p << 'X';
				nmm+=i-p;
				continue;
			}
			while (query[i] == '-' and target[i] != '-') {	i++;}
			if (i > p) {
				cigarstrm << i-p << 'D';
				ndel+=i-p;
				continue;
			}
			while (query[i] != '-' and target[i] == '-') {	i++;}
			if (i > p) {
				cigarstrm << i-p << 'I';
				nins+=i-p;
				continue;
			}
		}
		cigar=cigarstrm.str();
	}
	void CalculateStatistics() {

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
		if (strand == 1) {
			SetReverse();
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
		
		while (i < queryString.size()) {
			int end = min((int) queryString.size(), i+50);
			string qsub = queryString.substr(i,end-i);
			out.width(10);
			out << q + GetQStart() << " q: " << qsub << endl;
			q+= qsub.size() - count(qsub.begin(),qsub.end(),'-');
			out << "              " << alignString.substr(i,end-i) << endl;
			string tsub = refString.substr(i,end-i);
			out.width(10);
			out << t + GetTStart() << " t: " << tsub << endl;
			t+= tsub.size() - count(tsub.begin(), tsub.end(),'-');
			cout <<endl;
			i=end;
		}
	}

	void PrintBed(ostream &out) {
		out  << chrom << "\t" 
				 << tStart << "\t" 
				 << tEnd << "\t"
				 << (int) mapqv << "\t" 
				 << readName << "\t" << readLen << "\t"
				 << nm << "\t" << endl;
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
				if (strand == 0) {
					subStr=string(forward, blocks[0].qPos, blocks[last-1].qPos + blocks[last-1].length);
				}
				else {
					int start = readLen - (blocks[last-1].qPos + blocks[last-1].length);
					int end   = readLen - blocks[0].qPos;
					subStr=string(forward, start, end);
				}

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
	}
};


#endif
