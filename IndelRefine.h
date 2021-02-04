#ifndef INDEL_REFINE_H_
#define INDEL_REFINE_H_

#include "Alignment.h"
#include "Genome.h"
#include "Read.h"
#include "Options.h"
#include "AffineOneGapAlign.h"
#include <algorithm>

void PrintMat(vector<int> &mat, vector<int> &qS, vector<int> &qE) {
	int m=0;
	int s=qS[0];
	int w=4;
	cout << setw(w) << ".";
	for (int i=0; i < qE[qE.size()-1]; i++) {
		cout << setw(w) << i;
	}
	cout << endl;
	for (int i=0; i < qS.size(); i++) {
		cout << setw(w) << i;
		for (int j=0; j < qS[i] - s; j++) { cout << setw(w) << "."; }
		for (int j=0; j < qE[i] - qS[i] + 1; j++) { cout << setw(w) << mat[m]; m++;} cout << endl;
	}
}

void IndelRefineAlignment(Read &read, 
													Genome &genome, 
													Alignment &alignment, 
													Options &opts) {

	int startBlock=0, endBlock=0;
	int k=10;
	int maxGap=k-1;
	long qPos,tPos;

	vector<Block> refined;

	//
	// No modification on empty or ungapped alignment
	//
	if (alignment.blocks.size() == 0 or alignment.blocks.size() == 1) { return; }
	int st=0,sq=0;
	int tbe=min(10,(int) alignment.blocks.size());
	/*
	for (int tb=0;tb<tbe; tb++) {
		
		cout << "tb\t" << alignment.blocks[tb].qPos - alignment.blocks[0].qPos << "\t"
				 << alignment.blocks[tb].tPos - alignment.blocks[0].tPos << "\t"
				 << alignment.blocks[tb].length << "\t" << sq << "\t" << sq << endl;
	}
	*/
	//	cout << "ql: " << alignment.blocks[tbe-1].qPos + alignment.blocks[tbe-1].length - alignment.blocks[0].qPos << endl;
	//	cout << "tl: " << alignment.blocks[tbe-1].tPos + alignment.blocks[tbe-1].length - alignment.blocks[0].tPos << endl;

	while (endBlock < alignment.blocks.size()) {
		
		//
		// Indel is measured from next qpos/tpos to current position
		//
		long qStart = alignment.blocks[startBlock].qPos;
		long tStart = alignment.blocks[startBlock].tPos;
		int blockLen=alignment.blocks[startBlock].length;
		qPos=alignment.blocks[startBlock].qPos + blockLen;
		tPos=alignment.blocks[startBlock].tPos + blockLen;
		
		while (endBlock < alignment.blocks.size() - 1 and 
					 alignment.blocks[endBlock+1].qPos - qPos < maxGap and
					 alignment.blocks[endBlock+1].tPos - tPos < maxGap) {
			
			endBlock++;
			int blockLen=alignment.blocks[endBlock].length;
			qPos=alignment.blocks[endBlock].qPos + blockLen;
			tPos=alignment.blocks[endBlock].tPos + blockLen;
		}
		
		//
		// Just one block, do not refine.
		//
		if (endBlock == startBlock) {
			refined.push_back(alignment.blocks[startBlock]);
		}
		else {
			//
			// Multiple blocks, need to do a banded alignment
			//		
			long qEnd   = alignment.blocks[endBlock].qPos + alignment.blocks[endBlock].length;
			long tEnd   = alignment.blocks[endBlock].tPos + alignment.blocks[endBlock].length;

			long qLen=qPos-qStart;
			long tLen=tPos-tStart;

			vector<int> qS(tLen, -1), qE(tLen, -1), tS(qLen, -1), tE(qLen, -1);
			long t, q;
			int tOff=0, qOff=0;
			t=alignment.blocks[startBlock].tPos;
			q=alignment.blocks[startBlock].qPos;
			for (int b=startBlock; b <= endBlock; b++) {
				int qGap=0, tGap=0;
				int commonGap   = 0;
				int blockLength = alignment.blocks[b].length;
				if (b < endBlock) {
					qGap = alignment.blocks[b+1].qPos - 
						(alignment.blocks[b].qPos + blockLength);
					tGap = alignment.blocks[b+1].tPos - 
						(alignment.blocks[b].tPos + blockLength);

					if (qGap > 0 and tGap > 0) {
							commonGap=min(qGap,tGap);
							qGap-=commonGap;
							tGap-=commonGap;
							blockLength+=commonGap;
					}
				}
				//
				// Process contiguous alignment.
				//
				for (int bi = 0; bi < blockLength; tOff++, bi++, q++, t++) {
					if (qS[tOff] == -1) {
						qS[tOff] = max(q-k, qStart);
					}
					else {
						qS[tOff] = min((long)qS[tOff], max(q-k, (long)qStart));
						assert(qS[tOff] >= 0);
					}
					if (qE[tOff] == -1 or qE[tOff] < q+k) {
						qE[tOff] = min(qEnd-1, (long)(q+k));
						assert(qE[tOff] < qEnd);
					}
					for (int ki=0; ki < k; ki++) {
						//
						// ensure future entries process to this point k above and below.
						if (tOff - ki >= 0) {
							if (qE[tOff - ki]  < q) { 
								qE[tOff-ki] = q;
							}
						}
						if (tOff + ki < qS.size()) {
							if (qS[tOff+ki] == -1 or qS[tOff+ki] > q) {
								qS[tOff+ki] = q;
							}
						}
					}
				}
				//
				// Advance gap
				//
				if (qGap > tGap) {
					assert(tGap == 0);
					for (int qi=0; qi < qGap; qi++, q++) {
						for (int ki=0; ki < k; ki++) {
							if (tOff - ki >= 0) {
								if (qE[tOff - ki]  < q) { qE[tOff-ki] = q;}
							}
							if (tOff + ki < qS.size()) {
								if (qS[tOff+ki] == 0 or qS[tOff+ki] > q) {
									qS[tOff+ki] = q;
								}
							}
						}
					}
				}
				if (tGap > qGap) {
					assert(qGap == 0);
					for (int ti=0; ti < tGap; tOff++, ti++, t++) {
						qS[tOff] = max(q-k, qStart);
						qE[tOff] = min(qEnd-1,q+k);
					}
				}
			}
			
			long matSize=0;
			for (int qi=qS.size(); qi >1; qi--) {
				if (qS[qi-1] < qS[qi-2]) {
					qS[qi-2] = qS[qi-1];
				}
			}
			for (int qi=0; qi < qS.size()-1; qi++) {
				if (qE[qi] > qE[qi+1]) { qE[qi+1] = qE[qi]; }
			}
			for (int qi=0; qi < qS.size(); qi++) {
				matSize+=qE[qi]-qS[qi]+1;
			}

			vector<int> scoreMat(matSize,0);
			vector<int> pathMat(matSize,0);
			vector<int> indexMat(matSize,-1);
			vector<int> cuMatSize(qS.size(), 0);

			char *tSeq=genome.seqs[alignment.chromIndex];
			char *qSeq=alignment.read;
			long tSeqLen = tEnd-tStart;
			long qSeqLen = qEnd-qStart;
			int gap=opts.localIndel;
			int match=opts.localMatch;
			int mismatch=opts.localMismatch;
			
			if (tSeqLen < k or qSeqLen < k) {

				Alignment aln;
				AffineAlignBuffers buff;
				string qStr(&qSeq[qStart], qSeqLen);
				string tStr(&tSeq[tStart], tSeqLen);
				AffineOneGapAlign(qStr, qSeqLen, tStr, tSeqLen, 
													match, mismatch, gap, k, aln, buff);
				for (int afb=0; afb < aln.blocks.size(); afb++) {
					aln.blocks[afb].qPos += qStart;
					aln.blocks[afb].tPos += tStart;
				}
				refined.insert(refined.end(), aln.blocks.begin(), aln.blocks.end());
			}
			else {
				 for (int qi=1; qi < qS.size(); qi++) {
					 cuMatSize[qi] = cuMatSize[qi-1] + qE[qi-1] - qS[qi-1] + 1;
			 }
			 //
			 // The first base is always aligned here. 
			 int rowStart=0, rowEnd=-1;
			 int BAD=-99999999;
			 int diag=0;
			 int left=1;
			 int down=2;
			 int bound=3;
			 for (int ti=0; ti < tLen; ti++) {
				 int rowLen=qE[ti]-qS[ti]+1;
				 rowEnd=rowStart + rowLen-1;
				 if (rowStart > 0) {
					 scoreMat[rowStart] = BAD;
					 pathMat[rowStart] = bound;
				 }
				 else {
					 for (int qi = 1; qi < rowEnd; qi++ ){
						 scoreMat[qi] = scoreMat[qi-1] + gap;
						 pathMat[qi]  = left;
						 indexMat[qi] = qi-1;
					 }
				 }
				 if (ti < tLen-1) {
					 scoreMat[rowEnd] = BAD;
					 pathMat[rowEnd]  = bound;
				 }
				 rowStart+=rowLen;
			 }
			 //
			 // Now run dp
			 //
			 int curRowStart=qE[0]-qS[0]+1;
			 int prevRowStart=0;
			 int prevRowLen=qE[0]-qS[0]+1;
			 for (int ti=1; ti < tLen; ti++) {
				 int curRowLen=qE[ti] - qS[ti] + 1;
				 assert(qS[ti-1] <= qS[ti]);
				 int curRowOffset=qS[ti] - qS[ti-1];
				 // Iterate on current row, skipping boundaries.
				 // prevRowPos is set to the cell immediately below the current cell, so the diagonal cell is prevRowPos-1
				 // The +1 offset is to skip past boundary cells.
				 int prevRowPos=prevRowStart + curRowOffset +1;
				 int curRowPos=curRowStart + 1;
				 //				cout << "Row offsets\t" << ti << "\t" << prevRowPos << "\t" << curRowPos << endl;
				 int rowEnd;
				 // Last row is solved to final cell which should end on a match.
				 if (ti == tLen-1) {
					 rowEnd=curRowLen;
				 }
				 else {
					 rowEnd=curRowLen-1;
				 }
				 for (int qi=1; qi < rowEnd; qi++, curRowPos++, prevRowPos++) {
					 int matchScore, insScore, delScore;
					 int matchIndex, insIndex, delIndex;
					 if (qE[ti-1] >= qi + qS[ti]) {
						 if (tSeq[ti+tStart] == qSeq[qi+qS[ti]]) {
							 matchScore = scoreMat[prevRowPos-1] + match;
						 }
						 else { 
							 matchScore = scoreMat[prevRowPos-1] + mismatch;
						 }
					 }
					 else {
						 matchScore = BAD;
					 }

					 insScore = scoreMat[curRowPos-1] + gap;

					 if (qE[ti-1] >= qi+qS[ti]) {
						 delScore = scoreMat[prevRowPos] + gap;
					 }
					 else {
						 delScore =BAD;
					 }
					 int maxScore = max(matchScore,max(insScore, delScore));
					 scoreMat[curRowPos] = maxScore;
					 if (maxScore == matchScore) {
						 pathMat[curRowPos]  = diag;
						 assert(prevRowPos-1 <= cuMatSize[ti]);
						 indexMat[curRowPos] = prevRowPos-1;
					 }
					 else if (maxScore == insScore) {
						 pathMat[curRowPos]  = left;
						 indexMat[curRowPos] = curRowPos-1;
					 }
					 else {
						 pathMat[curRowPos]  = down;
						 assert(prevRowPos <= cuMatSize[ti]);
						 indexMat[curRowPos] = prevRowPos;
					 }
				 }

				 prevRowStart += prevRowLen;
				 curRowStart  += curRowLen;
				 prevRowLen   =  curRowLen;
			 }
			 //
			 // Trace back.
			 // 
			 vector<int> path;
			 int curMatPos=pathMat.size()-1;

			 while (curMatPos > 0) {
				 path.push_back(pathMat[curMatPos]);
				 pathMat[curMatPos] = 90+pathMat[curMatPos];
				 curMatPos=indexMat[curMatPos];				
			 }
			 path.push_back(diag);

			 int ci=0;
			 int totalI=0;
			 //			if (startBlock == 0) {

							 //			}

			 for (int r=0; r < qS.size(); r++) {
				 int np=0;
				 for (int ri=0; ri < qE[r] - qS[r] + 1; ri++) {
					 if (pathMat[ci] == 90 || pathMat[ci] == 92) { np++; totalI++;}
					 ++ci;

				 }
				 if (r > 0 and np != 1) {
					 PrintMat(pathMat, qS, qE);
					 cout << "ERROR at " << r << "\t" << np << endl;

				 }
				 assert(r == 0 || np == 1);
			 }


			 //			PrintMat(pathMat, qS, qE);
			 //			PrintMat(scoreMat, qS, qE);
			 long qPath=qStart;
			 long tPath=tStart;
			 int pi=path.size()-1;
			 reverse(path.begin(), path.end());
			 assert(path[0] == diag);
			 //			assert(path[pi] == diag);

			 pi=0;
			 while (pi < path.size()) {
				 int blockLen=0;
				 while (pi < path.size() and path[pi] == diag) {
					 blockLen++;
					 pi++;
				 }
				 int tgapLen=0;
				 int qgapLen=0;

				 if (path[pi] == left) {
					 while (pi < path.size() and path[pi] == left) {
						 qgapLen++;
						 pi++;
					 }
				 }
				 else if (path[pi] == down) {
					 while (pi < path.size() and path[pi] == down) {
						 tgapLen++;
						 pi++;
					 }
				 }

				 refined.push_back(Block(qPath,tPath, blockLen));
				 qPath+=blockLen+qgapLen;
				 tPath+=blockLen+tgapLen;
			 }				
			
			 if (endBlock < alignment.blocks.size()-1) {
				 assert(tPath <= alignment.blocks[endBlock+1].tPos);
				 assert(qPath <= alignment.blocks[endBlock+1].qPos);
			 }
			}

			//cerr << "ratio " << scoreMat.size() / (float)qS.size() << endl;			
		}
		endBlock++;
		startBlock=endBlock;		
	}
	alignment.blocks = refined;
}


#endif
