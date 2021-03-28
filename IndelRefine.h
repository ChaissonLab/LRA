#ifndef INDEL_REFINE_H_
#define INDEL_REFINE_H_

#include "Alignment.h"
#include "Genome.h"
#include "Read.h"
#include "Options.h"
#include "AffineOneGapAlign.h"
#include <algorithm>

void FlatPrintMat(vector                        <int> &mat, vector<int> &qS, vector<int> &qE) {
	int mi=0;
	for (int r=0; r< qS.size(); r++) {
		cout << r << "\t" << qS[r] << "\t" << qE[r] << "\t";
		for (int c=qS[r]; c <= qE[r]; c++, mi++) {
			cout << setw(3) << mat[mi];
		}
		cout << endl;
	}
}

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
class IndelRefineBuffers {
 public:

	vector<int> cuMatSize;
	vector<int> qS, qE;

	vector<int> insScoreMat, insPathMat, insIndexMat;
	vector<int> delScoreMat, delPathMat, delIndexMat;
	
	vector<int> scoreMat;
	vector<int> pathMat;
	vector<int> indexMat;
	

};

void IndelRefineAlignment(Read &read, 
													Genome &genome, 
													Alignment &alignment, 
													Options &opts,
													IndelRefineBuffers &buffers) {

	int startBlock=0, endBlock=0;
	int k=opts.refineBand;
	int maxGap=k-1;
	long qPos,tPos;

	vector<Block> refined;

	vector<int> &cuMatSize=buffers.cuMatSize;
	vector<int> &qS=buffers.qS;
	vector<int> &qE=buffers.qE;

	vector<int> &insScoreMat=buffers.insScoreMat, &insPathMat=buffers.insPathMat, &insIndexMat=buffers.insIndexMat;
	vector<int> &delScoreMat=buffers.delScoreMat, & delPathMat=buffers.delPathMat, & delIndexMat=buffers.delIndexMat;
	
	vector<int> &scoreMat=buffers.scoreMat;
	vector<int> &pathMat=buffers.pathMat;
	vector<int> &indexMat=buffers.indexMat;
	//
	// No modification on empty or ungapped alignment
	//
	if (alignment.blocks.size() == 0 or alignment.blocks.size() == 1) { return; }
	int st=0,sq=0;
	int tbe=min(10,(int) alignment.blocks.size());

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


			qS.resize(tLen, -1);
			qE.resize(tLen, -1);
			fill(qS.begin(), qS.end(),-1);
			fill(qE.begin(), qE.end(), -1);
			long t, q;
			int tOff=0, qOff=0;
			t=alignment.blocks[startBlock].tPos;
			q=alignment.blocks[startBlock].qPos;
			//			cout << "Block boundaries " << endl;
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
						//						cout << "START " << startBlock << "\tROW " << tOff << "\tqS " << qS[tOff];
					}
					else {
						qS[tOff] = min((long)qS[tOff], max(q-k, (long)qStart));
						//						cout << "START " << startBlock << "\tROW " << tOff << "\tqS " << qS[tOff] << endl;
						assert(qS[tOff] >= 0);
					}
					if (qE[tOff] == -1 or qE[tOff] < q+k) {
						qE[tOff] = min(qEnd-1, (long)(q+k));
						//						cout << "START " << startBlock << "\tROW " << tOff << "\tqE " << qE[tOff] << endl;
						assert(qE[tOff] < qEnd);
					}
					for (int ki=0; ki < k; ki++) {
						//
						// ensure future entries process to this point k above and below.
						if (tOff - ki >= 0) {
							if (qE[tOff - ki]  < q) { 
								//								cout << "START " << startBlock << "\tROW " << tOff-ki << "\tResetting qE " << tOff-ki << "\t" << qE[tOff-ki] << "\t" << q << endl;
								qE[tOff-ki] = q;
							}
						}
						if (tOff + ki < qS.size()) {
							if (qS[tOff+ki] == -1 or qS[tOff+ki] > q) {
								//								cout << "START " << startBlock << "\tROW " << tOff+ki << "\tResetting qS " << tOff+ki << "\t" << qS[tOff+ki] << "\t" << q << endl;
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
								if (qE[tOff - ki]  < q) { 
									//									cout << "START " << startBlock << "\tROW " << tOff-ki << "\tgap qe " << tOff-ki << "\t" << qE[tOff-ki] << "\t" << q << endl;
									qE[tOff-ki] = q;									
								}
							}
							if (tOff + ki < qS.size()) {
								if (qS[tOff+ki] == 0 or qS[tOff+ki] > q) {
									//									cout << "START " << startBlock << "\tROW " << tOff+ki << "\tgap qs " << tOff+ki << "\t" << qE[tOff+ki] << "\t" << q << endl;
									qS[tOff+ki] = q;
								}
							}
						}
					}
				}
				if (tGap > qGap) {
					assert(qGap == 0);
					for (int ti=0; ti < tGap; tOff++, ti++, t++) {
						//						cout << "START " << startBlock << "\tROW " << tOff << " qs gap " << qS[tOff] << "\t" << q-k << "\t" << qStart << endl;
						qS[tOff] = max(q-k, qStart);
						//						cout << "START " << startBlock << "\tROW " << tOff << " qe gap " << qE[tOff] << "\t" << q+k << "\t" << qEnd -1 << endl;
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
			for (int qi=1; qi < qS.size(); qi++) {
				assert(qS[qi] >= qS[qi-1]);
				assert(qE[qi] >= qE[qi-1]);
				assert(qS[qi] < qE[qi]);
			}
			char *tSeq=genome.seqs[alignment.chromIndex];
			char *qSeq=alignment.read;
			long tSeqLen = tEnd-tStart;
			long qSeqLen = qEnd-qStart;
			int gap=opts.localIndel;
			int gapOpen=opts.localIndel*2+1;
			int gapExtend=0;//opts.localMismatch+1;
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
				cuMatSize.resize(qS.size(), 0);
				fill(cuMatSize.begin(), cuMatSize.end(), 0);
				for (int qi=1; qi < qS.size(); qi++) {
					cuMatSize[qi] = cuMatSize[qi-1] + qE[qi-1] - qS[qi-1] + 1;
				}

				//
			 // The first base is always aligned here. 
			 int rowStart=0, rowEnd=-1;
			 int BAD=-999999999;
			 int diag=0;
			 int left=1;
			 int down=2;
			 int bound=3;
			 int done=20;
			 int delOpen=4;
			 int delExtend=5;
			 int delClose=6;

			 int insOpen=7;
			 int insExtend=8;
			 int insClose=9;


			 scoreMat.resize(matSize, 0);
			 fill(scoreMat.begin(), scoreMat.end(), 0);

			 pathMat.resize(matSize, bound);
			 fill(pathMat.begin(), pathMat.end(), bound);

			 indexMat.resize(matSize, -1);
			 fill(indexMat.begin(), indexMat.end(), -1);

			 delScoreMat.resize(matSize, BAD);
			 fill(delScoreMat.begin(), delScoreMat.end(), BAD);
			 delPathMat.resize(matSize, 0);
			 fill(delPathMat.begin(), delPathMat.end(), bound);
			 delIndexMat.resize(matSize, -1);
			 fill(delIndexMat.begin(), delIndexMat.end(), -1);


			 insScoreMat.resize(matSize, BAD);
			 fill(insScoreMat.begin(), insScoreMat.end(), BAD);
			 insPathMat.resize(matSize, 0);
			 fill(insPathMat.begin(), insPathMat.end(), bound);
			 insIndexMat.resize(matSize, -1);
			 fill(insIndexMat.begin(), insIndexMat.end(), -1);
			 
			 indexMat[0] = 0;
			 pathMat[0] = done;
			 for (int ti=0; ti < tLen; ti++) {
				 int rowLen=qE[ti]-qS[ti]+1;
				 rowEnd=rowStart + rowLen-1;
				 if (rowStart > 0) {
					 scoreMat[rowStart] = BAD;
					 pathMat[rowStart] = bound;
					 insPathMat[rowStart] = bound;

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
			 int curRowStart  = qE[0]-qS[0]+1;
			 int prevRowStart = 0;
			 int prevRowLen   = qE[0]-qS[0]+1;
			 for (int ti=1; ti < tLen; ti++) {
				 int curRowLen = qE[ti] - qS[ti] + 1;
				 assert(qS[ti-1] <= qS[ti]);
				 int curRowOffset = qS[ti] - qS[ti-1];
				 // Iterate on current row, skipping boundaries.
				 // prevRowPos is set to the cell immediately below the current cell, so the diagonal cell is prevRowPos-1
				 // The +1 offset is to skip past boundary cells.
				 int curRowPos  = curRowStart + 1;
				 int prevRowPos = prevRowStart + curRowOffset +1;
				 int rowEnd;
				 // Last row is solved to final cell which should end on a match.
				 if (ti == tLen-1) {
					 rowEnd = curRowLen;
				 }
				 else {
					 rowEnd = curRowLen - 1;
				 }
				 int *scoreMatPrevRowPosPtr=&scoreMat[prevRowPos];
				 int *scoreMatCurRowPosPtr=&scoreMat[curRowPos];
				 int *scoreMatDiagPosPtr=&scoreMat[prevRowPos-1];
				 int *scoreMatInsPosPtr=&scoreMat[curRowPos-1];

				 int *pathMatPrevRowPosPtr=&pathMat[prevRowPos];
				 int *pathMatCurRowPosPtr=&pathMat[curRowPos];
				 int *pathMatDiagPosPtr=&pathMat[prevRowPos-1];
				 
				 int *indexMatCurRowPosPtr = &indexMat[curRowPos];

				 int *delScoreMatPrevRowPosPtr=&delScoreMat[prevRowPos];
				 int *delScoreMatCurRowPosPtr=&delScoreMat[curRowPos];
				 int *delIndexMatCurRowPosPtr=&delIndexMat[curRowPos];
				 int *delPathMatCurRowPosPtr=&delPathMat[curRowPos];

				 int *insScoreMatInsPosPtr=&insScoreMat[curRowPos-1];
				 int *insScoreMatCurRowPosPtr=&insScoreMat[curRowPos];
				 int *insIndexMatCurRowPosPtr=&insIndexMat[curRowPos];
				 int *insPathMatCurRowPosPtr=&insPathMat[curRowPos];

				 char tChar=tSeq[ti+tStart];
				 char *qCharPtr=&qSeq[1+qS[ti]];
				 int qEPrev = qE[ti-1];
				 int qECur  = qE[ti];
				 int qSCur  = qS[ti];

				 for (int qi=1; qi < rowEnd; qi++, curRowPos++, prevRowPos++, ++scoreMatPrevRowPosPtr, ++scoreMatCurRowPosPtr, ++scoreMatDiagPosPtr, ++scoreMatInsPosPtr, ++pathMatPrevRowPosPtr, ++pathMatCurRowPosPtr, ++pathMatDiagPosPtr, ++delScoreMatPrevRowPosPtr,++delScoreMatCurRowPosPtr, ++insScoreMatInsPosPtr, ++insScoreMatCurRowPosPtr, ++delIndexMatCurRowPosPtr, ++insIndexMatCurRowPosPtr, ++qCharPtr, ++delPathMatCurRowPosPtr, ++insPathMatCurRowPosPtr, ++indexMatCurRowPosPtr) {
					 int matchScore, insScore, delScore;
					 int matchIndex, insIndex, delIndex;
					 int delOpenScore, delExtendScore;
					 int insOpenScore, insExtendScore;
					 //
					 // Del matrix (down)
					 //
					 //					 if (qE[ti-1] >= qi+qS[ti] && pathMat[prevRowPos] != bound) {
					 if (qE[ti-1] >= qi+qS[ti] && *pathMatPrevRowPosPtr != bound) {
						 //						 delOpenScore = scoreMat[prevRowPos] + gapOpen;
						 delOpenScore = *scoreMatPrevRowPosPtr + gapOpen;
						 
						 //						 delExtendScore = delScoreMat[prevRowPos] + gapExtend;
						 delExtendScore = *delScoreMatPrevRowPosPtr + gapExtend;
						 //						 assert(pathMat[prevRowPos] != bound);
					 }
					 else {
						 delOpenScore = BAD;
						 delExtendScore = BAD;
					 }
					 int maxScore=max(delOpenScore, delExtendScore);
					 if (maxScore == delOpenScore) {						 
						 //						 delPathMat[curRowPos] = delOpen;
						 *delPathMatCurRowPosPtr = delOpen;
						 //						 delIndexMat[curRowPos] = prevRowPos;
						 *delIndexMatCurRowPosPtr = prevRowPos;
					 }
					 else {
						 //						 delPathMat[curRowPos] = delExtend;
						 *delPathMatCurRowPosPtr = delExtend;
						 
						 //						 delIndexMat[curRowPos] = prevRowPos;
						 *delIndexMatCurRowPosPtr = prevRowPos;
					 }
					 //					 delScoreMat[curRowPos]=maxScore;
					 *delScoreMatCurRowPosPtr = maxScore;
					 //
					 // Insertion matrix.
					 //
					 //					 insOpenScore=scoreMat[curRowPos-1] + gapOpen;
					 insOpenScore = *scoreMatInsPosPtr + gapOpen;
					 

					 //					 insExtendScore=insScoreMat[curRowPos-1] + gapExtend;
					 insExtendScore = *insScoreMatInsPosPtr + gapExtend;
					 maxScore=max(insOpenScore, insExtendScore);
						 if (maxScore == insOpenScore) {							 
							 //							 insPathMat[curRowPos] = insOpen;
							 *insPathMatCurRowPosPtr = insOpen;
							 //							 insIndexMat[curRowPos]= curRowPos-1;
							 *insIndexMatCurRowPosPtr = curRowPos-1;
						 }
						 else {
							 //							 insPathMat[curRowPos] = insExtend;
							 *insPathMatCurRowPosPtr = insExtend;
							 //							 insIndexMat[curRowPos] = curRowPos-1;
							 *insIndexMatCurRowPosPtr = curRowPos-1;
						 }
						 assert(insIndexMat[curRowPos] >= 0);
						 //						 insScoreMat[curRowPos] = maxScore;
						 *insScoreMatCurRowPosPtr = maxScore;

							 

					 
					 if (qE[ti-1] >= qi + qS[ti] && pathMat[prevRowPos-1] != bound) {
						 //						 if (tSeq[ti+tStart] == qSeq[qi+qS[ti]]) {
						 if (tChar == *qCharPtr) {
							 //							 matchScore = scoreMat[prevRowPos-1] + match;
							 matchScore = *scoreMatDiagPosPtr + match;
							 //							 assert(pathMat[prevRowPos-1] != bound);
						 }
						 else { 
							 //							 matchScore = scoreMat[prevRowPos-1] + mismatch;
							 matchScore = *scoreMatDiagPosPtr + mismatch;
							 //							 assert(pathMat[prevRowPos-1] != bound);
						 }
					 }
					 else {
						 matchScore = BAD;
					 }
					 //					 insScore = scoreMat[curRowPos-1] + gap;
					 insScore = *scoreMatInsPosPtr + gap;
					 //					 if (qE[ti-1] >= qi+qS[ti] && pathMat[prevRowPos] != bound) {
					 if (qEPrev >= qi+qSCur && *pathMatPrevRowPosPtr != bound) {
						 //						 delScore = scoreMat[prevRowPos] + gap;
						 delScore = *scoreMatPrevRowPosPtr + gap;
						 //						 assert(pathMat[prevRowPos] != bound);
					 }
					 else {
						 delScore = BAD;
					 }
					 //					 int delCloseScore = delScoreMat[curRowPos];
					 int delCloseScore = *delScoreMatCurRowPosPtr;
					 //					 int insCloseScore = insScoreMat[curRowPos];
					 int insCloseScore = *insScoreMatCurRowPosPtr;
					 maxScore = max(matchScore,max(insScore, max(delScore, max(delCloseScore, insCloseScore))));
					 assert(curRowPos < scoreMat.size());
					 //					 scoreMat[curRowPos] = maxScore;
					 *scoreMatCurRowPosPtr = maxScore;
					 if (maxScore == matchScore) {
						 //						 pathMat[curRowPos]  = diag;
						 *pathMatCurRowPosPtr = diag;
						 //						 assert(prevRowPos-1 < cuMatSize[ti]);
						 //						 indexMat[curRowPos] = prevRowPos-1;
						 *indexMatCurRowPosPtr = prevRowPos-1;
					 }
					 else if (maxScore == insScore) {
						 //						 pathMat[curRowPos]  = left;
						 *pathMatCurRowPosPtr = left;
						 //						 indexMat[curRowPos] = curRowPos-1;
						 *indexMatCurRowPosPtr = curRowPos-1;
						 //						 assert(curRowPos-1 > cuMatSize[ti]);
					 }
					 else if (maxScore == delScore) {
						 //						 pathMat[curRowPos]  = down;
						 *pathMatCurRowPosPtr = down;

						 //						 assert(prevRowPos <= cuMatSize[ti]);
						 //						 indexMat[curRowPos] = prevRowPos;
						 *indexMatCurRowPosPtr = prevRowPos;
					 }
					 else if (maxScore == delCloseScore) {
						 //						 pathMat[curRowPos] = delClose;
						 *pathMatCurRowPosPtr = delClose;
						 //						 indexMat[curRowPos] = curRowPos;
						 *indexMatCurRowPosPtr = curRowPos;
					 }
					 else if (maxScore == insCloseScore) {
						 //						 pathMat[curRowPos] = insClose;
						 *pathMatCurRowPosPtr = insClose;
						 //indexMat[curRowPos] = curRowPos;
						 *indexMatCurRowPosPtr = curRowPos;
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
			 int matchMat=0, delMat=1, insMat=2;
			 int curMat=0;
			 int curMatPos=pathMat.size()-1;

			 while (curMatPos > 0) {
				 if (curMat == matchMat) {
					 
					 if (pathMat[curMatPos] == delClose) {
						 //						 cout << "Hopping to del " << endl;
						 curMat=delMat;
					 }
					 else if (pathMat[curMatPos] == insClose) {
						 //						 cout << "Hopping to ins" << endl;
						 curMat=insMat;
					 }
					 else {
						 assert(curMatPos < pathMat.size());
						 path.push_back(pathMat[curMatPos]);
						 pathMat[curMatPos] = 90+pathMat[curMatPos];
					 }
					 curMatPos=indexMat[curMatPos];
				 }
				 else if (curMat == delMat) {
					 path.push_back(down);
					 if (delPathMat[curMatPos] == delOpen) {
						 curMat=matchMat;
					 }
					 else {
						 curMat=delMat;
					 }
					 curMatPos=delIndexMat[curMatPos];
					 assert(curMatPos >= 0);
				 }
				 else {
					 assert(curMat == insMat);
					 path.push_back(left);
					 if (insPathMat[curMatPos] == insOpen) {
						 curMat = matchMat;
					 }
					 else {
						 curMat=insMat;
					 }
					 curMatPos=insIndexMat[curMatPos];
					 assert(curMatPos >= 0);
				 }					 						 
				 assert(curMatPos >= 0);
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
				 /*
				 if (r > 0 and np != 1) {
					 PrintMat(pathMat, qS, qE);
					 cout << "ERROR at " << r << "\t" << np << endl;

					 }
				 assert(r == 0 || np == 1);
				 */
			 }

			 //			 cout<< "Start " << startBlock << " end " << endBlock << endl;
			 //			 FlatPrintMat(pathMat, qS, qE);
			 //			PrintMat(scoreMat, qS, qE);
			 long qPath=qStart;
			 long tPath=tStart;
			 int pi=path.size()-1;
			 int nm=0, ni=0,nd=0;
			 for (int p=0;p<path.size();p++) {
				 assert(path[p] >= 0&& path[p] <= 2);
				 if (path[p] == 0) { nm++;}
				 if (path[p] == 1) { ni++;}
				 if (path[p] == 2) { nd++;}
			 }

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
				 if (pi < path.size()) {

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
				 }

				 refined.push_back(Block(qPath,tPath, blockLen));
				 qPath+=blockLen+qgapLen;
				 tPath+=blockLen+tgapLen;
			 }				
			
			 if (endBlock < alignment.blocks.size()-1) {
				 assert(tPath <= alignment.blocks[endBlock+1].tPos);
				 assert(qPath <= alignment.blocks[endBlock+1].qPos);
			 }
			
			 if (qPath != qSeqLen + qStart || tPath != tSeqLen + tStart) {
				 //				 cout << "ERROR on " << read.name << endl;		
				 assert(qPath == qSeqLen + qStart);
				 assert(tPath == tSeqLen + tStart);
			 }
			}
			//cerr << "ratio " << scoreMat.size() / (float)qS.size() << endl;			
		}
		endBlock++;
		startBlock=endBlock;		
	}
	alignment.blocks = refined;
	if (alignment.blocks.size() > 1) {
		for (int b=0; b < alignment.blocks.size()-1; b++) {
			if (alignment.blocks[b].qPos + alignment.blocks[b].length > alignment.blocks[b+1].qPos or
					alignment.blocks[b].tPos + alignment.blocks[b].length > alignment.blocks[b+1].tPos ) {
				cout << "ERROR with alignment consistency of " << read.name << endl;
				cout << "block " << b << endl;
				cout << "q: " << alignment.blocks[b].qPos + alignment.blocks[b].length << "\t" << alignment.blocks[b+1].qPos << endl;
				cout << "t: " << alignment.blocks[b].tPos + alignment.blocks[b].length << "\t" << alignment.blocks[b+1].tPos << endl;
				assert(0);
			}
		}
	}
}


#endif
