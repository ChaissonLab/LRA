#ifndef AFFINE_ONE_GAP_ALIGN_H_
#define AFFINE_ONE_GAP_ALIGN_H_

#include <limits.h>
#include "SeqUtils.h"
#include "Alignment.h"

// i, j index from string space into matrix space. 
// The strings are 1-based, and so they range from [1,qLen], [1,tLen]
// The optimization resolves from  [d-k, d+k]
// The matrix holds from [
int PreToIndex(int i, int j, int k, int band) {
	int d=i-j;
	// The extra + 1 is because of the rail on the side of each k-band
	return j*band+d+k+1;
}

int SuffToIndex(int ii, int jj, int is, int js, int k, int band) {
	//	return jj*band + ii;
	int i = ii - is;
	int j = jj - js;
	int d=i-j;
	//	cout << "sti " << ii << " " << jj << " " << i << " " << j << " " << d << " " << j*band+d+k+1 << endl;
	int res= j*band+d+k+1;
	return res;
}

#define MISSING -3200000
template<typename T>
void PrintMat(string &qSeq, string &tSeq, vector<T> &mat, int qLen, int tLen, int k, int R, int w) {
	vector<int> row(qLen+1);
	int diag=min(qLen,tLen);
	cout << "  ";
	cout.width(w);
	cout<< "s";
	for (int ii=0; ii < qLen+1; ii++) {
		cout.width(w);
		if (ii == 0) {
			cout << "-";
		}
		else {
			cout << qSeq[ii-1];
		}
	}
	cout << endl;
	
	cout << "  ";
	cout.width(w);
	cout << "p";
	for (int ii=0; ii < qLen+1; ii++) {
		cout.width(w);
		cout << ii;
	}
	cout << endl;
	int ri=0;
	for (int jj=0; jj < diag+k;jj++) {
		fill(row.begin(), row.end(), -1);
		for (int ii=max(0,jj-k-1); ii < min(qLen+1,jj+k+2); ii++) {
			row[ii] = mat[PreToIndex(ii,jj,k,R)];
		}
		if (jj == 0) {
			cout << "- ";
		}
		else {
			cout << tSeq[jj-1] << " ";
		}
		cout.width(w);
		cout <<ri;
		++ri;
		for (int r=0; r< row.size(); r++) {
			cout.width(w);
			cout << row[r];
		}
		cout << endl;
	}
}
template <typename T>
void PrintSuffMat(string q, string t,
									vector<T > &mat, int qLen, int tLen, int qLow, int tLow, int k, int R, int w=3) {
	vector<int> row(qLen+1);
	int diag=min(qLen,tLen);
	cout << "  ";
	cout.width(w);
  cout << "+";
	cout.width(w);
	cout << "-";
	for (int ii=0; ii < qLen; ii++) {
		cout.width(w);
		cout << q[ii];
	}
	cout << endl;
	cout << "  ";
	cout.width(w);
	cout<< "p";
	for (int ii=0; ii < qLow; ii++) {
		cout.width(w);
		cout << " ";
	}
	
	for (int ii=max(qLow,0); ii < qLen+1; ii++) {
		cout.width(w);
		cout << ii;
	}
	cout << endl;
	int ri=0;
	
	for (int jj=tLow; jj < tLen+1;jj++) {
		fill(row.begin(), row.end(), -1);
		int doff=diag-(tLen-jj);
		if (jj > 0) {
			cout << t[jj-1] << " ";
		}
		else {
			cout << "  ";
		}
		for (int ii=max(0,max(qLow,qLow+doff-k-1)); ii < min(qLen+1,qLow+doff+k+2); ii++) {
		//		for (int ii =0; ii < qLen; ii++) {
			assert(SuffToIndex(ii,jj,qLow, tLow, k,R) < mat.size());
			assert(SuffToIndex(ii,jj,qLow, tLow, k,R) >=0);
			assert(ii< row.size());
			row[ii] = mat[SuffToIndex(ii,jj,qLow, tLow, k,R)];
		}
		cout.width(w);
		cout << jj;
		++ri;
		for (int r=0; r< row.size(); r++) {
			cout.width(w);
			cout << row[r];
		}
		cout << endl;
	}
}

int AffineOneGapAlign(string &qSeq, string &tSeq, int m, int mm, int indel, int k,  Alignment &aln) 
{
	int qLen = qSeq.size();
	int tLen = tSeq.size();
	int diag = min(qLen, tLen);
	int doneAr=0;
	int leftAr=1;
	int downAr=2;
	int diagAr=3;
	int borderAr=4;
	int gapLeftAr=5;
	int gapDownAr=6;


	vector<int> qInt;
	qInt.resize(qLen+1,0);
	vector<int> tInt;
	tInt.resize(tLen+1,0);
	for (int s=0; s < qLen; s++) {
		qInt[s+1] = seqMapN[qSeq[s]];
	}
	for (int s=0; s < tLen; s++) {
		tInt[s+1] = seqMapN[tSeq[s]];
	}
	vector<int> upperDiagonalMax(diag+1, MISSING), 
		upperDiagonalIndex(diag+1, 0), 
		lowerDiagonalMax(diag+1, MISSING), 
		lowerDiagonalIndex(diag+1, borderAr);

	k = min(diag, k);
	bool alignTop = true;
	if (diag + 2*k > max(qLen, tLen)) {
		k=max(qLen,tLen);
		alignTop= false;
	}
		
		
		
	int band=k*2+1; // k on each side, and middle 
	
	int R=band+2;
	int matSize = (3+k+diag)*R;
	vector<long> pScore(matSize, MISSING);
	vector<int>	 pPath(matSize,-1);
	vector<long> sScore(matSize, MISSING); //(qLen+1)*(tLen+1), -1);
	vector<int>	 sPath(matSize, -1); //(qLen+1)*(tLen+1),-1);
	//	vector<long> sScore((max(R,qLen+2))*(tLen+2), -1);
	//	vector<int>	 sPath((max(R,qLen+2))*(tLen+2),-1);

	//
	// First fill out the prefix matrices
	//
	//	sPath[PreToIndex(0,0,k,R)] = doneAr;
	int i,j;
	//	cout << "Lower bound " << endl;
	for (i=1; i < k+1; i++) {
		assert(i+k+1 < sScore.size());
		pScore[PreToIndex(i,0,k,R)] = indel*i;
		pPath[PreToIndex(i,0,k,R)] = leftAr;
	}
	//	cout << "left bound" << endl;
	for (j=1; j <= k+1; j++) {
		assert(PreToIndex(0, j, k, R) < sScore.size());
		pScore[PreToIndex(0, j, k, R) ] = indel*j;
		pPath[PreToIndex(0, j, k, R)] = downAr;
	}
	pScore[PreToIndex(0,0,k,R)] = 0;
	pPath[PreToIndex(0,0,k,R)] = doneAr;
		
#ifdef _MAT_PRINT_
	cout <<"k-boundaries " << endl;
	PrintMat(qSeq, tSeq, pPath, qLen, tLen, k, R,4);
	PrintMat(qSeq, tSeq, pScore, qLen, tLen, k, R,4);
#endif
	if (qLen >= tLen) {
		// Left diagonal is a rail
		// If not rail, then store colmax
		//		cout << "lower diagonal " << endl;
		for (i=0; i <= diag-k-1; i++) {
			assert(PreToIndex(i,i+k+1,k,R) < pScore.size());
			pScore[PreToIndex(i,i+k+1,k,R)] = MISSING;
			pPath[PreToIndex(i,i+k+1,k,R)] = borderAr;
#ifdef _MAT_PRINT_
			cout << "init upper " << i << endl; 
			PrintMat(qSeq, tSeq, pPath, qLen, tLen, k, R,4);
			PrintMat(qSeq, tSeq, pScore, qLen, tLen, k, R,4);
#endif

		}
		for (i=1; i < diag+k-1; i++) {
			assert(PreToIndex(i+k+1,i,k,R) < pScore.size());
			pScore[PreToIndex(i+k+1,i,k,R)] = MISSING;
			pPath[PreToIndex(i+k+1,i,k,R)] = borderAr;
#ifdef _MAT_PRINT_
			cout << "init upper " << i << endl; 
			PrintMat(qSeq, tSeq, pPath, qLen, tLen, k, R,4);
			PrintMat(qSeq, tSeq, pScore, qLen, tLen, k, R,4);
#endif

		}

		lowerDiagonalMax[0] = 0;
		lowerDiagonalIndex[0] = 0;
	}
#ifdef _MAT_PRINT_
  cout << "initialized left " << endl;
	PrintMat(qSeq, tSeq, pPath, qLen, tLen, k, R,4);
	PrintMat(qSeq, tSeq, pScore, qLen, tLen, k, R,4);
#endif
			
	if (qLen <= tLen) {
		// Right diagonal is a rail.
		//		cout << "qlen < tlen " << endl;
		//		cout << "upper diagonal "<< endl;
		for (j=0; j < diag-1; j++) {
			assert(PreToIndex(j+k+1,j,k,R) < pScore.size());
			pScore[PreToIndex(j+k+1,j,k,R)] = MISSING;
			pPath[PreToIndex(j+k+1,j,k,R)] = borderAr;

		}
		for (j=1; j < diag+k; j++) {
			assert(PreToIndex(j-k-1,j,k,R) < pScore.size());
			pScore[PreToIndex(j-k-1,j,k,R)] = MISSING;
			pPath[PreToIndex(j-k-1,j,k,R)] = borderAr;
		}

		upperDiagonalMax[0] = 0;
		upperDiagonalIndex[0] = 0;
#ifdef _MAT_PRINT_
	PrintMat(qSeq, tSeq, pPath, qLen, tLen, k, R,4);
	PrintMat(qSeq, tSeq, pScore, qLen, tLen, k, R,4);
#endif
	}

	int d;
	int qBoundary = min(diag+k, qLen+1);
	int tBoundary = min(diag+k, tLen+1);
	int diagBoundary = max(qBoundary, tBoundary);

	for (j=1; j < tBoundary; j++ ) { // j = diagonal
		for (i = max(1,j-k); i < min(qBoundary, j+k+1); i++) {

			assert(PreToIndex(i,j-1,k,R) < pScore.size());
			assert(PreToIndex(i,j-1,k,R)>=0);
			long sIns = pScore[PreToIndex(i-1,j,k,R)] + indel;
			long sDel = pScore[PreToIndex(i,j-1,k,R)] + indel;
			long sMat;
			if (qInt[i] == tInt[j]) {
				sMat = pScore[PreToIndex(i-1,j-1,k, R)] + m;
			}
			else {
				sMat = pScore[PreToIndex(i-1,j-1,k,R)] + mm;
			}

			long maxScore = max(sIns,max(sDel,sMat));
			pScore[PreToIndex(i,j,k,R)] = maxScore;
			if (maxScore == sIns) {
				pPath[PreToIndex(i,j,k,R)] = leftAr;
			}
			else if (maxScore == sDel) {
				pPath[PreToIndex(i,j,k,R)] = downAr;
			}
			else {
				pPath[PreToIndex(i,j,k,R)] = diagAr;
			}

#ifdef _MAT_PRINT_
			cout <<"iter  " << i << " " << j << endl;
			PrintMat(qSeq, tSeq, pPath, qLen, tLen, k, R,4);
			PrintMat(qSeq, tSeq, pScore, qLen, tLen, k, R,4);
#endif

			if (i < qLen - k) {
				if (pScore[PreToIndex(i,j,k,R)] >= lowerDiagonalMax[j]) {
					lowerDiagonalMax[j] = pScore[PreToIndex(i,j,k,R)];
					lowerDiagonalIndex[j] = i;
				}
			}
			if (j < tLen and i < diag+1) {
				assert(PreToIndex(i,j,k,R) < pScore.size());
				assert(i < upperDiagonalMax.size());
				if (pScore[PreToIndex(i,j,k,R)] > upperDiagonalMax[i]) {
					upperDiagonalMax[i] = pScore[PreToIndex(i,j,k,R)];
					upperDiagonalIndex[i] = j;
				}
			}
		}		
	}
#ifdef _MAT_PRINT_
	cout << "Prefix " << endl;
	PrintMat(qSeq, tSeq, pPath, qLen, tLen, k, R,3);
	cout <<"  ";
	cout.width(4);
	cout << "udi";
	for (i=0;i<upperDiagonalIndex.size();i++) {
		cout.width(4);
		cout << upperDiagonalIndex[i];
	}
	cout << endl;
	cout <<"  ";
	cout.width(4);

	cout << "ldi";
	for (i=0;i<lowerDiagonalIndex.size();i++) {
		cout.width(4);
		cout << lowerDiagonalIndex[i];
	}
	cout << endl;
	cout <<"  ";
	cout.width(4);
	cout << " ";
	for (i=0;i<lowerDiagonalMax.size();i++) {
		cout.width(4);
		cout << lowerDiagonalMax[i];
	}
	cout << endl;

	PrintMat(qSeq, tSeq, pScore, qLen, tLen, k, R,4);
#endif

	vector<int> lengths;
	vector<int> ops;

	int maxAlnScore=-1;
	if (alignTop) {
		//
		// This is just standard affine alignment
		//


		//
		// First fill out the suffix matrices
		//
	
		// Boundary conditions, qStart is the position of the first row in the matrix
		int qStart = max(0, qLen - diag);
		int qEnd   = qLen+1;

		// tStart includes 0, so sequence alignment starts at tStart+1
		int tStart = max(0, tLen - diag );
		int tLow   = max(0, tLen - diag - k - 1 - 1);
		int qLow   = max(0, qLen - diag - k - 1);
		int tEnd   = tLen + 1;
		//		int R=qLen+1;
		int X=qLen+1;
		if (qLen >= tLen) {
			assert(tStart == 0);
			for (i = qLow, j=0; i < qStart+k+1; i++) {
				sScore[SuffToIndex(i, j, qLow, tLow, k, R)] = lowerDiagonalMax[j];
				sPath[SuffToIndex(i, j, qLow, tLow, k, R)] = gapLeftAr;
			}		

			for (i = qLow, j=1; i < qLow+diag; i++, j++) {
				sScore[SuffToIndex(i, j, qLow, tLow, k, R)] = lowerDiagonalMax[j];
				sPath[SuffToIndex(i, j, qLow, tLow, k, R)] = gapLeftAr;
			}		

			for (j=tStart+1, i=qStart; j < tEnd-k; i++,j++) {
	      sScore[SuffToIndex(i+k+1, j, qLow, tLow, k, R)] = MISSING;
				sPath[SuffToIndex(i+k+1, j, qLow, tLow, k, R)] = borderAr;			
			}
#ifdef _MAT_PRINT_
		PrintSuffMat(qSeq, tSeq, sPath, qLen, tLen, qLow, tLow, k, R);
		PrintSuffMat(qSeq, tSeq, sScore, qLen, tLen, qLow, tLow, k, R, 4);
#endif
			
		}
		if (qLen <= tLen) {
			assert(qStart == 0);
			// Init y axis
			//			cout << "Setting x axis upper grid " << endl;
			for (j = tLow, i=qStart; j < tStart + k+2; j++) {
				//				cout << " x uppergrid " << i << " " << j << endl;
				assert(SuffToIndex(i, j, qLow, tLow, k, R) < sScore.size());
				sScore[SuffToIndex(i, j, qLow, tLow, k, R)] = upperDiagonalMax[0];
				sPath[SuffToIndex(i, j, qLow, tLow, k, R)] = gapDownAr;
			}			
			// Init bottom diagonal to gap close
			//			cout << "Scoring lower diagonal upper grid" << endl;
			for (j = tStart+1, i = qStart+1; j < tEnd; i++, j++) {
				assert(SuffToIndex(i,j-k-1,qLow, tLow, k, R) < sScore.size());
				assert(i < upperDiagonalMax.size());
				sScore[SuffToIndex(i,j-k-1,qLow, tLow, k, R)] = upperDiagonalMax[i];
				sPath[SuffToIndex(i,j-k-1,qLow, tLow, k, R)] = gapDownAr;
			}
			// Init top diagonal to boundary
			//			cout << "Top diaonal border upper grid" << endl;
			for (j=tStart,i=qStart; j< tEnd-k-1; i++,j++) {
				//				cout << " x uppergrid upper diag " << i << " " << j -k+1<< endl;
				assert(SuffToIndex(i,j+k+1,qLow, tLow, k, R) < sScore.size());
				sScore[SuffToIndex(i,j+k+1,qLow, tLow, k, R)]= MISSING;
				sPath[SuffToIndex(i,j+k+1,qLow, tLow, k, R)] = borderAr;
			}
		}
#ifdef _MAT_PRINT_
		cout << "boundaries "<<endl;
		PrintSuffMat(qSeq, tSeq, sPath, qLen, tLen, qLow, tLow, k, R);
		PrintSuffMat(qSeq, tSeq, sScore, qLen, tLen, qLow, tLow, k, R, 4);
#endif

		for (j=tLow+1; j < tEnd; j++) {
			int doff=diag + 1 - (tEnd - j);
			for (i=max(qLow+1, qStart + doff  - k); i < min(qEnd, qStart + doff+k +1); i++) {
				long delClose=MISSING;
				long insClose=MISSING;
				if (qLen >= tLen) {
					assert(tStart == 0);
					assert(j < lowerDiagonalMax.size());
					delClose=lowerDiagonalMax[j];
				}
				if (tLen > qLen) {
					assert(qStart == 0);
					assert(i<upperDiagonalMax.size());
					insClose=upperDiagonalMax[i];
				}
				assert(SuffToIndex(i,j-1,qLow, tLow, k,R) < sScore.size());
				long sIns = sScore[SuffToIndex(i-1,j,qLow, tLow, k,R)] + indel;
				long sDel = sScore[SuffToIndex(i,j-1,qLow, tLow, k,R)] + indel;
				long sMat;
				if (qInt[i] == tInt[j]) {
					sMat = sScore[SuffToIndex(i-1,j-1,qLow, tLow, k, R)] + m;
				}
				else {
					sMat = sScore[SuffToIndex(i-1,j-1,qLow, tLow, k,R)] + mm;
				}
				long maxScore = max(delClose, max(insClose, max(sIns, max(sDel, sMat))));
				assert(SuffToIndex(i,j,qLow, tLow, k,R) < sScore.size());
				sScore[SuffToIndex(i,j,qLow, tLow, k,R)] = maxScore;
				if (maxScore == sIns) {
					sPath[SuffToIndex(i,j,qLow, tLow, k,R)] = leftAr;
				}
				else if (maxScore == sDel) {
					sPath[SuffToIndex(i,j,qLow, tLow, k,R)] = downAr;
				}
				else if (maxScore == sMat) {
					sPath[SuffToIndex(i,j,qLow, tLow, k,R)] = diagAr;
				}
				else if (maxScore == delClose) {
					sPath[SuffToIndex(i,j,qLow, tLow, k,R)] = gapLeftAr;
				}
				else if (maxScore == insClose) {
					sPath[SuffToIndex(i,j,qLow, tLow, k,R)] = gapDownAr;
				}
			}
		}
#ifdef _MAT_PRINT_
		PrintSuffMat(qSeq, tSeq, sPath, qLen, tLen, qLow, tLow, k, R);
		PrintSuffMat(qSeq, tSeq, sScore, qLen, tLen, qLow, tLow, k, R, 4);
#endif
		i=qLen;
		j=tLen;
		int arrow=sPath[SuffToIndex(i,j,qLow,tLow,k,R)];
	  maxAlnScore=sScore[SuffToIndex(i,j,qLow,tLow,k,R)];
		assert(arrow >= 0);
		//		cout << arrow << endl;
		while (arrow != doneAr and 
					 arrow != gapDownAr and arrow != gapLeftAr and 
					 i >= 0 and j >=0) {
			if (lengths.size() == 0 or ops[ops.size()-1] != arrow) {
				lengths.push_back(1);
				ops.push_back(arrow);
			}
			else {
				lengths[lengths.size()-1]++;
			}

			if (arrow == diagAr) {
				//				cout << i << "\t" << j << " diag "  << endl;
				i--;
				j--;
			}
			else if (arrow == leftAr) {
				//				cout << i << "\t" << j << " down "  << endl;			
				i--;
			}
			else if (arrow == downAr) {
				//				cout << i << "\t" << j << " left "  << endl;			
				j--;
			}
			else if (arrow == gapLeftAr) {
				//				cout << i << "\t" << j << " gap left "  << endl;			
				break;
				i--;
			}
			else if (arrow == gapDownAr) {
				//				cout << i << "\t" << j << " gap down "  << endl;			
				break;
				j--;
			}
			if (i >= 0 and j >= 0) {
				arrow=sPath[SuffToIndex(i,j,qLow,tLow,k,R)];		
			}
			assert(arrow >= 0);
		}
		if (arrow == gapDownAr) {
			
			lengths.push_back(j-upperDiagonalIndex[i]);
			ops.push_back(arrow);
			j=upperDiagonalIndex[i];
			//			cout << "will continue on " << j << endl;
		}
		if (arrow == gapLeftAr) {
			lengths.push_back(i-lowerDiagonalIndex[j]);
			ops.push_back(arrow);
			i=lowerDiagonalIndex[j];
			//			cout << "will continue on query " << i << "," << j << endl;
		}
	}
	else {
		i=qLen;
		j=tLen;
		maxAlnScore=pScore[PreToIndex(i,j,k,R)];
	}
	
	int arrow=pPath[PreToIndex(i,j,k,R)];
	
	//	cout << arrow << endl;
	while (arrow != borderAr and arrow != doneAr and i >= 0 and j >= 0) {
		assert(arrow != -1);
		assert( arrow != gapDownAr and arrow != gapLeftAr);
		//
		// start gap.
		if (lengths.size() == 0 or ops[ops.size()-1] != arrow) {
			lengths.push_back(1);
			ops.push_back(arrow);
		}
		else {
			lengths[lengths.size()-1]++;
		}
		
		if (arrow == diagAr) {
			//			cout << i << "\t" << j << " pre diag "  << endl;
			i--;
			j--;
		}
		else if (arrow == leftAr) {
			//			cout << i << "\t" << j << " pre left "  << endl;			
			i--;
		}
		else if (arrow == downAr) {
			//			cout << i << "\t" << j << " pre down "  << endl;			
			j--;
		}
		else if (arrow == gapLeftAr) {
			//			cout << i << "\t" << j << " pre gap left "  << endl;			
			break;
		}
		else if (arrow == gapDownAr) {
			//			cout << i << "\t" << j << " gap down "  << endl;			
			break;
		}
		arrow=pPath[PreToIndex(i,j,k,R)];		
	}
	int qPos=0;
	int tPos=0;
	for (i=lengths.size(); i > 0; i--) {
		int op=ops[i-1];
		int len=lengths[i-1];

		if (op == leftAr or op == gapLeftAr) {
			qPos+= len;
		}
		else if (op == downAr or op == gapDownAr) {
			tPos+= len;
		}
		else if (op == diagAr) {
			aln.blocks.push_back(Block(qPos, tPos, len));
			qPos+=len;
			tPos+=len;
		}
	}
	return maxAlnScore;
}		
	

#endif
