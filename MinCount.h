#ifndef MIN_COUNT_H_
#define MIN_COUNT_H_
#include "TupleOps.h"
#include "SeqUtils.h"
#include "htslib/kseq.h"
KSEQ_INIT(gzFile, gzread)



template <typename TupPos> void SimpleMinimizers(kseq_t *seq, int k, int w, vector<TupPos> &minimizers) {


	Tuple cur, curRC, minTuple, can;
	int minPos;
	if (seq->seq.l < k) {
		return;
	}
	int p = 0;
	int curMinimizer = -1;

	for (p = 0; p < seq->seq.l - w - k + 1; p++) {

		int i;
		minTuple = mask;
		int minPos = 0;
		for (i = p; i < p+w; i++) {
			StoreTuple(seq->seq.s, i, k, cur);
			TupleRC(cur, curRC, k);
			can = min(cur, curRC);
			if (can < minTuple) {
				minPos = i;
				minTuple = can;
			}
		}
		if (minPos != curMinimizer) {
			assert(minPos > curMinimizer);
			GenomeTuple tup;
			tup.t=minTuple;
			tup.pos = minPos;
			minimizers.push_back(tup);
			curMinimizer = minPos;
		}
	}
}

template <typename TupPos, typename Tup> void StoreMinimizers(char *seq, 
																															int seqLen,
																															int k, 
																															int w,
																															vector<TupPos> &minimizers) {
	//
	// Initialize first.
	//

	if (seqLen < k) {
		return;
	}
	Tup cur, curRC, minTuple, can;
	int minPos;

	int p = 0;
	Tup m;
	InitMask(m, k);
	StoreTuple(seq, p, k, cur);
	TupleRC(cur, curRC, k);

	//
	// Initialize the first minimzer.
	//
	can = min(cur, curRC);
	minPos = 0;
	TupPos activeMinimizer, curMinimizer;
	activeMinimizer.t = can;
	activeMinimizer.pos = 0;
	//	priority_queue<GenomeTuple, vector<GenomeTuple>, GenomeTupleComp > pQueue;
	vector<TupPos> curTuples(w);
	curTuples[0] = activeMinimizer;

	// 
	// Find the active minimizer in this window
	int nMinimizers=1;

	for (p = 1; p< w && p < seqLen-k+1 ; p++) {
		ShiftOne(seq, p+k-1, m, cur);
		ShiftOneRC(seq, p+k-1, k, curRC);
		/*
		Tuple test, testrc;
		StoreTuple(seq->seq.s, p, k, test);
		TupleRC(test, testrc, k);
		assert(test == cur);
		assert(testrc == curRC);
		*/
		curMinimizer.pos   = p;
		curMinimizer.t = min(cur, curRC);
		
		if (curMinimizer.t < activeMinimizer.t) {
			activeMinimizer.t = curMinimizer.t;
			activeMinimizer.pos = p;
		}
		curTuples[p%w] = curMinimizer;
	}
	minimizers.push_back(activeMinimizer);
	// Now scan the chromosome
	minTuple=mask;
	for (p=w; p < seqLen-k+1; p++) {
		
		// Check if past current active minimzier

		ShiftOne(seq, p+k-1, m, cur);
		ShiftOneRC(seq, p+k-1, k, curRC);
#ifdef _TESTING_
		Tuple test, testrc;
		StoreTuple(seq, p, k, test);
		TupleRC(test, testrc, k);

		assert(test == cur);
		assert(testrc == curRC);
#endif

		curMinimizer.t = min(cur, curRC);
		curMinimizer.pos   = p;

		curTuples[p%w] = curMinimizer;
		if (p - w >= activeMinimizer.pos) {
			activeMinimizer = curTuples[0];
			for (int j =1; j < w; j++) {
				if (curTuples[j].t < activeMinimizer.t) {
					activeMinimizer = curTuples[j];
				}
			}
			minimizers.push_back(activeMinimizer);
			nMinimizers+=1;
		}
		else {
			if (curMinimizer.t < activeMinimizer.t) {
				activeMinimizer = curMinimizer;
				minimizers.push_back(activeMinimizer);
				nMinimizers++;		
			}
		}		
		if (p + 1 % 10000 == 0) {
			cerr << p +1 << endl;
		}
	}
}

template <typename TupPos> void StoreAll(kseq_t *seq, int k, vector<TupPos> &tuples) {
	//
	// Initialize first.
	//

	Tuple cur, curRC, minTuple, can;
	int minPos;
	if (seq->seq.l < k) {
		return;
	}
	int p = 0;
	StoreTuple(seq->seq.s, p, k, cur);
	TupleRC(cur, curRC, k);

	//
	// Initialize the first minimzer.
	//
	can = min(cur, curRC);
	GenomeTuple gt;
	gt.t = can;
	gt.pos = 0;

	tuples.push_back(gt);
	Tuple m;
	InitMask(m,  k);
	for (p = 1; p < seq->seq.l - k + 1; p++) {
		ShiftOne(seq->seq.s, p+k-1, m, cur);
		ShiftOneRC(seq->seq.s, p+k-1, k, curRC);
		can=min(cur,curRC);
		//	priority_queue<GenomeTuple, vector<GenomeTuple>, GenomeTupleComp > pQueue;
		gt.t = can;
		gt.pos = p;

		tuples.push_back(gt);
	}
}	



#endif

