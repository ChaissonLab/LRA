#ifndef MIN_COUNT_H_
#define MIN_COUNT_H_
#include "TupleOps.h"
#include "SeqUtils.h"
#include "htslib/kseq.h"
KSEQ_INIT(gzFile, gzread)




template <typename TupPos, typename Tup> void StoreMinimizers(char *seq, 
																															int seqLen,
																															int k, 
																															int w,
																															vector<TupPos> &minimizers, bool canonical=true) {
	//
	// Initialize first.
	//

	if (seqLen < k) {
		return;
	}
	TupPos cur, curRC, minTuple, can;
	GenomePos minPos;

	GenomePos p = 0;
	TupPos m;
	InitMask(m, k);
	StoreTuple(seq, p, k, cur);
	TupleRC(cur, curRC, k);

	//
	// Initialize the first minimzer.
	//
	if (canonical) {
		can.t = min(cur.t, curRC.t);
	}
	else { can.t = cur.t; }
	minPos = 0;
	TupPos activeMinimizer, curMinimizer;
	activeMinimizer.t = can.t;
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
		if (canonical) {
			curMinimizer.t = min(cur.t, curRC.t);
		}
		else {
			curMinimizer.t = cur.t;
		}

		
		if (curMinimizer.t < activeMinimizer.t) {
			activeMinimizer.t = curMinimizer.t;
			activeMinimizer.pos = p;
		}
		curTuples[p%w] = curMinimizer;
	}
	minimizers.push_back(activeMinimizer);
	// Now scan the chromosome
	minTuple.t=m.t;
	for (p=w; p < seqLen-k+1; p++) {
		
		// Check if past current active minimzier

		ShiftOne(seq, p+k-1, m, cur);
		ShiftOneRC(seq, p+k-1, k, curRC);
#ifdef _TESTING_
		TupPos test, testrc;
		StoreTuple(seq, p, k, test);
		TupleRC(test, testrc, k);

		assert(test.t == cur.t);
		assert(testrc.t == curRC.t);
#endif
		if (canonical) {
			curMinimizer.t = min(cur.t, curRC.t);
		}
		else {
			curMinimizer.t = cur.t;
		}
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

#endif

