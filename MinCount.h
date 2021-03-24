#ifndef MIN_COUNT_H_
#define MIN_COUNT_H_
#include "TupleOps.h"
#include "SeqUtils.h"
#include "htslib/kseq.h"

template <typename TupPos, typename Tup> 
void StoreMinimizers(char *seq, GenomePos seqLen, int k, int w, vector<TupPos> &minimizers, bool Global, bool canonical = true) {
	//
	// Initialize first.
	//
	if (seqLen < k) {
		return;
	}
	TupPos cur, curRC, minTuple, can;
	GenomePos minPos;
	int windowSpan=w+k-1;
	GenomePos p = 0;
	TupPos m;
	//
	// Skip N's as a start
	InitMask(m, k);
	int nextValidWindowEnd=0;
	int nextValidWindowStart=0;
	bool valid=false;
	while (nextValidWindowStart < seqLen - windowSpan and !valid) {
		valid=true;
		for (int n=nextValidWindowStart; valid and n < nextValidWindowStart+windowSpan; n++ ) {
			if (seqMapN[seq[n]] > 3) {
				nextValidWindowStart = n+1;
				valid=false;
			}
		}
	}
	// all n
	if (valid == false) {
		return;
	}
	nextValidWindowEnd = nextValidWindowStart + windowSpan;

	StoreTuple(seq,p,k,cur);
	TupleRC(cur, curRC, k);
	//
	// Initialize the first minimzer.
	// Store canonical information in the rightest bit of can.t; The last bit = 1 ==> reverse strand
	// If canonical == True, for_mask = 0111...11 --> minimizer & for_mask = 0minimizer; rev_mask = 1000...00 --> minimizer | rev_mask = 1minimizer
	// Else for_mask = 111...11 --> minimizer & for_mask = minimizer, rev_mask = 000...00 --> minimizer | rev_mask = minimizer
	//
	Tup for_mask = TupPos::for_mask_s;	
	Tup rev_mask = TupPos::rev_mask_s;
	Tup mask = 0;
	if (!canonical and Global) {
		rev_mask = (rev_mask & mask); // 0000...00 64 bits
		for_mask = ~rev_mask; // 111...11 64 bits
	}

	if (canonical) { 
		if ((cur.t & for_mask) < (curRC.t & for_mask)) can.t = (cur.t & for_mask); //can.t = min(cur.t, curRC.t);
		else can.t = (curRC.t | rev_mask); 
	}
	else { can.t = cur.t; }
	minPos = 0;
	TupPos activeMinimizer, curMinimizer;
	activeMinimizer.t = can.t;
	activeMinimizer.pos = 0;
	vector<TupPos> curTuples(w);
	curTuples[0] = activeMinimizer;

	// 
	// Find the active minimizer in this window
	//
	int nMinimizers=1;

	for (p = 1; p < w && p < seqLen-k+1 ; p++) {

		ShiftOne(seq, p+k-1, m, cur);
		ShiftOneRC(seq, p+k-1, k, curRC);
		/*
		Tuple test, testrc;
		StoreTuple(seq->seq.s, p, k, test);
		TupleRC(test, testrc, k);
		assert(test == cur);
		assert(testrc == curRC);
		*/
		curMinimizer.pos = p;
		if ((cur.t & for_mask) < (curRC.t & for_mask)) curMinimizer.t = (cur.t & for_mask);
		else curMinimizer.t = (curRC.t | rev_mask); 
		if (curMinimizer.t < activeMinimizer.t) {  
			activeMinimizer.t = curMinimizer.t;
			activeMinimizer.pos = p;
		}	
		curTuples[p%w] = curMinimizer;
	}
	//
	// Only store the first minimizer if the first window starts at the beginning of the sequence.
	//
	if (nextValidWindowEnd == windowSpan ) {
		minimizers.push_back(activeMinimizer);
	}
	// Now scan the chromosome
	minTuple.t=m.t;
	for (p = w; p < seqLen-k+1; p++) {
		// If the next valid window ends at the next nucleotide, check to see if 
		// it is a valid window (no N's). If so, bump by one.
		// Otherwise, search for the next valid window end.
		if ( nextValidWindowEnd == p+k-1)  {
			if ( seqMapN[seq[p+k-1]] <= 3 ) {
				nextValidWindowEnd++;
			}
			else {
				nextValidWindowStart = p+k;
				valid=false;			
				while (nextValidWindowStart < seqLen - windowSpan and not valid) {
					valid=true;
					for (int n=nextValidWindowStart; valid and n < nextValidWindowStart+windowSpan; n++ ) {
						if (seqMapN[seq[n]] > 3) {
							nextValidWindowStart = n+1;
							valid=false;
						}
					}			
					// all n
					if (valid == false) {
						return;
					}
					nextValidWindowEnd = nextValidWindowStart + windowSpan;
				}
			}
		}		

		ShiftOne(seq, p+k-1, m, cur);
		ShiftOneRC(seq, p+k-1, k, curRC);
#ifdef _TESTING_
		TupPos test, testrc;
		StoreTuple(seq, p, k, test);
		TupleRC(test, testrc, k);

		assert(test.t == cur.t);
		assert(testrc.t == curRC.t);
#endif
		if ((cur.t & for_mask) < (curRC.t & for_mask)) curMinimizer.t = (cur.t & for_mask);
		else curMinimizer.t = (curRC.t | rev_mask); 
		curMinimizer.pos = p;
		curTuples[p%w] = curMinimizer;
		if (p - w >= activeMinimizer.pos) {
			activeMinimizer = curTuples[0];
			for (int j =1; j < w; j++) {
				if ((curTuples[j].t & for_mask) < (activeMinimizer.t & for_mask)) { 
					activeMinimizer = curTuples[j];
				}		
			}
			if (nextValidWindowEnd == p+k) {
				minimizers.push_back(activeMinimizer);			
				nMinimizers+=1;
			}
			else {
				cout << "skipping at " << p << endl;
			}
		}
		else {
			if ((curMinimizer.t & for_mask) < (activeMinimizer.t & for_mask)) { //TODO(Jingwen)
				activeMinimizer = curMinimizer;
				if (nextValidWindowEnd == p+k) {
					minimizers.push_back(activeMinimizer);
					nMinimizers++;
				}
				else {
					cout << "dotopart skipping at " << p << endl;
				}
			}		
		}		
		if (p + 1 % 10000 == 0) {
			cerr << p +1 << endl;
		}
	}
}

template <typename TupPos, typename Tup> 
void StoreMinimizers_noncanonical(char *seq, GenomePos seqLen, int k, int w, vector<TupPos> &minimizers, bool Global) {
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
	// TupleRC(cur, curRC, k);
	//
	// Initialize the first minimzer.
	// Store canonical information in the rightest bit of can.t; The last bit = 1 ==> reverse strand
	// If canonical == True, for_mask = 0111...11 --> minimizer & for_mask = 0minimizer; rev_mask = 1000...00 --> minimizer | rev_mask = 1minimizer
	// Else for_mask = 111...11 --> minimizer & for_mask = minimizer, rev_mask = 000...00 --> minimizer | rev_mask = minimizer
	//
	Tup for_mask = TupPos::for_mask_s;	
	Tup rev_mask = TupPos::rev_mask_s;
	Tup mask = 0;
	if (Global) {
		rev_mask = (rev_mask & mask); // 0000...00 64 bits
		for_mask = ~rev_mask; // 111...11 64 bits
	}
	// cerr << "Global: " << Global << endl;
	// cerr << "for_mask: " << for_mask << endl;
	// cerr << "rev_mask: " << rev_mask << endl;

	can.t = cur.t;
	minPos = 0;
	TupPos activeMinimizer, curMinimizer;
	activeMinimizer.t = can.t;
	activeMinimizer.pos = 0;
	vector<TupPos> curTuples(w);
	curTuples[0] = activeMinimizer;

	// 
	// Find the active minimizer in this window
	//
	int nMinimizers=1;

	for (p = 1; p < w && p < seqLen-k+1 ; p++) {
		ShiftOne(seq, p+k-1, m, cur);
		// ShiftOneRC(seq, p+k-1, k, curRC);
		/*
		Tuple test, testrc;
		StoreTuple(seq->seq.s, p, k, test);
		TupleRC(test, testrc, k);
		assert(test == cur);
		assert(testrc == curRC);
		*/
		curMinimizer.pos = p;
		curMinimizer.t = (cur.t & for_mask);
		// if ((cur.t & for_mask) < (curRC.t & for_mask)) curMinimizer.t = (cur.t & for_mask);
		// else curMinimizer.t = (curRC.t | rev_mask); 
		if (curMinimizer.t < activeMinimizer.t) {  
			activeMinimizer.t = curMinimizer.t;
			activeMinimizer.pos = p;
		}	
		curTuples[p%w] = curMinimizer;
	}
	minimizers.push_back(activeMinimizer);
	// Now scan the chromosome
	minTuple.t=m.t;
	for (p = w; p < seqLen-k+1; p++) {
		// Check if past current active minimzier
		ShiftOne(seq, p+k-1, m, cur);
		// ShiftOneRC(seq, p+k-1, k, curRC);

		curMinimizer.t = (cur.t & for_mask);
		// if ((cur.t & for_mask) < (curRC.t & for_mask)) curMinimizer.t = (cur.t & for_mask);
		// else curMinimizer.t = (curRC.t | rev_mask); 
		curMinimizer.pos = p;
		curTuples[p%w] = curMinimizer;
		if (p - w >= activeMinimizer.pos) {
			activeMinimizer = curTuples[0];
			for (int j =1; j < w; j++) {
				if ((curTuples[j].t & for_mask) < (activeMinimizer.t & for_mask)) { 
					activeMinimizer = curTuples[j];
				}		
			}
			minimizers.push_back(activeMinimizer);
			nMinimizers+=1;
		}
		else {
			if ((curMinimizer.t & for_mask) < (activeMinimizer.t & for_mask)) { //TODO(Jingwen)
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

