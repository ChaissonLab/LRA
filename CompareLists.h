#ifndef COMPARE_LISTS_H_
#define COMPARE_LISTS_H_
#include <algorithm>
#include "Options.h"
#include "Types.h"

// Compare minimizers from reference and reads without frequency
template<typename tup, typename Tup> 
void CompareLists(typename vector<tup>::iterator qBegin, typename vector<tup>::iterator qEnd, typename vector<tup>::iterator tBegin, typename vector<tup>::iterator tEnd, 
						vector<pair<tup, tup>> &result, const Options &opts, bool Global, int64_t maxDiagNum = 0, int64_t minDiagNum = 0, bool canonical=true) {
	//
	// If canonical == True, for_mask = 0111...11 --> minimizer & for_mask = 0minimizer.
	// Else, for_mask = 111...11 --> minimizer & for_mask = minimizer
	// For LocalTuple, for_mask = 111...11 always.
	//
	Tup for_mask = tup::for_mask_s;
	if (!canonical and Global) {
		for_mask = ~(for_mask & 0);
	}
    // cerr << "canonical: " << canonical << "  for_mask: " << for_mask << endl;
	int nAdded=0;
	int origSize=result.size();
	long qs = 0, qe = qEnd - qBegin - 1;
	long ts = 0, te = tEnd - tBegin;
	typename vector<tup>::iterator slb;

	if (qBegin == qEnd or tBegin == tEnd) {
		result.clear();
		return;
	}
	long maxFreq;
	if ( Global ) { maxFreq = opts.globalMaxFreq;} else { maxFreq = opts.localMaxFreq;}
	
#ifdef _TESTING_
	vector<tup> isect;
	cout << "comparing " << qEnd-qBegin << " and " << te-ts << " lists" << endl;
  std::set_intersection(qBegin, qEnd, tBegin, tEnd, back_inserter(isect));
	cout << "Matched " << isect.size() << " slowly." << endl;
#endif

	long nMatch=0;
	long iter=0;
	do {
		tup startGap, endGap;
		++iter;
		// ts is at the beginning of the matches in t. Anything before that in qBegin can't match, sl skip past that
		while (qs <= qe and (qBegin[qs].t & for_mask) < (tBegin[ts].t & for_mask)) {
			qs++;
		}
		// no possible matches, exit.
		if (qs >= qe) {
		  return;
		}
		// startGap.t is how far below the tuple in t is the tuple in q
		if (qs < qe) {
		  startGap.t = (qBegin[qs].t & for_mask) - (tBegin[ts].t & for_mask);
		}
		if (qs == qe) {
			endGap = startGap;
		}
		else {
			// Move past any entries guaranteed to not be in target
			while (qe > qs and te > ts and (qBegin[qe].t & for_mask) > (tBegin[te-1].t & for_mask)) {
				qe--;
			}
			// End gap  is how much less the last element in qBegin is than TBegin
			endGap.t = (tBegin[te-1].t & for_mask) - (qBegin[qe].t & for_mask);
		}
	if (startGap.t == 0 or startGap > endGap) {
			//
			// Find entry in t that could match qs
			//
		    typename vector<tup>::iterator lb;
			long tsOrig=ts;
			long qsOrig=qs;
			lb = lower_bound(tBegin + ts, tBegin + te, qBegin[qs]);
			ts = lb - tBegin;
			if ((tBegin[ts].t & for_mask) == (qBegin[qs].t & for_mask)) {
				GenomePos tsStart = ts;
				GenomePos tsi = ts;
				while (tsi != te and (qBegin[qs].t & for_mask) == (tBegin[tsi].t & for_mask)) { tsi++; }
				GenomePos qsStart = qs;
				while (qs < qe and (qBegin[qs+1].t & for_mask) == (qBegin[qs].t & for_mask)) { qs++; }				
				for (GenomePos ti = tsStart; ti != tsi; ti++) {
					if (qs - qsStart < maxFreq) {
						for (GenomePos qi = qsStart; qi <= qs; qi++) {
				      		if (maxDiagNum != 0 and minDiagNum != 0) {
								int64_t Diag = (int64_t) tBegin[ti].pos - (int64_t) qBegin[qi].pos;
								if (Diag <= maxDiagNum and Diag >= minDiagNum) {
									result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
								}
							}
				     		 else {								
								result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
				      		}
				    	}
					}
				}
			}

			while (ts < te and tBegin[ts].t == tBegin[tsOrig].t) { ts++; }
			while (qs < qe and qBegin[qs].t == qBegin[qsOrig].t)  { qs++;}
		}
		else {
			//
			// End gap is greater than start gap, search from the other direction
			//
			typename vector<tup>::iterator ub;
			assert(te > ts);
			if (tBegin + te != tEnd and (tBegin[te-1].t & for_mask) == (qBegin[qe].t & for_mask)) {
				// pass
			} 
			else {
				ub = upper_bound(tBegin+ts, tBegin+te, qBegin[qe]);
				// *ub is > qBegin[qe]
				te = ub - tBegin;
			}
			GenomePos teStart=te, tei=te;
			while (tei > ts and (tBegin[tei-1].t & for_mask) == (qBegin[qe].t & for_mask)) {
				tei--;
			}
			if (tei < teStart and teStart > 0) {
				GenomePos qeStart=qe;
				while (qe > qs and (qBegin[qe].t & for_mask) == (qBegin[qe-1].t & for_mask)) { qe--;}
				for (GenomePos ti = tei; ti < teStart; ti++) {
					if (qeStart - qe < maxFreq) {
					    for (GenomePos qi = qe; qi <= qeStart; qi++) {
					      	if (maxDiagNum != 0 and minDiagNum != 0) {
								int64_t Diag = (int64_t) tBegin[ti].pos - (int64_t) qBegin[qi].pos;
								if (Diag <= maxDiagNum and Diag >= minDiagNum) {
									result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
								}
					      	}
					      	else {
								result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
					     	}
					    }
					}
				}
			}
			te=tei;
		}
	} while (qs < qe and ts < te);
	nAdded=result.size() - origSize;

}

template<typename tup, typename Tup> 
void CompareLists(vector<tup> &query, vector<tup> &target, vector<pair<tup, tup> > &result, const Options &opts, bool Global) {
	CompareLists<tup, Tup>(query.begin(), query.end(), target.begin(), target.end(), result, opts, Global);
}

#endif
