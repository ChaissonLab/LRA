#ifndef COMPARE_LISTS_H_
#define COMPARE_LISTS_H_
#include <algorithm>
#include "Options.h"


template<typename tup> void CompareLists(typename vector<tup>::iterator qBegin, typename vector<tup>::iterator qEnd, typename vector<tup>::iterator tBegin,
											typename vector<tup>::iterator tEnd, vector<pair<tup, tup> > &result, Options &opts) {
	int qs = 0;
	int qe = qEnd-qBegin - 1;
	int ts = 0, te = tEnd - tBegin;
	typename vector<tup>::iterator slb;
	if (qBegin == qEnd or tBegin == tEnd) {
		result.clear();
		return;
	}
	
#ifdef _TESTING_
	vector<tup> isect;
	cout << "comparing " << qEnd-qBegin << " and " << te-ts << " lists" << endl;
  std::set_intersection(qBegin, qEnd,
												tBegin, tEnd, back_inserter(isect));
	cout << "Matched " << isect.size() << " slowly." << endl;
#endif
	int nMatch=0;
	int iter=0;
	do {
		tup startGap, endGap;
		++iter;
		while (qs <= qe and qBegin[qs].t < tBegin[ts].t) {
			qs++;
		}
		startGap.t = qBegin[qs].t - tBegin[ts].t;
		if (qs == qe) {
			endGap = startGap;
		}
		else {
			// Move past any entries guaranteed to not be in target
			while (qe > qs and te > ts and qBegin[qe].t > tBegin[te-1].t) {
				qe--;
			}
			endGap.t = tBegin[te-1].t - qBegin[qe].t;
		}
		if (startGap > endGap) {
			//
			// Find entry in t that could match qs
			//
			typename vector<tup>::iterator lb;
			lb = lower_bound(tBegin+ts, tBegin+te, qBegin[qs]);
			ts=lb-tBegin;
			if (tBegin[ts].t == qBegin[qs].t) {
				GenomePos tsStart=ts;
				GenomePos tsi=ts;
				while (tsi != te && 
							 qBegin[qs].t == tBegin[tsi].t) {
					tsi++;
				}
				GenomePos qsStart=qs;
				while (qs < qe and qBegin[qs+1].t == qBegin[qs].t) { qs++; }
				
				if (tsi - tsStart < opts.globalMaxFreq) {
					for(GenomePos ti=tsStart; ti != tsi; ti++) {
						for (GenomePos qi=qsStart; qi <= qs; qi++) {
							result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
						}
					}
				}
			}
		}
		else {
			//
			// End gap is greater than start gap, search from the other direction
			//
			typename vector<tup>::iterator ub;
			assert(te > ts);
			if (tBegin+ te != tEnd and tBegin[te-1].t == qBegin[qe].t) {
				// pass
			} 
			else {
				ub = upper_bound(tBegin+ts, tBegin+te, qBegin[qe]);
				// *ub is > qBegin[qe]
				te = ub - tBegin;
			}
			GenomePos teStart=te, tei=te;
			while (tei > ts and tBegin[tei-1].t == qBegin[qe].t) {
				tei--;
			}
			if (tei < teStart and teStart > 0) {
				GenomePos qeStart=qe;
				while (qe > qs and qBegin[qe].t == qBegin[qe-1].t) { qe--;}

				if (teStart - te < opts.globalMaxFreq) {
					for (GenomePos ti=tei; ti < teStart; ti++) {
						for (GenomePos qi=qe; qi <= qeStart; qi++) {
							result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
						}
					}
				}
			}
			te=tei;
		}
	} while (qs < qe and ts < te);
}

template<typename tup> void CompareLists(vector<tup> &query, vector<tup> &target, vector<pair<tup, tup> > &result, Options &opts) {
	CompareLists(query.begin(), query.end(),							 
							 target.begin(), target.end(), result, opts);
}

#endif
