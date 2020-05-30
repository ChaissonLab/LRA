#ifndef COMPARE_LISTS_H_
#define COMPARE_LISTS_H_
#include <algorithm>
#include "Options.h"
#include "Types.h"



template<typename tup, typename Tup> 
void CompareLists(typename vector<tup>::iterator qBegin, typename vector<tup>::iterator qEnd, 
						typename vector<tup>::iterator tBegin, typename vector<tup>::iterator tEnd, 
						vector<pair<tup, tup> > &result, Options &opts, long long int maxDiagNum = 0,
						 long long int minDiagNum = 0, bool canonical=true) {
    Tup Bi=1; Tup Ai=1;
    int nOfBits=0;
    while (Bi != 0) {
        Bi = Bi << 1;
        nOfBits++;
    }
	Tup for_mask = ~(Ai << (nOfBits-1));
	if (!canonical) {
		for_mask = ~(for_mask & 0);
	}
    cerr << "nOfBits: " << nOfBits << endl;
    cerr << "canonical: " << canonical << "  for_mask: " << for_mask << endl;

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
  std::set_intersection(qBegin, qEnd, tBegin, tEnd, back_inserter(isect));
	cout << "Matched " << isect.size() << " slowly." << endl;
#endif
	int nMatch=0;
	int iter=0;
	do {
		tup startGap, endGap;
		++iter;
		while (qs <= qe and qBegin[qs].t < (tBegin[ts].t & for_mask)) {
			qs++;
		}
		startGap.t = qBegin[qs].t - (tBegin[ts].t & for_mask);
		if (qs == qe) {
			endGap = startGap;
		}
		else {
			// Move past any entries guaranteed to not be in target
			while (qe > qs and te > ts and qBegin[qe].t > (tBegin[te-1].t & for_mask)) {
				qe--;
			}
			endGap.t = (tBegin[te-1].t & for_mask) - qBegin[qe].t;
		}
		if (startGap > endGap) {
			//
			// Find entry in t that could match qs
			//
			typename vector<tup>::iterator lb;
			lb = lower_bound(tBegin+ts, tBegin+te, qBegin[qs]);
			ts=lb-tBegin;
			if ((tBegin[ts].t & for_mask) == qBegin[qs].t) {
				GenomePos tsStart=ts;
				GenomePos tsi=ts;
				while (tsi != te and qBegin[qs].t == (tBegin[tsi].t & for_mask)) {
					tsi++;
				}
				GenomePos qsStart=qs;
				while (qs < qe and qBegin[qs+1].t == qBegin[qs].t) { qs++; }
				
				if (tsi - tsStart < opts.globalMaxFreq) {
					for(GenomePos ti=tsStart; ti != tsi; ti++) {
						for (GenomePos qi=qsStart; qi <= qs; qi++) {
							if (maxDiagNum != 0 and minDiagNum != 0) {
								long long int Diag = (long long int) tBegin[ti].pos - (long long int) qBegin[qi].pos;
								if (Diag <= maxDiagNum and Diag >= minDiagNum) result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
							}
							else result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
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
			if (tBegin + te != tEnd and (tBegin[te-1].t & for_mask) == qBegin[qe].t) {
				// pass
			} 
			else {
				ub = upper_bound(tBegin+ts, tBegin+te, qBegin[qe]);
				// *ub is > qBegin[qe]
				te = ub - tBegin;
			}
			GenomePos teStart=te, tei=te;
			while (tei > ts and (tBegin[tei-1].t & for_mask) == qBegin[qe].t) {
				tei--;
			}
			if (tei < teStart and teStart > 0) {
				GenomePos qeStart=qe;
				while (qe > qs and qBegin[qe].t == qBegin[qe-1].t) { qe--;}

				if (teStart - te < opts.globalMaxFreq) {
					for (GenomePos ti=tei; ti < teStart; ti++) {
						for (GenomePos qi=qe; qi <= qeStart; qi++) {
							if (maxDiagNum != 0 and minDiagNum != 0) {
								long long int Diag = (long long int) tBegin[ti].pos - (long long int) qBegin[qi].pos;
								if (Diag <= maxDiagNum and Diag >= minDiagNum) result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
							}
							else result.push_back(pair<tup,tup>(qBegin[qi], tBegin[ti]));
						}
					}
				}
			}
			te=tei;
		}
	} while (qs < qe and ts < te);
}

template<typename tup, typename Tup> 
void CompareLists(vector<tup> &query, vector<tup> &target, vector<pair<tup, tup> > &result, Options &opts) {
	CompareLists<tup, Tup>(query.begin(), query.end(), target.begin(), target.end(), result, opts);
}

#endif
