#ifndef COMPARE_LISTS_H_
#define COMPARE_LISTS_H_
#include <algorithm>
#include "Options.h"


template<typename tup> void CompareLists(typename vector<tup>::iterator qBegin,
																				 typename vector<tup>::iterator qEnd,
																				 typename vector<tup>::iterator tBegin,
																				 typename vector<tup>::iterator tEnd,
																				 vector<pair<tup, tup> > &result, 
																				 Options &opts) {
	int qs = 0;
	int ts = 0, te = tEnd - tBegin;
	typename vector<tup>::iterator slb;

#ifdef _TESTING_
	vector<tup> isect;
	cout << "comparing " << qEnd-qBegin << " and " << te-ts << " lists" << endl;
  std::set_intersection(qBegin, qEnd,
												tBegin, tEnd, back_inserter(isect));
		cout << "Matched " << isect.size() << " slowly." << endl;
#endif
	int nMatch=0;
	typename vector<tup>::iterator qi;
	
	for (qi = qBegin; qi != qEnd; ++qi) {
		slb = lower_bound(tBegin+ts, tEnd, *qi);
		ts=slb-tBegin;
		if (slb->t == qi->t) {
			typename vector<tup>::iterator slb0 = slb;
			while (slb != tEnd && 
						 qi->t == slb->t) {
				slb++;
			}

			if (slb - slb0 < opts.maxFreq) {
				for(slb0; slb0 != slb; ++slb0) {
					result.push_back(pair<tup,tup>(*qi, *slb0));
				}
			}
		}
	}	
}

template<typename tup> void CompareLists(vector<tup> &query, 
																				 vector<tup> &target,
																				 vector<pair<tup, tup> > &result, 
																				 Options &opts) {
	CompareLists(query.begin(), query.end(),							 
							 target.begin(), target.end(), result, opts);
}

#endif
