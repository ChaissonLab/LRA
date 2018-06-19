#ifndef COMPARE_LISTS_H_
#define COMPARE_LISTS_H_
#include <algorithm>
#include "Options.h"

template<typename tup> void CompareLists(vector<tup> &query, 
																				 vector<tup> &target,
																				 vector<pair<tup, tup> > &result, 
																				 Options &opts) {
	int qs = 0;
	int ts = 0, te = target.size();
	typename vector<tup>::iterator slb;

	//
	// Long form, testing.
	//

	vector<tup> isect;
#ifdef _TESTING_
	cout << "comparing " << query.size() << " and " << target.size() << " lists" << endl;
  std::set_intersection(query.begin(), query.end(),
	target.begin(), target.end(), back_inserter(isect));
		cout << "Matched " << isect.size() << " slowly." << endl;
#endif
	int nMatch=0;
	while (qs < query.size() && ts < te ) {
		slb = lower_bound(target.begin()+ts, target.end(), query[qs]);
		ts=slb-target.begin();
		if (slb->t == query[qs].t) {
			int ts0 = ts;
			while (ts < target.size() && 
						 query[qs].t == target[ts].t) {
				ts++;
			}

			if (ts - ts0 < opts.maxFreq) {
				for(int i =ts0; i < ts; i++) {
					result.push_back(pair<tup,tup>(query[qs], target[i]));
				}
			}
		}
		
		qs+=1;
	}	
}

#endif
