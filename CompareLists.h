#ifndef COMPARE_LISTS_H_
#define COMPARE_LISTS_H_
#include <algorithm>

template<typename tup> void CompareLists(vector<tup> &query, 
																				 vector<tup> &target,
																				 vector<pair<tup, tup> > &result) {
	int qs = 0;
	int ts = 0, te = target.size();
	typename vector<tup>::iterator slb;

	//
	// Long form, testing.
	//

	vector<tup> isect;
	/*
  std::set_intersection(query.begin(), query.end(),
	target.begin(), target.end(), back_inserter(isect));
		cout << "Matched " << isect.size() << " slowly." << endl;
	*/
	int nMatch=0;
	while (qs < query.size() && ts < te ) {
		slb = lower_bound(target.begin()+ts, target.end(), query[qs]);
		if (slb->tuple == query[qs].tuple) {
			ts=slb-target.begin();
			while (ts < target.size() && 
						 query[qs].tuple == target[ts].tuple) {
				result.push_back(pair<tup,tup>(query[qs], target[ts]));
				ts++;
			}
		}
		
		qs+=1;
	}	
}

#endif
