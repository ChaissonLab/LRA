#ifndef GLOBAL_CHAIN_H_
#define GLOBAL_CHAIN_H_
using namespace std;
#include "TupleOps.h"
#include <algorithm>
#include "PrioritySearchTree.h"
#include <assert.h>
#include "Fragment.h"
#include <iostream>
class Endpoint {
 public:
	int x;
	int y;
	int fragment;

	int score;
	int chainPrev;
	typedef unsigned int KeyType;

	enum WhichEnd {Start, End};
	WhichEnd side;
 Endpoint(int _x, int _y, int _f, WhichEnd _s) : x(_x), y(_y), fragment(_f), side(_s) { chainPrev=-1;}
	int GetX() const {
		return x;
	}
	int GetY() const {
		return y;
	}
	void SetScore(int s) {
		score = s;
	}
	WhichEnd GetSide() const {
		return side;
	}
	bool operator()(const Endpoint &rhs) const {
		if (x != rhs.x) {
			return x <rhs.x;
		}
		else {
			return y < rhs.y;
		}
	}
	Endpoint(){}
	class LessThan {
	public:
		int operator()(const Endpoint &lhs, const Endpoint &rhs) const {
			if (lhs.x != rhs.x) {
				return lhs.x < rhs.x;
			}
			else {
				return lhs.y < rhs.y;
			}
		}
	};

	int GetKey() const { 
		return y;
	}

};



using namespace std;
template<typename T_Fragment,typename T_Endpoint>
	void FragmentSetToEndpoints(vector<T_Fragment> &fragments, vector<T_Endpoint> &endpoints) {
	
	endpoints.resize(fragments.size()*2);
	
	int i;
	int ep = 0;
	for (i = 0; i < fragments.size(); i++) {
		endpoints[ep].x = fragments[i].xl;
		endpoints[ep].y = fragments[i].yl;
		endpoints[ep].side = T_Endpoint::Start;
		endpoints[ep].fragment = i;
		ep++;

		endpoints[ep].x = fragments[i].xh;
		endpoints[ep].y = fragments[i].yh;
		endpoints[ep].side = T_Endpoint::End;
		endpoints[ep].fragment = i;
		ep++;
	}
}

template<typename T_Fragment, typename T_Endpoint>
	int GlobalChain( vector<T_Fragment> &fragments,
									 vector<int> &optFragmentChainIndices,
									 vector<T_Endpoint> &endpoints) {
	

	//
	// Initialize the fragment score to be the length of each fragment.
	//
	if (fragments.size() == 0) {
		return 0;
	}

	//
	// Add the start/end points of each fragment. This allows separate scoring
	// of start points and activation of endpoints.
	//
				
	
	FragmentSetToEndpoints<T_Fragment, T_Endpoint>(fragments, endpoints);

	//
	// The Starting points of all the fragmements are in order, 
	// but not necessarily all of the end endpoints, so
	// the list must be resorted.
	//
	std::sort(endpoints.begin(), endpoints.end(), typename T_Endpoint::LessThan());
	
	PrioritySearchTree<T_Endpoint> pst;

	pst.CreateTree(endpoints);
	//	pst.Print();
	unsigned int p;
	unsigned int maxScoringEndpoint = 0;
	bool maxScoringEndpointFound = false;
	for (p = 0; p < endpoints.size(); p++) {
		int x = endpoints[p].x;
		int y = endpoints[p].y;
		if (endpoints[p].GetSide() == T_Endpoint::Start) {
			int maxPointIndex=0;
			if (pst.FindIndexOfMaxPoint(endpoints, endpoints[p].y, maxPointIndex)) {
				assert(endpoints[maxPointIndex].fragment != endpoints[p].fragment);
				int fPrev = endpoints[maxPointIndex].fragment;
				fragments[endpoints[p].fragment].prev = fPrev;
				int score = fragments[endpoints[maxPointIndex].fragment].score + fragments[endpoints[p].fragment].score;				
				/*
				cerr << "Score at " << endpoints[p].x << "\t" << endpoints[p].y << "\t" << fragments[fPrev].xl  << "\t" 
						 << fragments[fPrev].yl  << "\t" 
						 << fragments[fPrev].xh  << "\t" 
						 << fragments[fPrev].yh  << "\t" 
						 << fragments[fPrev].score  << "\tscore:\t" << score << endl;
				*/
				fragments[endpoints[p].fragment].score = score;
				//				pst.Print();
			}
			else {
				fragments[endpoints[p].fragment].prev = -1;
			}
		}	else {
			assert(endpoints[p].GetSide() == T_Endpoint::End);
			// 
			// The score of the fragment should be already set.  So simply activate
			// it here (make the point be visible in a search).
			//
			endpoints[p].score =  fragments[endpoints[p].fragment].score;
			pst.Activate(endpoints, p);
			if (maxScoringEndpointFound == false or
					fragments[endpoints[maxScoringEndpoint].fragment].score < fragments[endpoints[p].fragment].score) {
				maxScoringEndpoint = p;
				maxScoringEndpointFound = true;
			}
		}
	}
	
	// 
	// Now compute the chain of optimum fragments
	//
	T_Fragment *optFragmentPtr;
	if (maxScoringEndpointFound == false) {
		// 
		// Null case, no endpoints have been processed.
		//
		return 0;
	}

	int prev = endpoints[maxScoringEndpoint].fragment;
	unsigned int numIter = 0;
	while (prev != -1 ) {
		optFragmentChainIndices.push_back(prev);

		prev = fragments[prev].prev;
		// 
		// Do a sanity check to make sure this loop is finite -- the optimal
		// fragment chain should never contain more fragments than what are
		// input.
		//
		assert(numIter < fragments.size());
		++numIter;
	}
  reverse(optFragmentChainIndices.begin(), optFragmentChainIndices.end());
	return optFragmentChainIndices.size();

}

#endif
