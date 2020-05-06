#ifndef SUB_ROUTING_H_
#define SUB_ROUTING_H_

// This program implements O(mlogn + n) Minimization which is slightly differernt from the algorithm in Giancarlo paper
// we have D array which contains points in part A of some subproblems, which E array containing points in part B.


#include <iostream>	// std::cout 
#include <ostream>	// std::endl
#include <algorithm>	// std::lower_bound std::sort
#include <vector>	// std::vector
#include <stack>          // std::stack
#include <utility>      // std::pair
#include <string>	// std::string
#include <fstream>      // std::fstream
#include <tuple>        // std::tuple, std::get
#include <cmath>	// std::log std::pow
#include <numeric>	// std::iota
#include <limits>
#include <iterator>
#include <cstdlib> // std::labs
#include "Types.h"
#include "Sorting.h"

using std::cerr;
using std::cout;
using std::endl;


// w function
float
w (long int i, long int j, const std::vector<float> & LookUpTable, Options &opts, bool step) { // step == 1 means the second sdp; step == 0 means the first sdp
	long int x = labs(j-i) + 1;

	if (opts.LookUpTable) {
		if (x < 501) {
	    	return - opts.coefficient*std::log(x) - 1;   
		}
		else if (x <= 10001){
			// check LookUpTable
			// TODO(Jingwen): finish the code here
			float f = std::floor((x-501)/5);
			int a = static_cast<int> (f);
			return - opts.coefficient*LookUpTable[a] - 1;
		}
		else {
			return -1000; //-800
			//if (step == 1) {return -10000;}
			//else {return -800;}
		}
	}
	else {
		if (x < 501) {
			return - opts.coefficient*std::log(x) - 1;  
		}
		else if (x <= 10001) {
			return - opts.coefficient*std::log(x) - 1;  
		}
		else {
			if (step == 1) {return - 10000;}
			else {return - opts.coefficient*std::log(x) - 1;}
		}
	}
}  



typedef std::pair<unsigned int, unsigned int> Pair;
typedef std::pair<long int, long int> LPair;


// Find the first LPair s in [first, last) with s.second > val 
std::vector<LPair>::iterator
UPPERbound (std::vector<LPair>::iterator first, std::vector<LPair>::iterator last, unsigned int val) {
	
	std::vector<LPair>::iterator it;
	unsigned int count, step;
	count = std::distance(first, last);
	while (count > 0) {
		it = first; step = count/2; std::advance(it, step);
		if (val >= it->second) {
			first = ++it;
			count -= step + 1;
		}
		else count = step;
	}
	return first;
}



// TODO(Jingwen): Change this to first get Ev for all the points. Only use "lower_bound" to retrieve the index
void
FindValueInBlock (long int ForwardDiag, std::stack<LPair> & S_1, std::vector<long int> & Ei, std::vector<LPair> & Block, unsigned int & i1, unsigned int & i2) {

	if (i1 >= Block.back().second and i1 < S_1.top().second) {
		i2 = S_1.top().first;
	}
	else {
		std::vector<LPair>::iterator it2 = UPPERbound(Block.begin(), Block.end(), i1); // Find the best candidate index for point Ei[i1]
		i2 = it2->first;
	}
}




// Using Binary search to find the first index in [first, last) that a is worse than b
unsigned int
FindBoundary (unsigned int first, unsigned int last, unsigned int a, unsigned int b, std::vector<long int> & Di, std::vector<float> & Dv, std::vector<long int> & Ei, 
					const std::vector<float> & LookUpTable, Options &opts, bool step) {

	if (b != -1) {
		unsigned int it;
		unsigned int count, step;
		count = last - first;
		while (count > 0) {
			it = first; step = count/2; it += step;
			if (Dv[a] + w(Di[a], Ei[it], LookUpTable, opts, step) > Dv[b] + w(Di[b], Ei[it], LookUpTable, opts, step)) { // if a is better than b
				first = ++it;
				count -= step + 1;
			}
			else count = step;
		}	
	}
	else {
		first = Ei.size();
	}

	return first;
}


void
Maximization (unsigned int & now, long int & last, std::vector<long int> & Di, std::vector<long int> & Ei, std::vector<float> & Dv, std::vector<long int> & Db, 
					std::vector<std::pair<long int, long int>> & Block, std::stack<LPair> & S_1, const std::vector<float> & LookUpTable, Options &opts, bool step) { // last and now are both index

 	unsigned int m = Di.size();
 	unsigned int n = Ei.size();

	//get the block from this function

	std::vector<long int> E(n);
	std::iota(E.begin(), E.end(), 0);

	//cerr << "hallelujah\n";
	//cerr << "last: " << last << " now: " << now << "\n";

	for (unsigned int i = last + 1; i <= now; ++i) {

		//cerr << "i: " << i << endl;
		
		LPair v;

		if (Db[i] == -1) { // means there is no lower bound in array Ei for Di[i]
			//cerr << "there is no lower bound in array Ei for Di[i]\n";
			break;
		}
		else {
			if (S_1.top().second == n + 1) { // for trivial case where S_1 only has the dummy pair
				LPair dummy_pair = std::make_pair(-1, Db[i]); // This dummy pair helps to sets the start of the block to Db[i]
				Block.push_back(dummy_pair);
				LPair z = std::make_pair(i, n);
				S_1.push(z);
				//cerr << "the trivial case happens.  push dummy_pair (-1, " << Db[i] << ") to Block and push pair (" << i << ", " << n << ") to S_1\n";
 			}

			while (Db[i] >= S_1.top().second) { // trivial case
				//cerr << "while" << endl;
				v = S_1.top();
				Block.push_back(v);
				S_1.pop();
				//cerr << "while  push v in Block and pop v from S_1.  " << "v: (" << v.first << ", " << v.second << ")" << endl;
				//cerr << "Block: " << Block << ",  S_1: " << S_1 << endl;
			}



			// Update the blocks 
			long int l = S_1.top().first; 

			if (Dv[i] + w(Di[i], Ei[Db[i]], LookUpTable, opts, step) > Dv[l] + w(Di[l], Ei[Db[i]], LookUpTable, opts, step)) { // Di[i] is better than Di[l] at Db[i]
				
				//cerr << "Di[i] is better than Di[l] at Db[i]\n";

				if (Db[i] < S_1.top().second and !Block.empty() and Db[i] > (Block.back()).second) {
					LPair b = std::make_pair(S_1.top().first, Db[i]);
					Block.push_back(b);
					//cerr << "the non-trivial case happens. Push pair (" << S_1.top().first << ", " << Db[i] << ") to Block." << endl;
					//cerr << "S_1: " << S_1 << ", Block: " << Block << endl;
				}

				LPair cur = S_1.top();
				LPair prev = S_1.top();
				//cerr << "S_1: " << S_1 << endl;
				while (!S_1.empty() and Dv[i] + w(Di[i], Ei[cur.second - 1], LookUpTable, opts, step) > Dv[cur.first] + w(Di[cur.first], Ei[cur.second - 1], LookUpTable, opts, step)) {
					//cerr << "t " << endl;
					S_1.pop();
					prev = cur;
					cur = S_1.top();
					if (cur.second == n + 1) break;
				}
				//cerr << "prev: " << prev << endl;
				//cerr << "cur: " << cur << endl;
				unsigned int h = FindBoundary(prev.second, cur.second, i, cur.first, Di, Dv, Ei, LookUpTable, opts, step);
				//cerr << "h: " << h << endl;
				LPair e =  std::make_pair(i, h);
				S_1.push(e);
				//cerr << "push pair " << e << "to the S_1" << endl;
			}	

		}

	}


	unsigned int v = S_1.size();
	if (now == m - 1) { // array Di is already filled up
		// push all the elements in S_1 to Block
		while (S_1.top().second != n + 1) {
			LPair s = S_1.top();
			Block.push_back(s);	
			S_1.pop();
		}
	}
	else {
		while (Db[now + 1] >= S_1.top().second) {
			LPair s = S_1.top();
			Block.push_back(s);	
			S_1.pop();			
		}
	}

	last = now;
	//cerr << "Block: " << Block << endl;
	//cerr << "S_1: " << S_1 << endl;

}


#endif