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

#define NUMPWL 20
#define MAXPWL 20000 
static int firstGapCeiling;
static int secondGapCeiling;
static long STOPS[NUMPWL];
static float INTER[NUMPWL];
static float SLOPE[NUMPWL];


float nroot(float x, float root) {
  return std::pow(x,1/root);
}

void InitPWL(float intercept, float scalar, float root, int gapCeiling1=1500, int gapCeiling2=1500) { // 2000; 50000
  //
  // Determine the spacing.
  //
firstGapCeiling=gapCeiling1;
secondGapCeiling=gapCeiling2;
STOPS[0] = 0;

  int width=MAXPWL/NUMPWL;
  float vals[NUMPWL];
  vals[0] = 0;
  STOPS[1]=5;
  STOPS[2]=10;
  STOPS[3]=20;
  STOPS[4]=40;
  STOPS[5]=80;
  STOPS[6]=100;
  STOPS[7]=200;
  STOPS[8]=300;
  STOPS[9]=500;
  STOPS[10]=1000;
  STOPS[11]=2000;
  STOPS[12]=3000;
  STOPS[13]=4000;
  STOPS[14]=5000;
  STOPS[15]=6000;
  STOPS[16]=7000;
  STOPS[17]=8000;
  STOPS[18]=9000;
  STOPS[19]=20000;
  /*
  for (int i=1; i < NUMPWL;i++) {
    STOPS[i] = i*width; //floor(nroot(i*width, root));
    // intercept and scalar can be added later, but why not use now to not forget
  }
  */
  for (int i=1; i < NUMPWL;i++) {
  	if (i <= 2) intercept = 0;
    vals[i] = intercept+scalar*nroot(STOPS[i], root);	
  }
  for (int i=0; i < NUMPWL-1; i++) {		
    float slope=(vals[i+1]-vals[i])/(STOPS[i+1]-STOPS[i]);
		if (STOPS[i] <= 10) {
			SLOPE[i] = 0;
			INTER[i] = 0;
		}
		else {
			SLOPE[i] = slope;
			INTER[i] = vals[i]-STOPS[i]*slope+intercept;
		}
  }
}

float PWL_w(long x, long minX=0) {
  // no gap is no penalty
	minX=2;
  long penalty;
int bound=0;
  if (x <= minX) { penalty=0;}
	else {
     bound=std::upper_bound(&STOPS[0], &STOPS[NUMPWL-1], x)-&STOPS[0];
	   penalty=SLOPE[bound-1]*x + INTER[bound-1];
     //
     // Add two steps to the end of the gap function, one for medium length gaps, and one for long
     //
     if (penalty >= firstGapCeiling && penalty <secondGapCeiling) {
        penalty= firstGapCeiling;
     }
     else if (penalty > secondGapCeiling) {
        penalty = secondGapCeiling;
     }
}


return penalty;

}
  

// // w function 
// float
// w (long int i, long int j, const std::vector<float> & LookUpTable, const Options &opts, bool &step_sdp) {  
// // step_sdp == 0 means the first sdp; step_sdp == 1 means the second sdp;
// 	long int x = labs(j - i); 
// 	if (x == 0) return 0;
// 	int a = (int) floor(x/5);
// 	if (step_sdp == 0) {
// 		if (x <= 100) {
// 			return - opts.firstcoefficient * (0.5f * x + opts.gapopen);
// 		}
// 		else if (x <= 500){
// 			return - opts.firstcoefficient * (0.375f * x + 12.5f + opts.gapopen) - 300;
// 		}
// 		else if (x <= 1000) {
// 			return - opts.firstcoefficient * (0.2f * x + 100 + opts.gapopen) - 300;
// 		}
// 		else if (x <= 5000) {
// 			return - opts.firstcoefficient * (0.175f * x + 125 + opts.gapopen) - 300;
// 		}
// 		else if (x <= 20000) {
// 			return - opts.firstcoefficient * (0.1f * x + 500 + opts.gapopen) - 300;
// 		}
// 		else if (x <= 50000) {
// 			return - opts.firstcoefficient * (0.05f * x + 1500 + opts.gapopen) - 300;
// 		}
// 		else if (x <= 100000) {
// 			return - opts.firstcoefficient * (0.04f * x + 2000 + opts.gapopen) - 300;
// 		}
// 		else {
// 			return - opts.firstcoefficient * 6000 - 300;
// 		}
// 	}
// 	else {
// 		if (x <= 100) {
// 			return - opts.secondcoefficient * (0.5f * x + opts.gapopen);
// 		}
// 		else if (x <= 500){
// 			return - opts.secondcoefficient * (0.375f * x + 12.5f + opts.gapopen) - 300;
// 		}
// 		else if (x <= 1000) {
// 			return - opts.secondcoefficient * (0.2f * x + 100 + opts.gapopen) - 300;
// 		}
// 		else if (x <= 5000) {
// 			return - opts.secondcoefficient * (0.175f * x + 125 + opts.gapopen) - 300;
// 		}
// 		else if (x <= 20000) {
// 			return - opts.secondcoefficient * (0.1f * x + 500 + opts.gapopen) - 300;
// 		}
// 		else if (x <= 50000) {
// 			return - opts.secondcoefficient * (0.05f * x + 1500 + opts.gapopen) - 300;
// 		}
// 		else if (x <= 100000) {
// 			return - opts.secondcoefficient * (0.04f * x + 2000 + opts.gapopen) - 300;
// 		}
// 		else {
// 			return - opts.secondcoefficient * 6000 - 300;
// 		}
// 	}
// }


// w function 
float
w (long int i, long int j, const std::vector<float> & LookUpTable, const Options &opts, bool &step_sdp) {  // step_sdp == 0 means the first sdp; step_sdp == 1 means the second sdp;
	long int x = labs(j - i) + 1; 
	if (x == 1) return 0;
	int a = (long) floor((x-1)/5);
	//	float exact=opts.gapextend*nroot(x,opts.root)+opts.gapopen;
	//	float pwl=PWL_w(x);
	return -PWL_w(x, opts.freeGap);
	if (step_sdp == 0) {
		if (opts.LookUpTable) {
			if (x <= 20) {
				return - x - opts.gapopen;
			}
			else if (x <= 10001){
				return - opts.firstcoefficient*LookUpTable[a] - opts.gapopen;
			}
			else if (x <= 500001) {
				return -2000 - opts.gapopen;
			}
			else if (x <= 100001){
				return -4000 - opts.gapopen;
			}
			else {
				return -6000 - opts.gapopen;
			}
		}
		else {
			if (x < 501) {
				return - opts.firstcoefficient*logf(x) - opts.gapopen;  
			}
			else if (x <= 10001) {
				return - opts.firstcoefficient*logf(x) - opts.gapopen;  
			}
			else if (x <= 100001) {
				return - opts.firstcoefficient*logf(x) - opts.gapopen;
			}
			else {
				return - opts.firstcoefficient*logf(x) - opts.gapopen;	
			}
		}
	}
	else {
		if (opts.LookUpTable) {
			if (x <= 20) {
				return -x- opts.gapopen;
			}
			else if (x <= 10001){
				return - opts.secondcoefficient*LookUpTable[a] - opts.gapopen;
			}
			else if (x <= 500001) {
				return -2000 - opts.gapopen;
			}
			else if (x <= 100001){
				return -4000 - opts.gapopen; //-800
			}
			else {
				return -6000 - opts.gapopen;
			}
		}
		else {
			if (x < 501) {
				return - opts.secondcoefficient*logf(x) - opts.gapopen;  
			}
			else if (x <= 10001) {
				return - opts.secondcoefficient*logf(x) - opts.gapopen;  
			}
			else if (x <= 100001) {
				return - opts.secondcoefficient*logf(x) - opts.gapopen;
			}
			else {
				return - opts.secondcoefficient*logf(x) - opts.gapopen;	
			}
		}		
	}
}

	// original gap penalty
	// 	if (x < 501) {
	//     	return - opts.coefficient*logf(x) - 1;   
	// 	}
	// 	else if (x <= 10001){
	// 		// check LookUpTable
	// 		// TODO(Jingwen): finish the code here
	// 		float f = std::floor((x-501)/5);
	// 		int a = (int) (f);
	// 		return - opts.coefficient*LookUpTable[a] - 1;
	// 	}
	// 	else if (x <= 100001){
	// 		if (step_sdp == 1) return -1000; //-800
	// 		else return -2000;
	// 	}
	// 	else {
	// 		if (step_sdp == 1) return -2000;
	// 		else return -4000;
	// 	}
	// }




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
FindValueInBlock (long int ForwardDiag, std::stack<LPair> & S_1, std::vector<long int> & Ei, std::vector<LPair> & Block, 
					unsigned int & i1, unsigned int & i2) {

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
FindBoundary (unsigned int first, unsigned int last, unsigned int a, unsigned int b, std::vector<long int> & Di, 
				std::vector<float> & Dv, std::vector<long int> & Ei, const std::vector<float> & LookUpTable, const Options &opts, bool &step_sdp) {
	if (b != -1) {
		unsigned int it;
		unsigned int count, step;
		count = last - first;
		while (count > 0) {
			it = first; step = count/2; it += step;
			if (Dv[a] + w(Di[a], Ei[it], LookUpTable, opts, step_sdp) > Dv[b] + w(Di[b], Ei[it], LookUpTable, opts, step_sdp)) { // if a is better than b
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
Maximization (unsigned int & now, long int & last, std::vector<long int> & Di, std::vector<long int> & Ei, std::vector<float> & Dv, 
				std::vector<long int> & Db, std::vector<std::pair<long int, long int>> & Block, std::stack<LPair> & S_1, 
				const std::vector<float> & LookUpTable, const Options &opts, bool &step_sdp) { // last and now are both index

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

			if (Dv[i] + w(Di[i], Ei[Db[i]], LookUpTable, opts, step_sdp) > Dv[l] + w(Di[l], Ei[Db[i]], LookUpTable, opts, step_sdp)) { // Di[i] is better than Di[l] at Db[i]
				
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
				while (!S_1.empty() and Dv[i] + w(Di[i], Ei[cur.second - 1], LookUpTable, opts, step_sdp) > Dv[cur.first]
								 + w(Di[cur.first], Ei[cur.second - 1], LookUpTable, opts, step_sdp)) {
					//cerr << "t " << endl;
					S_1.pop();
					prev = cur;
					cur = S_1.top();
					if (cur.second == n + 1) break;
				}
				//cerr << "prev: " << prev << endl;
				//cerr << "cur: " << cur << endl;
				unsigned int h = FindBoundary(prev.second, cur.second, i, cur.first, Di, Dv, Ei, LookUpTable, opts, step_sdp);
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
