#ifndef NAIVE_D_P
#define NAIVE_D_P


#include <iostream> //std::cout  
#include <cmath>        // std::labs
#include <cstdio>    // std::FILE std::perror
#include <set>

#include "seqan/basic.h" // Triple
#include "MergeSplit.h" // beginPosition
#include "seqan/seeds.h" 
#include "IndexedSeed.h"
/*
typedef seqan::Seed<IndexedSeed> IndSeed;
typedef seqan::SeedSet<IndSeed> IndSeedSet;
*/



// Gap cost function: log 
int64_t
GapCost (unsigned int i, unsigned int j, unsigned int i_prime, unsigned int j_prime) { // end_x, end_y, start_x, start_y (x cordinate is read, y cordinate is genome)
	// some function about j-i - (j_prime - i_prime)
	int64_t ii = (int64_t) i;
    int64_t jj = (int64_t) j;
    int64_t ii_prime = (int64_t) i_prime;
    int64_t jj_prime = (int64_t) j_prime;

	int64_t t = (jj - ii) - (jj_prime - ii_prime);
   	double a = floor(8*log(abs(t) + 1) + 2);
    int64_t b = (int64_t)a;
    return b;  
}  

//-------------debug
int64_t
Gaplength (unsigned int i, unsigned int j, unsigned int i_prime, unsigned int j_prime) { // end_x, end_y, start_x, start_y (x cordinate is read, y cordinate is genome)
	// some function about j-i - (j_prime - i_prime)
	int64_t ii = (int64_t) i;
    int64_t jj = (int64_t) j;
    int64_t ii_prime = (int64_t) i_prime;
    int64_t jj_prime = (int64_t) j_prime;

	int64_t t = abs((jj - ii) - (jj_prime - ii_prime));
 
    return t;  
}  


/*
int64_t
GapCost (unsigned int i, unsigned int j, unsigned int i_prime, unsigned int j_prime) { // end_x, end_y, start_x, start_y (x cordinate is read, y cordinate is genome)
	// some function about j-i - (j_prime - i_prime)
	int64_t ii = (int64_t) i;
    int64_t jj = (int64_t) j;
    int64_t ii_prime = (int64_t) i_prime;
    int64_t jj_prime = (int64_t) j_prime;

	int64_t t = (jj - ii) - (jj_prime - ii_prime);
   	double a = floor(cbrt(abs(t) + 1));
    int64_t b = (int64_t)a;
    return b;  
}  
*/
/*
int64_t
GapCost (unsigned int i, unsigned int j, unsigned int i_prime, unsigned int j_prime) { // end_x, end_y, start_x, start_y (x cordinate is read, y cordinate is genome)
    // some function about j-i - (j_prime - i_prime) !!!!!!!!!!!!!!!!!!!
    int64_t ii = (int64_t) i;
    int64_t jj = (int64_t) j;
    int64_t ii_prime = (int64_t) i_prime;
    int64_t jj_prime = (int64_t) j_prime;
	int64_t t = (jj - ii) - (jj_prime - ii_prime);

    float gap_score = 2; // gap openning penalty
    for (unsigned y = 0; y < t; y++) {
			gap_score = gap_score + max(1.00, 2.00 - 0.15*y);       	
    }
    int64_t b = (int64_t)gap_score;
    return b;  
}  

*/


template<typename TSeedSet, typename TSeed>
void NaiveDP (TSeedSet &seedSet, seqan::String<TSeed> &chain) {
	seqan::String<TSeed> seeds;
	seqan::resize(seeds, seqan::length(seedSet));
	std::copy(seedSet._seeds.begin(), seedSet._seeds.end(), seqan::begin(seeds, seqan::Standard()));
	//
	//std::cout << "length(seeds): " << length(seeds) << std::endl;
	//---------------------------------------------------------------------------------------------
	// Step 1: generate the sorted list of interval points
	// --------------------------------------------------------------------------------------------

	typedef seqan::Triple<GenomePos, bool, unsigned> TIntervalPoint;
	typedef seqan::String<TIntervalPoint> TIntervalPoints; 
	typedef typename seqan::Iterator<TIntervalPoints, seqan::Standard>::Type TIntervalPointsIterator;
	TIntervalPoints intervalPoints; //intervalPoints contains all the start/end points of seeds
	std::map<unsigned, int64_t> qualityOfChainEndingIn;
	std::map<unsigned, unsigned> predecessor;

	for (unsigned i = 0; i < seqan::length(seeds); ++i) {

		qualityOfChainEndingIn[i] = seqan::seedSize(seeds[i]);
		predecessor[i] = std::numeric_limits<unsigned>::max();
		seqan::appendValue(intervalPoints, TIntervalPoint(beginPositionH(seeds[i]), true, i));
		seqan::appendValue(intervalPoints, TIntervalPoint(endPositionH(seeds[i]), false, i));
	}
	std::sort(seqan::begin(intervalPoints, seqan::Standard()), seqan::end(intervalPoints, seqan::Standard())); // end goes before start if their positions are the same

//debug code 
	//cout << "length(seeds): " << length(seeds) << endl;
	//cout << "length(intervalPoints): " << length(intervalPoints) << endl;
/*
	for (TIntervalPointsIterator it = seqan::begin(intervalPoints, seqan::Standard()); it != seqan::end(intervalPoints, seqan::Standard()); ++it) {
		if (it->i2 == true) {
			cout << "*it-true:  " << *it << endl;
		}
		else {
			cout << "*it--false:  " << *it << endl;
		}
	}

*/
	// ---------------------------------------------------------------------------------------
	// Step 2: bulid the chain
	// ----------------------------------------------------------------------------------------
	// build a list of "intermediate solutions"
	// Each solution is represented by the triple (value of best chain so far, endPosition in dim(Genome), last seed of the chain)

	typedef seqan::Triple <int64_t, GenomePos, unsigned> TIntermediateSolution;  
	typedef std::multiset<TIntermediateSolution> TIntermediateSolutions; // Elements in multiset are in ascending order
	typedef typename TIntermediateSolutions::iterator TIntermediateSolutionsIterator;

	// For all interval points.....
	TIntermediateSolutions intermediateSolutions;


	for (TIntervalPointsIterator it_k = seqan::begin(intervalPoints); it_k != seqan::end(intervalPoints); ++it_k) {
		TSeed const & seed_k = seeds[it_k->i3];


	
		if (it_k->i2) { // It's a begin point
			// Find the seed j so that seed j's Genome.cordinate is less or equal to the beginPositionV of seed_k
			
			/*
			TIntermediateSolution referenceSolution(beginPositionV(seed_k), std::numeric_limits<int64_t>::max(), std::numeric_limits<unsigned>::max());
			TIntermediateSolutionsIterator it_q = intermediateSolutions.upper_bound(referenceSolution); // the beginPositionV of it_q is larger than the beginPositionV

			// STl gives us upper_bound which returns a pointer to the first one that compares greater than the reference one. 
			// Special case: If intermediateSolutions is empty or there is no chain that ends before seed_k begins
			if (intermediateSolutions.empty() || it_q == intermediateSolutions.begin()) {
				continue;
			}
			*/

			if (intermediateSolutions.empty()) {
				continue;
			}





			int64_t quality = qualityOfChainEndingIn[it_k->i3]; // quality stores the current maximum
			for (TIntermediateSolutionsIterator it_j = intermediateSolutions.begin(); it_j != intermediateSolutions.end(); ++it_j) { // it_j->i1 <= beginPositionV(seed_k)
				//cout << "endPositionH(seeds[it_j->i3]:   " << endPositionH(seeds[it_j->i3]) << " " << endPositionV(seeds[it_j->i3]) << " "
				//	<< beginPositionH(seed_k) << " " << beginPositionV(seed_k) << endl;


				if (beginPositionV(seed_k) >= it_j->i2 && quality <= qualityOfChainEndingIn[it_k->i3] + it_j->i1 -
					GapCost(endPositionH(seeds[it_j->i3]), endPositionV(seeds[it_j->i3]), beginPositionH(seed_k), beginPositionV(seed_k))) { // Jingwen: Is it "<=" or "<"

					quality = qualityOfChainEndingIn[it_k->i3] + it_j->i1 -
					GapCost(endPositionH(seeds[it_j->i3]), endPositionV(seeds[it_j->i3]), beginPositionH(seed_k), beginPositionV(seed_k));
					
					predecessor[it_k->i3] = it_j->i3;

				}
			}
			qualityOfChainEndingIn[it_k->i3] = quality;


		}
		else { // It's an end point
			TIntermediateSolution intermediate_k(qualityOfChainEndingIn[it_k->i3], endPositionV(seeds[it_k->i3]), it_k->i3);
			intermediateSolutions.insert(intermediate_k);
		}
	}

	// -------------------------------------------------------------------------------------------
	// Step 3: Write out the resulting chain
	// -------------------------------------------------------------------------------------------

	clear(chain);
    unsigned next = intermediateSolutions.rbegin()->i3;
    while (next != std::numeric_limits<unsigned>::max())
    {
        appendValue(chain, seeds[next]);
        next = predecessor[next];
    }
    reverse(chain);


/*
    //Debug-----print intermediateSolutions
    cout << "intermediateSolutions.size():  " << intermediateSolutions.size() << endl;
    for (TIntermediateSolutionsIterator it = intermediateSolutions.begin(); it != intermediateSolutions.end(); ++it) {
    	cout << "*it:    " << *it << endl;
    }
*/
}

#endif
