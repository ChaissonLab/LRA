// generating k-mers of fixed length + merging k-mers into longer fragments + preparing for spase dp
#ifndef MERGE_SPLIT_H_
#define MERGE_SPLIT_H_

#include <iostream> //std::cout 
#include <fstream>   
#include <cstdlib>   // std::labs, std::EXIT FAILURE, std::EXIT SUCCESS
#include <string> 
#include <cmath>        // std::labs
#include <cstdio>    // std::FILE std::perror
#include <vector>
#include <set>
#include "seqan/seeds.h"
#include "IndexedSeed.h"
#include "Clustering.h"

using namespace std;


typedef seqan::Seed<IndexedSeed> IndSeed;
typedef seqan::SeedSet<IndSeed> IndSeedSet;
typedef seqan::Iterator<IndSeedSet>::Type TIterator;


template <typename TConfig>
typename seqan::Position<seqan::Seed<IndexedSeed, TConfig> >::Type
	beginPositionH(seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._beginPositionH;
}


template <typename TConfig>
typename seqan::Position<seqan::Seed<IndexedSeed, TConfig> >::Type
beginPositionV(seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._beginPositionV;
}


template <typename TConfig>
typename seqan::Position<seqan::Seed<IndexedSeed, TConfig> >::Type
endPositionH(seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._endPositionH;
}


template <typename TConfig>
typename seqan::Position<seqan::Seed<IndexedSeed, TConfig> >::Type
endPositionV(seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._endPositionV;
}


// Debugging code
template <typename TStream, typename TConfig>
inline TStream &
operator<<(TStream & stream, seqan::Seed<IndexedSeed, TConfig> const & seed)
{
    return stream << "Seed<IndexedSeed, TConfig>(" << beginPositionH(seed)
                  << ", " << beginPositionV(seed) << ", "
                  << endPositionH(seed) << ", "
                  << endPositionV(seed) << ", lower diag = "
                  << lowerDiagonal(seed) << ", upper diag = "
                  << upperDiagonal(seed) << ")";
}


int64_t 
Diagonalnum (GenomePair & t) {
	return (int64_t)t.first.pos - (int64_t)t.second.pos;
}





// Debug code
void 
SaveOriginalSeed (Cluster &rCr, FILE* fh, int k) {
	if (fh == NULL) {
		perror("Eorror opening file: ");
	}
	else {
		for (vector<GenomePair>::iterator it = rCr.matches.begin(); it != rCr.matches.end(); ++it) {
			fprintf(fh, "%u %u %u %u\n", (*it).first.pos, (*it).second.pos, (*it).first.pos + k - 1, (*it).second.pos + k - 1);
		}
	}
}

void 
SaveLongSeed (vector<vector<GenomePair*>> &v, FILE* fp, int k) {
	if (fp == NULL) {
		perror("Error opening file: ");
	}
	else {
		for (unsigned int tt = 0; tt < v.size(); ++tt) {
			if (! v[tt].empty()) { // v[tt] is not empty
				fprintf(fp, "%u %u %u %u\n", (**(v[tt].begin())).first.pos, (**(v[tt].begin())).second.pos, (**(v[tt].end() - 1)).first.pos + k - 1, (**(v[tt].end() - 1)).second.pos + k - 1); // (read_start, genome_start,read_end, genome_end)
			}
		}
	}
}


void 
SaveSplit (vector<vector<IndSeed>> &splitVSeed, FILE* yy) {
	if (yy == NULL) {
		perror("Error opening file: ");
	}
	else {
		for (unsigned int tt = 0; tt < splitVSeed.size(); ++tt) {
			if (! splitVSeed[tt].empty()) { // v[tt] is not empty
				fprintf(yy, "%lu %lu %lu %lu\n", beginPositionH(*(splitVSeed[tt].begin())), beginPositionV(*(splitVSeed[tt].begin())), endPositionH(*(splitVSeed[tt].end() - 1)), endPositionV(*(splitVSeed[tt].end() - 1))); // (read_start, genome_start,read_end, genome_end)
			}
		}
	}
}

void
SaveseedSet (IndSeedSet &seedSet, FILE* fd) {
	if (fd == NULL) {
		perror("Error opening file: ");
	}
	else {
		for (TIterator tt = begin(seedSet, seqan::Standard()); tt != end(seedSet, seqan::Standard()); ++tt) {
		fprintf(fd, "%lu %lu %lu %lu\n", beginPositionH(*tt), beginPositionV(*tt), endPositionH(*tt), endPositionV(*tt)); // (read_start, genome_start,read_end, genome_end)
		}		
	}
}



void 
SaveSparse (seqan::String<IndSeed> &chain, FILE* fi) {
	if (fi == NULL) {
		perror("Error opening file: ");
	}
	else {
		for (unsigned i = 0; i < length(chain); ++i) {
			fprintf(fi, "%lu %lu %lu %lu\n", beginPositionH(chain[i]), beginPositionV(chain[i]), endPositionH(chain[i]), endPositionV(chain[i]));
		}
	}
}



// This function merges seeds which satisfying situation 1 & situation 2
void 
merge (Cluster &rCr, vector<vector<GenomePair*>> &v, int h, int k) { // h is the threshold of the distance btwn read/genome																							
/*
    // Debug code
	const string filename1("/home/cmb-16/mjc/jingwenr/lra/lra_test/seeds.txt");  // file to store original_seed
    const string filename2("/home/cmb-16/mjc/jingwenr/lra/lra_test/merge.txt"); // file to store merge result(after removal)
 
    // save the original seeds
    FILE *fh = fopen(filename1.c_str(), "w");
    SaveOriginalSeed (rCr, fh, k);
    fclose(fh);
 */  

	for (vector<GenomePair>::iterator it = rCr.matches.begin(); it != rCr.matches.end(); ++it) {

		vector<GenomePair*> r;
		GenomePair * l;

		if (v.size() == 0) {
			l = &(*it);
			r.push_back(l);
			v.push_back(r);
		}
		else{
			unsigned int ii = 0;

			// criteria to merge *it in some v[ii]: if (situation 1 || situation 2), then merge them
	 		// situation 1: (on the same diagonal) && (overlap || nonoverlap and difference(read) <= h);
	 		// situation 2: (on different diagonal and the difference(diagonal) < 0.25*k) && (difference(read) <= h, differemce(genome) <= h) (avoid merging two unrelated seeds)
	 		//while (!situation 1) && (!situation 2); ++ii


			while (ii <= v.size() - 1 && ((Diagonalnum(*it) != Diagonalnum(**(v[ii].end() - 1)) || (*it).first.pos > h + (**(v[ii].end() -1)).first.pos + k - 1)
			&& (labs(Diagonalnum(**(v[ii].end() - 1)) - Diagonalnum(*it)) > 0.25*k || ((*it).first.pos > h + (**(v[ii].end() -1)).first.pos + k - 1|| 
			(**(v[ii].end() -1)).first.pos - k + 1  >= (*it).first.pos) || ((*it).second.pos > h + (**(v[ii].end() - 1)).second.pos + k - 1 || 
			(**(v[ii].end() - 1)).second.pos - k + 1 >= (*it).second.pos)))) {
				++ii;
			}

			l = &(*it);
			if (ii != v.size()) { // merge *it in v[ii]
				v[ii].push_back(l);
			}
			else{ // cannot merge *it in any v[]
				r.push_back(l);
				v.push_back(r);
			}
		}
	}

////////////////Debug Code
/*
for (unsigned int i = 0; i < v.size(); ++i){
	cout << "i: "<< i << endl;
	for (vector<GenomePair*>::iterator it = v[i].begin(); it != v[i].end(); ++it) {
		cout << "kmer: " << (**it).first.pos << " " << (**it).second.pos << " "<< (**it).first.pos + k - 1 << " "<< (**it).second.pos + k - 1 << " diagonal: " << (**it).first.pos - (**it).second.pos<< endl;
	}	
}
*/

///////////////Debug code
/*
    // save the merge_result
    FILE *fp = fopen(filename2.c_str(), "w");
    SaveLongSeed(v, fp, k);
    fclose(fp);
*/
}


// This function removes very short seeds (distances on read & genome are less than 100bp) and seeds with very low density(less than parameter g)
void 
RemoveSomeSeed (vector<vector<GenomePair*>> &v, vector<vector<GenomePair*>> &v_prime, int g, unsigned int k){ // g is the threshold of the density

	for (unsigned int i = 0; i < v.size(); ++i) {
		if ((**(v[i].end() - 1)).first.pos + k - 1 <= 20 + (**v[i].begin()).first.pos && (**(v[i].end() - 1)).second.pos + k - 1 <= 20 + (**v[i].begin()).second.pos) { // if the seed is too short, delete it
			// do nothing
		}
		else if (v[i].size() < g) { // if the v[i] is less than g density, delete it. 
			// do nothing
		}
		else { // insert v[i] into v_prime
			v_prime.push_back(v[i]);
		}
	}


////////////////Debug Code
/*
for (unsigned int i = 0; i < v_prime.size(); ++i){
	cout << "i: "<< i << endl;
	for (vector<GenomePair*>::iterator it = v_prime[i].begin(); it != v_prime[i].end(); ++it) {
		cout << "kmer: " << (**it).first.pos << " " << (**it).second.pos << " "<< (**it).first.pos + k - 1 << " "<< (**it).second.pos + k - 1 << " diagonal: " << (**it).first.pos - (**it).second.pos<< endl;
	}	
}
*/

}




void 
SplitH (vector<vector<GenomePair*>> &v_prime, vector<vector<IndSeed>> &splitHSeed, unsigned int k) {

	std::set<unsigned int> s;  //elements in s are arranged in strictly increasing order. and no duplicates
	for (unsigned int q = 0; q < v_prime.size(); ++q) { // insert the start_read and end_read coordinates into set s
		s.insert((**v_prime[q].begin()).first.pos);
		s.insert((**(v_prime[q].end() - 1)).first.pos + k - 1);
	}

// debug code 
/*
cout << "set: " << endl;
for (set<unsigned int>::iterator it = s.begin(); it != s.end(); ++it) {
	cout << "*it: " << *it << endl;
}
*/



	for (unsigned int q = 0; q < v_prime.size(); ++q) {
		vector<GenomePair*>::iterator it = v_prime[q].begin();
		std::set<unsigned int>::iterator low_1,low_2;
		low_1 = s.lower_bound ((**v_prime[q].begin()).first.pos);           
 		low_2 = s.lower_bound ((**(v_prime[q].end() - 1)).first.pos + k - 1);


 		std::set<unsigned int>::iterator l = low_1;
 		vector<GenomePair*>::iterator ty = v_prime[q].begin(); // **ty is the seed across the *l
 		advance(low_1, 1);


 		if (low_1 == low_2) { // no splitting on v_prime[q]
 			vector<IndSeed> h;
 			for (vector<GenomePair*>::iterator tt = v_prime[q].begin(); tt != v_prime[q].end(); ++tt) {
 				h.push_back(IndSeed((**tt).first.pos, (**tt).second.pos, (**tt).first.pos + k - 1, (**tt).second.pos + k - 1, splitHSeed.size()));
 			}
 			splitHSeed.push_back(h);

 		}
 		else {

 			std::set<unsigned int>::iterator f = low_2;
 			advance(low_2, 1);

 			for (std::set<unsigned int>::iterator lt = low_1; lt != low_2; ++lt) {
 				vector<IndSeed> h;

 				if (lt != f) {
 					it = ty; 
 					while((**it).first.pos + k - 1 <= *lt) {
 						++it;
 					}
 				}
 				else{
 					it = v_prime[q].end() - 1;
 				}

 				if (ty == it) {
 					if ((**ty).first.pos <= *lt) {  
 						if ((**ty).first.pos >= *l) {
 							h.push_back(IndSeed((**ty).first.pos, (**ty).second.pos, *lt, (unsigned int)(*lt - Diagonalnum(**ty)), splitHSeed.size()));
 						}
 						else {
 							h.push_back(IndSeed(*l, (unsigned int)(*l -Diagonalnum(**ty)), *lt, (unsigned int)(*lt - Diagonalnum(**ty)), splitHSeed.size()));

 						}
 					 	splitHSeed.push_back(h);
 					}
 				}
 				else {

 					//insert the 1-st seed
 					if ((**ty).first.pos >= *l) {
 						h.push_back(IndSeed((**ty).first.pos, (**ty).second.pos, (**ty).first.pos + k - 1, (**ty).second.pos + k - 1, splitHSeed.size()));
 					}
 					else {
 						h.push_back(IndSeed(*l, (unsigned int)(*l -Diagonalnum(**ty)), (**ty).first.pos + k - 1, (**ty).second.pos + k - 1, splitHSeed.size()));
 					}

 					// insert the middle seeds
 					if (ty + 1 != it) {
 						for (std::vector<GenomePair*>::iterator y = ty + 1; y != it; ++y) {
 							h.push_back(IndSeed((**y).first.pos, (**y).second.pos, (**y).first.pos + k - 1, (**y).second.pos + k - 1, splitHSeed.size()));
 						}
 					}

 					// insert the last seed
 					if ((**it).first.pos <= *lt) {
 						h.push_back(IndSeed((**it).first.pos, (**it).second.pos, *lt, (unsigned int)(*lt - Diagonalnum(**it)), splitHSeed.size()));
 					}

 				 	advance(l, 1);
 					ty = it;
 					splitHSeed.push_back(h);
 				}
 			}
 		}

 	}


////////////////Debug Code
/*
cout << "splitHSeed.size(): " << splitHSeed.size() << endl;
for (unsigned int i = 0; i < splitHSeed.size(); ++i){
	cout << "i: "<< i << endl;
	for (vector<IndSeed>::iterator it = splitHSeed[i].begin(); it != splitHSeed[i].end(); ++it) {
		cout << "longseed: *it" << *it << endl;
	}	
}
*/
}




void
SplitV (vector<vector<IndSeed>> &splitHSeed, vector<vector<IndSeed>> &splitVSeed) {
	std::set<unsigned int> s;  //elements in s are arranged in strictly increasing order. and no duplicates
	for (unsigned int q = 0; q < splitHSeed.size(); ++q) { // insert the start_read and end_read coordinates into set s
		s.insert(beginPositionV(*splitHSeed[q].begin()));
		s.insert(endPositionV(*(splitHSeed[q].end() - 1)));
	}

	for (unsigned int q = 0; q < splitHSeed.size(); ++q) {
		std::vector<IndSeed>::iterator it = splitHSeed[q].begin();
		std::set<unsigned int>::iterator low_1,low_2;
		low_1 = s.lower_bound (beginPositionV(*splitHSeed[q].begin()));           
 		low_2 = s.lower_bound (endPositionV(*(splitHSeed[q].end() - 1)));

		std::set<unsigned int>::iterator l = low_1;
 		std::vector<IndSeed>::iterator ty = splitHSeed[q].begin(); // *ty is the seed across the *l
 		advance(low_1, 1);


 		if (low_1 == low_2) { // no spliting on splitHSeed[q]


 			vector<IndSeed> h;
 			for (std::vector<IndSeed>::iterator tt = splitHSeed[q].begin(); tt != splitHSeed[q].end(); ++tt) {
 				h.push_back(IndSeed(beginPositionH(*tt), beginPositionV(*tt), endPositionH(*tt), endPositionV(*tt), splitVSeed.size()));
 			}
 			splitVSeed.push_back(h);
 		}
 		else {

 			std::set<unsigned int>::iterator f = low_2;
 			advance(low_2, 1);
 			for (std::set<unsigned int>::iterator lt = low_1; lt != low_2; ++lt) {
 				vector<IndSeed> h;
 				if (lt != f) {
 					it = ty; 
 					while(endPositionV(*it) <= *lt) {
 						++it;
 					}
 				}
 				else{
 					it = splitHSeed[q].end() - 1;
 				}


				if (ty == it) {
 					if (beginPositionV(*ty) <= *lt) {  
 						if (beginPositionV(*ty) >= *l) {
 							h.push_back(IndSeed(beginPositionH(*ty), beginPositionV(*ty), *lt + beginDiagonal(*ty), *lt, splitVSeed.size()));
 						}
 						else {
 							h.push_back(IndSeed(*l + beginDiagonal(*ty), *l, *lt + beginDiagonal(*ty), *lt, splitVSeed.size()));

 						}
 					 	splitVSeed.push_back(h);
 					}
 				}
 				else {
 					//insert the 1-st seed
 					if (beginPositionV(*ty) >= *l) {
 						h.push_back(IndSeed(beginPositionH(*ty), beginPositionV(*ty), endPositionH(*ty), endPositionV(*ty), splitVSeed.size()));
 					}
 					else {
 						h.push_back(IndSeed(*l + beginDiagonal(*ty), *l, endPositionH(*ty), endPositionV(*ty), splitVSeed.size()));
 					}

 					// insert the middle seeds
 					if (ty + 1 != it) {
 						for (std::vector<IndSeed>::iterator y = ty + 1; y != it; ++y) {
 							h.push_back(IndSeed(beginPositionH(*y), beginPositionV(*y), endPositionH(*y), endPositionV(*y), splitVSeed.size()));
 						}
 					}

 					// insert the last seed
 					if (beginPositionV(*it) <= *lt) {
 						h.push_back(IndSeed(beginPositionH(*it), beginPositionV(*it), *lt + beginDiagonal(*it), *lt, splitHSeed.size()));
 					}

 					advance(l, 1);
 					ty = it;
 					splitVSeed.push_back(h);
 				}
 			}
 		}
	 }


////////////////Debug Code
/*
cout << "splitVSeed.size(): " << splitVSeed.size() << endl;
for (unsigned int i = 0; i < splitVSeed.size(); ++i){
	cout << "i: "<< i << endl;
	for (vector<IndSeed>::iterator it = splitVSeed[i].begin(); it != splitVSeed[i].end(); ++it) {
		cout << "longseed: *it" << *it << endl;
	}	
}
*/
/*
    // debug code
    const string filename3("/home/cmb-16/mjc/jingwenr/lra/lra_test/split.txt"); // file to store split result
    // save the split_result
    FILE *yy = fopen(filename3.c_str(), "w");
    SaveSplit (splitVSeed, yy);
    fclose(yy);
*/

}


void 
AddInset (vector<vector<IndSeed>> &splitVSeed, IndSeedSet &seedSet) {
	for (unsigned int i = 0; i != splitVSeed.size(); ++i) {
		seqan::addSeed(seedSet, IndSeed(beginPositionH(*splitVSeed[i].begin()), beginPositionV(*splitVSeed[i].begin()), endPositionH(*(splitVSeed[i].end() - 1)), endPositionV(*(splitVSeed[i].end() - 1)), i), seqan::Single());
	}
}


#endif
