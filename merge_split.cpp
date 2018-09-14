// generating k-mers of fixed length + merging k-mers into longer fragments + preparing for spase dp

#include <iostream> //std::cout 
#include <fstream>   
#include <cstdlib>   // std::labs, std::EXIT FAILURE, std::EXIT SUCCESS
#include <string> 
#include <cmath>        // std::labs
#include <cstdio>    // std::FILE std::perror
#include <algorithm> // std::max
#include <vector>
#include <set>
#include "seqan/seeds.h"
#include "seqan/align.h"


using namespace std;
using namespace seqan;


typedef seqan::Seed<Simple> TSeed;

struct IndexedSeed_;
typedef seqan::Tag<IndexedSeed_> IndexedSeed;

template<typename TConfig>
class seqan::Seed<IndexedSeed, TConfig> : public seqan::Seed<Simple, TConfig>{
public:
	int index;
	Seed(int i, int j, int k, int l) : seqan::Seed<Simple,TConfig>(i,j,k,l){index=-1;}
	Seed(int i, int j, int k) : seqan::Seed<Simple,TConfig>(i,j,k){index=-1;}
	Seed() : seqan::Seed<Simple,TConfig>(){index=-1;}	
	Seed(int i, int j, int k, int l, int idx) : seqan::Seed<Simple,TConfig>(i,j,k,l), index(idx) {}
};


typedef seqan::Seed<IndexedSeed> IndSeed;

typedef SeedSet<IndSeed> IndSeedSet;

template <typename TConfig>
typename Position<Seed<IndexedSeed, TConfig> >::Type
beginPositionH(Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._beginPositionH;
}


template <typename TConfig>
typename Position<Seed<IndexedSeed, TConfig> >::Type
beginPositionV(Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._beginPositionV;
}


template <typename TConfig>
typename Position<Seed<IndexedSeed, TConfig> >::Type
endPositionH(Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._endPositionH;
}


template <typename TConfig>
typename Position<Seed<IndexedSeed, TConfig> >::Type
endPositionV(Seed<IndexedSeed, TConfig> const & seed)
{
    return seed._endPositionV;
}


template <typename TStream, typename TConfig>
inline TStream &
operator<<(TStream & stream, Seed<IndexedSeed, TConfig> const & seed)
{
    return stream << "Seed<IndexedSeed, TConfig>(" << beginPositionH(seed)
                  << ", " << beginPositionV(seed) << ", "
                  << endPositionH(seed) << ", "
                  << endPositionV(seed) << ", lower diag = "
                  << lowerDiagonal(seed) << ", upper diag = "
                  << upperDiagonal(seed) << ")";
}


// read fasta file 
static void
read_fasta_file_single_sequence(const string &filename, string &T) {

	ifstream in(filename.c_str());
	string line;
	in >> line;
    while (in >> line)
    	T += line;
}


// generate k-mers and add them into seedSet
void
find_match (const string &read, const string &genome, const unsigned int &k, vector<TSeed> &seedVector) { 
	unsigned t = 0; // record the number of matches
	cout << "read length: " << read.length()<< endl; /////////
	cout << "genome length: " << genome.length() << endl; ////////////
	for (unsigned i = 0; i < read.length(); ++i) {
		for (unsigned j = 0; j < genome.length(); ++j) {
			t = 0;
			while((genome[j + t] == read[i + t] || genome[j + t] == read[i + t] + 32 || genome[j + t] + 32 == read[i + t] || genome[j + t] + 32 == read[i + t] + 32) && t < k) { 
				++t;
			}
			if (t == k) {  // a k-mer (i, j, i + k -1, j + k - 1) (read_start, genome_start, read_end, genome_end)
				seedVector.push_back(TSeed(i, j, k));
			}
		}
	}
}


// This function merges seeds which satisfying situation 1 & situation 2
void 
merge (vector<TSeed> &seedVector, vector<vector<TSeed*>> &v, const unsigned int &h, const unsigned int &k) { // h is distance threshold btwn read/genome
	 
	 for (vector<TSeed>::iterator it = seedVector.begin(); it != seedVector.end(); ++it) {
	 	vector<TSeed*> r;
	 	TSeed * l;
	 	if (v.size() == 0) {
	 		l = &(*it);
	 		r.push_back(l); 
	 		v.push_back(r);// make a copy of r
	 	}
	 	else { // for seed *it, try to find a mergepartner in vector<vector<TSeed>> v for it. 
	 		// find the first one in vector<vector<TSeed>> v, so that difference(diagonal) <= 0.25*k; difference(read) <= h; difference(genome) <= h
	 		// ah=and they a
	 		unsigned int ii = 0;
	 		// criteria to merge *it in some v[ii]: if (situation 1 || situation 2), then merge them
	 		// situation 1: on the same diagonal, and  labs(difference(read)) <= h;
	 		// situation 2: on different diagonal(the difference(diagonal) < 0.25.k), difference(read) <= h, differemce(genome) <= h
	 		// situation 1: same diagonal, overlap or non-overlap but within threshold h
	 		// situation 2: different diagonal, but the difference between different diagonals is within threshold 0.25*k. 
	 		// and the distance between reads and genomes are within threshold h. (avoid merging two unrelated seeds
	 		//while (!situation 1) && (!situation 2); ++ii

	 		while (ii <= v.size() - 1 && ((beginDiagonal(**(v[ii].end() - 1)) != beginDiagonal(*it) || labs(beginPositionH(*it) - endPositionH(**(v[ii].end() - 1))) > h) &&
	 		 (labs(beginDiagonal(**(v[ii].end() - 1)) - beginDiagonal(*it)) > 0.25*k || beginPositionH(*it) - endPositionH(**(v[ii].end() - 1)) > h ||
	 		 	beginPositionV(*it) - endPositionV(**(v[ii].end() - 1)) > h) )) {
	 			++ii;
	 		}
	 		l = &(*it);
	 		if (ii != v.size()) { // merge *it in v[ii]
	 			v[ii].push_back(l);
	 		}
	 		else { // cannot merge *it in any v[]
	 			r.push_back(l);
	 			v.push_back(r);	
	 		} 
	 	} 
	 }
}


// This function removes very short seeds (distances on read & genome are less than 100bp) and seeds with very low density(less than parameter g)
// TODO(Jingwen): delete vector<vector<TSeed*>> v, which won't be used anymore.
void 
RemoveSomeSeed (vector<vector<TSeed*>> &v, vector<vector<TSeed*>> &v_prime, const unsigned int &g){ // g is the threshold of the density
	for (unsigned int i = 0; i < v.size(); ++i) {
		if (beginPositionH(**(v[i].end() - 1)) - beginPositionH(**v[i].begin()) <= 100 && beginPositionV(**(v[i].end() - 1)) - beginPositionV(**v[i].begin()) <= 100) { // if the seed is too short, delete it
			// do nothing
		}
		else if (v[i].size() < g) { // if the v[i] is less than g density, delete it. 
			// do nothing
		}
		else { // insert v[i] into v_prime
			v_prime.push_back(v[i]);
		}
	}

}


void 
SplitH (vector<vector<TSeed*>> &v_prime, vector<vector<IndSeed>> &splitHSeed) {// m is the split length threshold

	std::set<unsigned int> s;  //elements in s are arranged in strictly increasing order. and no duplicates
	for (unsigned int q = 0; q < v_prime.size(); ++q) { // insert the start_read and end_read coordinates into set s
		s.insert(beginPositionH(**v_prime[q].begin()));
		s.insert(endPositionH(**(v_prime[q].end() - 1)));
	}

	for (unsigned int q = 0; q < v_prime.size(); ++q) {
		std::vector<TSeed*>::iterator it = v_prime[q].begin();
		std::set<unsigned int>::iterator low_1,low_2;
		low_1 = s.lower_bound (beginPositionH(**v_prime[q].begin()));           
 		low_2 = s.lower_bound (endPositionH(**(v_prime[q].end() - 1)));


 		std::set<unsigned int>::iterator l = low_1;
 		std::vector<TSeed*>::iterator ty = v_prime[q].begin(); // **ty is the seed across the *l
 		advance(low_1, 1);


 		if (low_1 == low_2) { // no splitting on v_prime[q]
 			vector<IndSeed> h;
 			for (std::vector<TSeed*>::iterator tt = v_prime[q].begin(); tt != v_prime[q].end(); ++tt) {
 				h.push_back(IndSeed(beginPositionH(**tt), beginPositionV(**tt), endPositionH(**tt), endPositionV(**tt), splitHSeed.size()));
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
 					while(endPositionH(**it) <= *lt) {
 						++it;
 					}
 				}
 				else{
 					it = v_prime[q].end() - 1;
 				}


 				if (ty == it) {
 					if (beginPositionH(**ty) <= *lt) {  
 						if (beginPositionH(**ty) >= *l) {
 							h.push_back(IndSeed(beginPositionH(**ty), beginPositionV(**ty), *lt, *lt - beginDiagonal(**ty), splitHSeed.size()));
 						}
 						else {
 							h.push_back(IndSeed(*l, *l -beginDiagonal(**ty), *lt, *lt - beginDiagonal(**ty), splitHSeed.size()));

 						}
 					 	splitHSeed.push_back(h);
 					}
 				}
 				else {

 					//insert the 1-st seed
 					if (beginPositionH(**ty) >= *l) {
 						h.push_back(IndSeed(beginPositionH(**ty), beginPositionV(**ty), endPositionH(**ty), endPositionV(**ty), splitHSeed.size()));
 					}
 					else {
 						h.push_back(IndSeed(*l, *l -beginDiagonal(**ty), endPositionH(**ty), endPositionV(**ty), splitHSeed.size()));
 					}

 					// insert the middle seeds
 					if (ty + 1 != it) {
 						for (std::vector<TSeed*>::iterator y = ty + 1; y != it; ++y) {
 							h.push_back(IndSeed(beginPositionH(**y), beginPositionV(**y), endPositionH(**y), endPositionV(**y), splitHSeed.size()));
 						}
 					}

 					// insert the last seed
 					if (beginPositionH(**it) <= *lt) {
 						h.push_back(IndSeed(beginPositionH(**it), beginPositionV(**it), *lt, *lt - beginDiagonal(**it), splitHSeed.size()));
 					}

 				 	advance(l, 1);
 					ty = it;
 					splitHSeed.push_back(h);
 				}
 			}
 		}

 	}
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
}


void 
AddInset (vector<vector<IndSeed>> &splitVSeed, IndSeedSet seedSet) {
	for (unsigned int i = 0; i != splitVSeed.size(); ++i) {
		addSeed(seedSet, IndSeed(beginPositionH(*splitVSeed[i].begin()), beginPositionV(*splitVSeed[i].begin()), endPositionH(*(splitVSeed[i].end() - 1)), endPositionV(*(splitVSeed[i].end() - 1)), i), Single());
	}
}


int 
main (int argc, const char * const argv[]) {

   if (argc != 5) {
      cerr << "NEEDS INPUTS: " << argv[0] << " <FASTA-FILE1>/read and <FASTA-FILE2>/genome and K-MER LENGTH and MERGING THRESHOLD h and DENSITY THRESHOLD g" << endl;
      return EXIT_FAILURE;
    }


   const string filename1(argv[1]);  // read
   const string filename2(argv[2]);  // genome

   const unsigned int k = atoi(argv[3]); // kmer length 
   const unsigned int h = atoi(argv[4]); //  merge threshold (distance btwn read and gemone) 1000
   const unsigned int g = atoi(argv[5]); // the density threshold of long seed 100

   std::string read;
   std::string genome;
   read_fasta_file_single_sequence(filename1, read);
   read_fasta_file_single_sequence(filename2, genome);

   vector<TSeed> seedVector;
   find_match(read, genome, k, seedVector);

   vector<vector<TSeed*>> v;
   merge(seedVector, v, h, k);


   vector<vector<TSeed*>> v_prime; // v_prime[i] means i-th long seed
   RemoveSomeSeed(v, v_prime, g);


   vector<vector<IndSeed>> splitHSeed;
   SplitH(v_prime, splitHSeed);

 
   vector<vector<IndSeed>> splitVSeed;
   SplitV(splitHSeed, splitVSeed); // splitVSeed stores all seeds we want

   IndSeedSet seedSet;
   AddInset(splitVSeed, seedSet);

  return 0;
}