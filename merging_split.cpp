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
#include <seqan/sequence.h>

using namespace std;
using namespace seqan;

typedef seqan::Seed<Simple> TSeed;

// read fasta file 
static void
read_fasta_file_single_sequence(const string &filename, string &T) {

	ifstream in(filename.c_str());
	string line;
	in >> line;
    while (in >> line)
    	T += line;
}

unsigned int
linefunc (unsigned int x_0, unsigned int y_0, unsigned int x_1, unsigned int y_1, unsigned int x, int b) {
	float a;
	float x_00 = static_cast<float>(x_0);
	float x_10 = static_cast<float>(x_1);
	float y_00 = static_cast<float>(y_0);
	float y_10 = static_cast<float>(y_1);
	float x0 = static_cast<float>(x);
	if (b == 1) { // input x is a read coordinate, return the genome coordinate


		a = floor(((y_10 - y_00)/(x_10 - x_00))*x0 - (x_00*y_10 - x_10*y_00)/(x_10 - x_00));
		unsigned int a_prime = static_cast<unsigned int>(a);

		return a_prime;
	}
	else { // input x is a genome coordinate, return the read coordinate
/*
////////////////////////
	cout << "(x_10 - x_00): " << (x_10 - x_00) << endl;
	cout << "(y_10 - y_00)): " << (y_10 - y_00) << endl;

	cout << "((x_10 - x_00)/(y_10 - y_00)): " << ((x_10 - x_00)/(y_10 - y_00)) << endl;
/////////////////////////
*/
		a = floor(((x_10 - x_00)/(y_10 - y_00))*(x0 + (x_00*y_10 - x_10*y_00)/(x_10 - x_00)));
		unsigned int a_prime = static_cast<unsigned int>(a);
		return a_prime;
	}
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
	 		//cout << "v.size() == 0" << endl;/////////////////
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
	 			//cout << "ii != v.size() " << endl;
	 			v[ii].push_back(l);
	 		}
	 		else { // cannot merge *it in any v[]
	 			//cout << "cannot merge *it in any v[], create a new v[]" << endl;
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



// This function 
void
SplitH (vector<vector<TSeed*>> &v_prime, vector<vector<TSeed>> &splitHSeed, const unsigned int &m) {// m is the split length threshold

	std::set<unsigned int> s; // elements in s are arranged in strictly increasing order. and no duplicates
	for (unsigned int p = 0; p < v_prime.size(); ++p) { // insert the start_read and end_read coordinates into set s
		s.insert(beginPositionH(**v_prime[p].begin()));
		s.insert(endPositionH(**(v_prime[p].end() - 1)));
	}

	for (unsigned int q = 0; q < v_prime.size(); ++q) {
		std::set<unsigned int>::iterator low_1,low_2;
		low_1 = s.lower_bound (beginPositionH(**v_prime[q].begin()));           
 		low_2 = s.lower_bound (endPositionH(**(v_prime[q].end() - 1)));
 		std::set<unsigned int>::iterator l = low_1;

 		vector<TSeed> h;
 		advance(low_1, 1);

 		if (low_1 == low_2) { // actually means low_1 + 1 == low_2
 			// no splitting on long seed *v_prime[q]

 			h.push_back(TSeed(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), endPositionH(**(v_prime[q].end() - 1)), endPositionV(**(v_prime[q].end() - 1))));
 			splitHSeed.push_back(h);
 		}
 		else {
 			int b = 1;
 			// insert the 1-st piece
 			if (*low_1 - *l >= m) { // m is the threshold of split-length
 				unsigned int endPV = linefunc(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), endPositionH(**(v_prime[q].end() - 1)), endPositionV(**(v_prime[q].end() - 1)), *low_1 - 1, b);
 				h.push_back(TSeed(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), *low_1 -1, endPV)); 
 				splitHSeed.push_back(h);
 			} 

 			advance(l, 1);
 			advance(low_1, 1);


 			if (low_1 != low_2) {

 				for (std::set<unsigned int>::iterator it = low_1; it != low_2; ++it) {
 					if (*it - *l >= m) {

 						if (splitHSeed.size() == q) { // splitHSeed[q] is empty						
 							unsigned int endPV = linefunc(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), endPositionH(**(v_prime[q].end() - 1)), endPositionV(**(v_prime[q].end() - 1)), *it -1, b);	
							unsigned int endPV_1 = linefunc(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), endPositionH(**(v_prime[q].end() - 1)), endPositionV(**(v_prime[q].end() - 1)), *l, b);
							h.push_back(TSeed(*l, endPV_1, *it -1, endPV));
							splitHSeed.push_back(h);
 						}
 						else {
 							if (*l == endPositionH(*(splitHSeed[q].end() - 1)) + 1) {
 								unsigned int endPV = linefunc(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), endPositionH(**(v_prime[q].end() - 1)), endPositionV(**(v_prime[q].end() - 1)), *it -1, b);
 								splitHSeed[q].push_back(TSeed(endPositionH(*(splitHSeed[q].end() - 1)) + 1, endPositionV(*(splitHSeed[q].end() - 1)) + 1, *it -1, endPV));
 							}
 							else {
 								unsigned int endPV = linefunc(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), endPositionH(**(v_prime[q].end() - 1)), endPositionV(**(v_prime[q].end() - 1)), *it -1, b);	
								unsigned int endPV_1 = linefunc(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), endPositionH(**(v_prime[q].end() - 1)), endPositionV(**(v_prime[q].end() - 1)), *l, b);
 								splitHSeed[q].push_back(TSeed(*l, endPV_1, *it -1, endPV));					
 							}
 						}
 					}
 					advance(l, 1);
 				}

 			} 

 			// insert the last piece
 			if (*low_2 - *l >= m) {

 				if (splitHSeed.size() == q) {
 					unsigned int endPV_1 = linefunc(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), endPositionH(**(v_prime[q].end() - 1)), endPositionV(**(v_prime[q].end() - 1)), *l, b);
 					h.push_back(TSeed(*l, endPV_1, *low_2, endPositionV(**(v_prime[q].end() - 1))));
 					splitHSeed.push_back(h);
 				}
 				else {
 					if (*l == endPositionH(*(splitHSeed[q].end() - 1)) + 1) {
 						 splitHSeed[q].push_back(TSeed(endPositionH(*(splitHSeed[q].end() - 1)) + 1, endPositionV(*(splitHSeed[q].end() - 1)) + 1, *low_2, endPositionV(**(v_prime[q].end() - 1))));
 					}
 					else {
 						unsigned int endPV_1 = linefunc(beginPositionH(**v_prime[q].begin()), beginPositionV(**v_prime[q].begin()), endPositionH(**(v_prime[q].end() - 1)), endPositionV(**(v_prime[q].end() - 1)), *l, b);
 						splitHSeed[q].push_back(TSeed(*l, endPV_1, *low_2, endPositionV(**(v_prime[q].end() - 1))));
 					}
 				}
 			}
 		}

 		// if splitHSeed[q] = NULL	
 		if (splitHSeed.size() == q) {
 			splitHSeed.push_back(h);
 		}

	}
}







void 
RowiseHSeed (vector<vector<TSeed>> &splitHSeed, vector<TSeed*> &HSeed) {
	for (unsigned int p = 0; p < splitHSeed.size(); ++p) {
		TSeed *l;
		for (vector<TSeed>::iterator it = splitHSeed[p].begin(); it != splitHSeed[p].end(); ++it) {
			l = &(*it);
			HSeed.push_back(l);
		}
	}
}





void
SplitV (vector<TSeed*> &HSeed, vector<vector<TSeed>> &splitVSeed, const unsigned int &m) {

	std::set<unsigned int> s; 
	for (unsigned int p = 0; p < HSeed.size(); ++p) { // insert the startV and endV coordinates into set s

		s.insert(beginPositionV(*(HSeed[p])));
		s.insert(endPositionV(*(HSeed[p])));
	}

	for (unsigned int q = 0; q < HSeed.size(); ++q) {
		std::set<unsigned int>::iterator low_1,low_2;
		low_1 = s.lower_bound (beginPositionV(*(HSeed[q])));               
 		low_2 = s.lower_bound (endPositionV(*(HSeed[q])));
 		std::set<unsigned int>::iterator l = low_1;

 		vector<TSeed> h;
 		advance(low_1, 1);

 		if (low_1 == low_2) { // actually means low_1 + 1 == low_2
 			// no splitting on long seed *v_prime[q]
 			h.push_back(TSeed(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPositionH(*(HSeed[q])), endPositionV(*(HSeed[q]))));
 			splitVSeed.push_back(h);
 		}
 		else {

 			// insert 1-st piece 
 			int b = 0;
 			if (*low_1 - *l >= m) {
 				unsigned int endPH = linefunc(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPositionH(*(HSeed[q])), endPositionV(*(HSeed[q])),*low_1 - 1, b);
 				h.push_back(TSeed(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPH, *low_1 - 1)); 
 				splitVSeed.push_back(h);
 			}

 			advance(l, 1);
 			advance(low_1,1);




 			if (low_1 != low_2) {

 				for (std::set<unsigned int>::iterator it = low_1; it != low_2; ++it) {

					if (*it - *l >= m) {

 						if (splitVSeed.size() == q) { // splitHSeed[q] is empty						
 							unsigned int endPH = linefunc(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPositionH(*(HSeed[q])), endPositionV(*(HSeed[q])), *it -1, b);	
							unsigned int endPH_1 = linefunc(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPositionH(*(HSeed[q])), endPositionV(*(HSeed[q])), *l, b);
							h.push_back(TSeed(endPH_1, *l, endPH, *it - 1));
							splitVSeed.push_back(h);
 						}
 						else {
 							if (*l == endPositionV(*(splitVSeed[q].end() - 1)) + 1) {
 								unsigned int endPH = linefunc(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPositionH(*(HSeed[q])), endPositionV(*(HSeed[q])), *it -1, b);
 								splitVSeed[q].push_back(TSeed(endPositionH(*(splitVSeed[q].end() - 1)) + 1, endPositionV(*(splitVSeed[q].end() - 1)) + 1, endPH, *it - 1));
 							}
 							else {
 								unsigned int endPH = linefunc(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPositionH(*(HSeed[q])), endPositionV(*(HSeed[q])), *it -1, b);	
								unsigned int endPH_1 = linefunc(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPositionH(*(HSeed[q])), endPositionV(*(HSeed[q])), *l, b);
 								splitVSeed[q].push_back(TSeed(endPH_1, *l, endPH, *it - 1));						
 							}
 						}
 					}

 				advance(l, 1);
 				}
 			} 
			// insert the last piece
 			if (*low_2 - *l >= m) {

 				if (splitVSeed.size() == q) { // splitVSeed[q] is NULL
 					unsigned int endPH_1 = linefunc(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPositionH(*(HSeed[q])), endPositionV(*(HSeed[q])), *l, b);
 					h.push_back(TSeed(endPH_1, *l, endPositionH(*(HSeed[q])), *low_2));
 					splitVSeed.push_back(h);
 				}
 				else {// splitVSeed[q] is not NULL
 					 // the previous one has been insert
 					if (*l == endPositionV(*(splitVSeed[q].end() - 1)) + 1) {
 						 splitVSeed[q].push_back(TSeed(endPositionH(*(splitVSeed[q].end() - 1)) + 1, endPositionV(*(splitVSeed[q].end() - 1)) + 1, endPositionH(*(HSeed[q])), *low_2));
 					}
 					else {// the previous one has not been inserted
 						unsigned int endPH_1 = linefunc(beginPositionH(*(HSeed[q])), beginPositionV(*(HSeed[q])), endPositionH(*(HSeed[q])), endPositionV(*(HSeed[q])), *l, b);
 						splitVSeed[q].push_back(TSeed(endPH_1, *l, endPositionH(*(HSeed[q])), *low_2));
 					}
 				}
 			}
 			
 		}	
	}

}




void 
SaveOriginalSeed (vector<TSeed> &v, FILE* fh) {  
	if (fh == NULL) {
		perror("Error opening file: ");
	}
	else {
		for (unsigned int tt = 0; tt < v.size(); ++tt) {
			fprintf(fh, "%lu %lu %lu %lu\n", beginPositionH(v[tt]), beginPositionV(v[tt]), endPositionH(v[tt]), endPositionV(v[tt])); // (read_start, genome_start,read_end, genome_end)
		}
	}
}


void 
SaveLongSeed (vector<vector<TSeed*>> &v, FILE* fp) {  
	if (fp == NULL) {
		perror("Error opening file: ");
	}
	else {
		for (unsigned int tt = 0; tt < v.size(); ++tt) {
			fprintf(fp, "%lu %lu %lu %lu\n", beginPositionH(**v[tt].begin()), beginPositionV(**v[tt].begin()), endPositionH(**(v[tt].end() - 1)), endPositionV(**(v[tt].end() - 1))); // (read_start, genome_start,read_end, genome_end)
		}
	}
}


void 
SaveSplit (vector<vector<TSeed>> &v, FILE* yy) {
	if (yy == NULL) {
		perror("Error opening file: ");
	}
	else {
		for (unsigned int tt = 0; tt < v.size(); ++tt) {
			for (vector<TSeed>::iterator it = v[tt].begin(); it != v[tt].end(); ++it) {
				fprintf(yy, "%lu %lu %lu %lu\n", beginPositionH(*it), beginPositionV(*it), endPositionH(*it), endPositionV(*it)); // (read_start, genome_start,read_end, genome_end)
			}
		}
	}
}



int 
main (int argc, const char * const argv[]) {

   if (argc != 10) {
      cerr << "NEEDS INPUTS: " << argv[0] << " <FASTA-FILE1>/read and <FASTA-FILE2>/genome and  <PATH TO STORE ORIGINAL SEEDS> and <PATH TO STORE MERGE RESULT> and <PATH TO STORE SPLIT SEEDS> and K-MER LENGTH and MERGING THRESHOLD h and DENSITY THRESHOLD g and SPLIT LENGTH THRESHOLD" << endl;
      return EXIT_FAILURE;
    }


   const string filename1(argv[1]);  // read
   const string filename2(argv[2]);  // genome
   const string filename3(argv[3]);  // file to store original_seed
   const string filename4(argv[4]); // file to store merge result(after removal)
   const string filename5(argv[5]); // file to store split result


   const unsigned int k = atoi(argv[6]); // kmer length 
   const unsigned int h = atoi(argv[7]); //  merge threshold (distance btwn read and gemone) 1000
   const unsigned int g = atoi(argv[8]); // the density threshold of long seed 100
   const unsigned int m = atoi(argv[9]);

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

   vector<vector<TSeed>> splitHSeed;
   SplitH(v_prime, splitHSeed, m);

   vector<TSeed*> HSeed;
   RowiseHSeed(splitHSeed, HSeed); // HSeed stores all seeds in splitHSeed

   vector<vector<TSeed>> splitVSeed;
   SplitV(HSeed, splitVSeed, m); // splitVSeed stores all seeds we want



// save the original seeds
   FILE *fh = fopen(filename3.c_str(), "w");
   SaveOriginalSeed (seedVector, fh);
   fclose(fh);


// save the merge_result
   FILE *fp = fopen(filename4.c_str(), "w");
   SaveLongSeed(v_prime, fp);
   fclose(fp);

// save the split_result
   FILE *yy = fopen(filename5.c_str(), "w");
   SaveSplit (splitVSeed, yy);
   fclose(yy);
 
  return 0;
}
