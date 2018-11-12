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

MergeSplit () {
	//------------------------------------------------------------------
   	// Step 1: merge matches with diagnol difference smaller than maxDiag
    //------------------------------------------------------------------
    vector<Cluster> v;
    vector<Cluster> vq;
    vector<Cluster> vt;
    opts.maxGap = -1;
    opts.maxDiag = 10;

    StoreDiagonalClusters(refinedClusters[r].matches, v, opts, true); 
    std::set<unsigned int> s1;     // s1 stores x boundary; s2 stores y boundary      
    std::set<unsigned int> s2;    //elements in s are arranged in strictly increasing order. and no duplicates

    //-------------- debug code
    //cout << "v.size(): " << v.size() << endl;
    /*
	// Debug code ----------- print out "orignal"
	if (r == 10) {
		seqan::SeedSet<IndSeed, seqan::Unordered> seedSet1;
	for (unsigned int i = 0; i != v.size(); ++i) {
		seqan::addSeed(seedSet1, IndSeed(v[i].qStart, v[i].tStart, v[i].qEnd, v[i].tEnd, i), seqan::Single());
	}  
	const string filename1("/home/cmb-16/mjc/jingwenr/lra/lra1/original.txt");  
	FILE *fr = fopen(filename1.c_str(), "w");
	SaveseedSet (seedSet1, fr); 
	fclose(fr);			
	}
	*/
			

    //---------------------------------------------------------------------------------
    // Step 2: sort matches by x and y coordinates in each Cluster v[i](merged matches)
    //---------------------------------------------------------------------------------
    for (unsigned int i = 0; i < v.size(); ++i) {

    	//insert x y boundaries into s1 and s2
    	s1.insert(v[i].qStart);
    	s1.insert(v[i].qEnd);
    	s2.insert(v[i].tStart);
    	s2.insert(v[i].tEnd);

    	// sort matches by x and y inside each merged matches cluster
    	vector<pair<GenomeTuple, GenomeTuple>>::iterator begin = refinedClusters[r].matches.begin();
    	if (v[i].start != 0) {
        	std::advance(begin , v[i].start);
    	} 

    	vector<pair<GenomeTuple, GenomeTuple>>::iterator end = refinedClusters[r].matches.begin();
    	if (v[i].end != 0) {
            std::advance(end, v[i].end);
    	}
    	CartesianSort<GenomeTuple>(begin, end);            
    }

    //----------------------------------------------------------------------------
    // Step 3:  split based on x boundaries
    //----------------------------------------------------------------------------
    for (unsigned int i = 0; i < v.size(); ++i) {
    	std::set<unsigned int>::iterator qs, qe; 
        qs = s1.lower_bound(v[i].qStart);
        qe = s1.lower_bound(v[i].qEnd);
        std::set<unsigned int>::iterator l = qs;
		std::advance(l, 1);
        unsigned int ii = v[i].start;
        unsigned int jj = ii;



        if (*l == *qe) {
            vq.push_back(Cluster(v[i].start, v[i].end, v[i].qStart, v[i].qEnd, v[i].tStart, v[i].tEnd, 0));
        }
        else{
            while (ii != v[i].end && l != qe) { // end is end + 1!!!!!!
            	while (refinedClusters[r].matches[ii].first.pos <= *l && ii != v[i].end) {
                    if (refinedClusters[r].matches[ii].first.pos + opts.globalK - 1 >= *l) {
                        // split 
                        if (jj + 1 <= ii) {
                            vq.push_back(Cluster(jj, ii, refinedClusters[r].matches[jj].first.pos, refinedClusters[r].matches[ii - 1].first.pos + opts.globalK,
                            	refinedClusters[r].matches[jj].second.pos, refinedClusters[r].matches[ii - 1].second.pos + opts.globalK, 0));
                        }
                            vq.push_back(Cluster(ii, ii + 1, refinedClusters[r].matches[ii].first.pos, refinedClusters[r].matches[ii].first.pos + opts.globalK,
                        		refinedClusters[r].matches[ii].second.pos, refinedClusters[r].matches[ii].second.pos + opts.globalK, 0));
                        	jj = ii + 1;
                    }
                    ++ii;
                }
                std::advance(l, 1);
            }	
           	if (ii != v[i].end && l == qe) {
               	vq.push_back(Cluster(jj, v[i].end, refinedClusters[r].matches[jj].first.pos, v[i].qEnd, refinedClusters[r].matches[jj].second.pos, v[i].tEnd, 0));
            }

        }
   	}

            //------------- debug code 
            //cout << "vq.size(): " << vq.size() << endl;

            /*
			// Debug code ----------- print out "seedSet"
			if (r == 10) {
			seqan::SeedSet<IndSeed, seqan::Unordered> seedSet2;
			for (unsigned int i = 0; i != vq.size(); ++i) {
				seqan::addSeed(seedSet2, IndSeed(vq[i].qStart, vq[i].tStart, vq[i].qEnd, vq[i].tEnd, i), seqan::Single());
			}  
			const string filename2("/home/cmb-16/mjc/jingwenr/lra/lra1/seedSet1.txt");  
			FILE *fm = fopen(filename2.c_str(), "w");
			SaveseedSet (seedSet2, fm); 
			fclose(fm);			
			}
			*/



    //----------------------------------------------------------------------------
   	// Step 4:  split based on y boundaries
   	//----------------------------------------------------------------------------
    for (unsigned int i = 0; i < vq.size(); ++i) {
        std::set<unsigned int>::iterator ts, te; 
        ts = s2.lower_bound(vq[i].tStart);
        te = s2.lower_bound(vq[i].tEnd);
        std::set<unsigned int>::iterator l = ts;
        std::advance(l, 1);
        unsigned int ii = vq[i].start;
        unsigned int jj = ii;

        if (*l == *ts) {
            vt.push_back(Cluster(vq[i].start, vq[i].end, vq[i].qStart, vq[i].qEnd, vq[i].tStart, vq[i].tEnd, 0));
        }
        else {
            while (ii != vq[i].end && l != te) { 
                while (refinedClusters[r].matches[ii].second.pos <= *l && ii != vq[i].end) {
                    if (refinedClusters[r].matches[ii].second.pos + opts.globalK - 1 >= *l) {
                        // split 
                        if (jj + 1 <= ii) {
                            vt.push_back(Cluster(jj, ii, refinedClusters[r].matches[jj].first.pos, refinedClusters[r].matches[ii - 1].first.pos + opts.globalK,
                            	refinedClusters[r].matches[jj].second.pos, refinedClusters[r].matches[ii - 1].second.pos + opts.globalK, 0));
                        }
                        vt.push_back(Cluster(ii, ii + 1, refinedClusters[r].matches[ii].first.pos, refinedClusters[r].matches[ii].first.pos + opts.globalK,
                        	refinedClusters[r].matches[ii].second.pos, refinedClusters[r].matches[ii].second.pos + opts.globalK, 0));
                        jj = ii + 1;
                    }
                	++ii;
                }
                std::advance(l, 1);
            }
            if (ii != vq[i].end && l == te) {
                		vt.push_back(Cluster(jj, vq[i].end, refinedClusters[r].matches[jj].first.pos, vq[i].qEnd, refinedClusters[r].matches[jj].second.pos, vq[i].tEnd, 0));
                	}

        }
    }
            //------------ debug
            //cout << "vt.size(): " << vt.size() << endl;


    //----------------------------------------------------------------------------
            // Step 5:  Store the result in seedSet
            //----------------------------------------------------------------------------
            // Not necessarily use seedSet as input of NaiveDp 
	for (unsigned int i = 0; i != vt.size(); ++i) {
		seqan::addSeed(seedSet, IndSeed(vt[i].qStart, vt[i].tStart, vt[i].qEnd, vt[i].tEnd, i), seqan::Single());
	}  



}

#endif
