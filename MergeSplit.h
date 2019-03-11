// generating k-mers of fixed length + merging k-mers into longer fragments + preparing for spase dp
#ifndef MERGE_SPLIT_H_
#define MERGE_SPLIT_H_

#include <iostream> //std::cout 
#include <fstream>   
#include <cstdlib>   // std::labs, std::EXIT FAILURE, std::EXIT SUCCESS
#include <cstdio>    // std::FILE std::perror
#include <vector>
#include <set>
#include "Sorting.h"
#include "Clustering.h"
#include "Types.h"

using namespace std;


void MergeClusters (Options & smallOpts, vector<Cluster> & refinedClusters, vector<Cluster> & vt, int r) {
 
    // add merge split code here
    //------------------------------------------------------------------
    // Step 1: merge matches with diagnol difference smaller than maxDiag
    //------------------------------------------------------------------
    vector<Cluster> v;
    vector<Cluster> vq;
    Options mergeOpts = smallOpts;
    mergeOpts.maxGap = 3500;
    mergeOpts.maxDiag = 20;

    StoreDiagonalClusters(refinedClusters[r].matches, v, mergeOpts, true); 


    if (v.size() != 0) {
        std::set<unsigned int> s1;     // s1 stores x boundary; s2 stores y boundary      
        std::set<unsigned int> s2;    //elements in s are arranged in strictly increasing order. and no duplicates
    
        /*
        // debug
        if (r == 2) {
        cout << "refinedClusters[r].matches.size(): " << refinedClusters[r].matches.size() << endl;
        }
        */

        /*
        // Debug code ----------- print out "orignal"
    
        seqan::SeedSet<IndSeed, seqan::Unordered> seedSet1;
        for (unsigned int i = 0; i != v.size(); ++i) {
            seqan::addSeed(seedSet1, IndSeed(v[i].qStart, v[i].tStart, v[i].qEnd, v[i].tEnd, i), seqan::Single());
        } 
        for (TIterator tt = begin(seedSet1, seqan::Standard()); tt != end(seedSet1, seqan::Standard()); ++tt) {
            cerr << *tt << endl;
        }
    
        
        cout << "v.szie(): " << v.size() << endl;
        const string filename7("/home/cmb-16/mjc/jingwenr/lra/lra_test/TEST/Merge.txt");  
        FILE *fr = fopen(filename7.c_str(), "w");
        SaveseedSet (seedSet1, fr); 
        fclose(fr);
        */  
    

        //---------------------------------------------------------------------------------
        // Step 2: sort matches by x first and then by y coordinates in each Cluster v[i](merged matches)
        //---------------------------------------------------------------------------------
        for (unsigned int i = 0; i < v.size(); ++i) {

            //insert x boundaries into s1 
            s1.insert(v[i].qStart);
            s1.insert(v[i].qEnd);

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
                vq.push_back(Cluster(v[i].start, v[i].end));
            }
            else{
                while (ii != v[i].end && l != qe) { // end is end + 1!!!!!!
                    while (refinedClusters[r].matches[ii].first.pos <= *l && ii != v[i].end) {
                        if (refinedClusters[r].matches[ii].first.pos + smallOpts.globalK >= *l) {
                            // split 
                            if (jj + 1 <= ii) {
                                assert(jj < refinedClusters[r].matches.size());
                                vq.push_back(Cluster(jj, ii));
                            }
                            vq.push_back(Cluster(ii, ii + 1));
                            jj = ii + 1;
                        }
                        ++ii;
                    }
                    assert(ii - 1 < refinedClusters[r].matches.size());
                    if (ii > jj && refinedClusters[r].matches[ii - 1].first.pos + smallOpts.globalK < *l) {
                        vq.push_back(Cluster(jj, ii));
                        jj = ii;
                    }
                    std::advance(l, 1);
                }   
                if (ii != v[i].end && l == qe) {
                    assert(i < v.size());
                    vq.push_back(Cluster(ii, v[i].end));
                }
            }
        }


        /*
        cerr << "length(seedSet2): " << length(seedSet2)<< endl;;
        for (TIterator tt = begin(seedSet2, seqan::Standard()); tt != end(seedSet2, seqan::Standard()); ++tt) {
            cerr << *tt << endl;
        }
        */
        /*
        //cout << "vq.szie(): " << vq.size() << endl;
        const string filename2("/home/cmb-16/mjc/jingwenr/lra/lra_test/TEST/seedSet1.txt");  
        FILE *fm = fopen(filename2.c_str(), "w");
        SaveseedSet (seedSet2, fm); 
        fclose(fm); 
        */      


        //----------------------------------------------------------------------------
        // Step 4: sort matches by y first and then by x coordinates in each Cluster vq[i](merged matches)  && split based on y boundaries
        //----------------------------------------------------------------------------
        for (unsigned int ij = 0; ij < vq.size(); ++ij) {
            int css = vq[ij].start; 
            int cee = vq[ij].start;
            cee = css + 1;
            GenomePos qStart = refinedClusters[r].matches[css].first.pos, 
            qEnd = refinedClusters[r].matches[css].first.pos + smallOpts.globalK, 
            tStart = refinedClusters[r].matches[css].second.pos, 
            tEnd = refinedClusters[r].matches[css].second.pos + smallOpts.globalK;  

            while(cee < vq[ij].end) {
                qStart = min(qStart, refinedClusters[r].matches[cee].first.pos);
                qEnd = max(qEnd, refinedClusters[r].matches[cee].first.pos + smallOpts.globalK);
                tStart = min(tStart, refinedClusters[r].matches[cee].second.pos);
                tEnd = max(tEnd, refinedClusters[r].matches[cee].second.pos + smallOpts.globalK);
                ++cee;
            }
            vq[ij].qStart = qStart;
            vq[ij].qEnd = qEnd;
            vq[ij].tStart = tStart;
            vq[ij].tEnd = tEnd;

            //insert y boundaries into s2
            s2.insert(vq[ij].tStart);
            s2.insert(vq[ij].tEnd);

            // sort matches by y and then by x inside each vq[i]
            vector<pair<GenomeTuple, GenomeTuple>>::iterator begin = refinedClusters[r].matches.begin();
            if (vq[ij].start != 0) {
                std::advance(begin , vq[ij].start);
            }    

            vector<pair<GenomeTuple, GenomeTuple>>::iterator end = refinedClusters[r].matches.begin();
            if (vq[ij].end != 0) {
                std::advance(end, vq[ij].end);
            }
            CartesianTargetSort<GenomeTuple>(begin, end);            
        }


        
        v.clear();
        for (unsigned int i = 0; i < vq.size(); ++i) {
            std::set<unsigned int>::iterator ts, te; 
            ts = s2.lower_bound(vq[i].tStart);
            te = s2.lower_bound(vq[i].tEnd);
            std::set<unsigned int>::iterator l = ts;
            std::advance(l, 1);
            unsigned int ii = vq[i].start;
            unsigned int jj = ii;

            /*
            // debug code
            cout << "ts: " << *ts << endl;
            cout << "te: " << *te << endl;
            cout << "l: " << *l << endl;
            */
            
            
            if (*l == *ts) {
                //debug code

                //cerr << "*l == *ts, push back" << "Cluster(" << vq[i].start << ", " << vq[i].end << ", " << vq[i].qStart << ", " << vq[i].qEnd << ", " << vq[i].tStart << ", " << vq[i].tEnd << ", 0)" << endl; 
                vt.push_back(Cluster(vq[i].start, vq[i].end));
            }
            else {
                while (ii != vq[i].end && l != te) { 
                    while (refinedClusters[r].matches[ii].second.pos <= *l && ii != vq[i].end) {
                        if (refinedClusters[r].matches[ii].second.pos + smallOpts.globalK >= *l) {
                            // split 
                            if (jj + 1 <= ii) {
                                assert(jj < ii);
                                vt.push_back(Cluster(jj, ii));
                            }
                            vt.push_back(Cluster(ii, ii + 1));
                            jj = ii + 1;
                        }
                        ++ii;
                    }
                    if (ii > jj && refinedClusters[r].matches[ii - 1].second.pos + smallOpts.globalK < *l) {
                        assert(jj < ii);
                        vt.push_back(Cluster(jj, ii));
                        jj = ii;
                    }
                    std::advance(l, 1);
                }
                if (ii != vq[i].end && l == te) {
                    assert(ii < vq[i].end);
                    vt.push_back(Cluster(ii, vq[i].end));
                }
            }
        }
    

        vq.clear();
        for (unsigned int ij = 0; ij < vt.size(); ++ij) {
            int css = vt[ij].start; 
            int cee = vt[ij].start;
            cee = css + 1;
            GenomePos qStart = refinedClusters[r].matches[css].first.pos, 
            qEnd = refinedClusters[r].matches[css].first.pos + smallOpts.globalK, 
            tStart = refinedClusters[r].matches[css].second.pos, 
            tEnd = refinedClusters[r].matches[css].second.pos + smallOpts.globalK;  

            while(cee < vt[ij].end) {
                qStart = min(qStart, refinedClusters[r].matches[cee].first.pos);
                qEnd = max(qEnd, refinedClusters[r].matches[cee].first.pos + smallOpts.globalK);
                tStart = min(tStart, refinedClusters[r].matches[cee].second.pos);
                tEnd = max(tEnd, refinedClusters[r].matches[cee].second.pos + smallOpts.globalK);
                ++cee;
            }
            vt[ij].qStart = qStart;
            vt[ij].qEnd = qEnd;
            vt[ij].tStart = tStart;
            vt[ij].tEnd = tEnd;
            assert(vt[ij].qStart < vt[ij].qEnd);
            assert(vt[ij].tStart < vt[ij].tEnd);
        }
    }
}

#endif
