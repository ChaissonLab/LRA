// Merge fragments 
#ifndef MERGE_H_
#define MERGE_H_

#include <iostream> //std::cout 
#include <fstream>   
#include <cstdlib>   // std::labs, std::EXIT FAILURE, std::EXIT SUCCESS
#include <string> 
#include <cmath>        // std::labs
#include <cstdio>    // std::FILE std::perror
#include <vector>
#include <numeric>
#include <utility>
#include <set>
#include <list>
#include "Sorting.h"
#include "Clustering.h"
#include "Options.h"


int 
gapdifference (int & ce, GenomePairs & matches) {
	int prevDiag = matches[ce-1].second.pos - matches[ce-1].first.pos;
	int curDiag = matches[ce].second.pos - matches[ce].first.pos;
	return std::abs(prevDiag - curDiag);
}



void 
mergeClusters (Options & opts, GenomePairs & matches, std::vector<Cluster> & v, int r, string & baseName) {
	// ---------- sort fragments(starting points) first by first.pos, then by second.pos ---------------
	// store the number of starting points in every window in ReadWindow
	CartesianSort<GenomeTuple>(matches.begin(), matches.end());     
	int ReadStart = matches[0].first.pos;
	int ch = matches.size();
	int ReadEnd = matches[ch - 1].first.pos;        
	int length, ind;
	if ((ReadEnd - ReadStart + 1)%50 != 0) {
		length = (ReadEnd - ReadStart + 1)/50 + 1;
		ind = 1;
	}
	else {
		length = (ReadEnd - ReadStart + 1)/50;
		ind = 0;
	}
	std::vector<int> ReadWindow(length, 0);
	
	int j = 0;
	for (int i = 0; i < ReadWindow.size(); ++i) { // The region that ReadWindow[i] represents is [ReadStart + i*50, ReadStart + (i+1)*50)
		while (j < matches.size() and matches[j].first.pos >= ReadStart + i*50 and matches[j].first.pos < ReadStart + (i+1)*50) {
			++ReadWindow[i];
			++j;
		}
	}

	for (int i = 0; i < ReadWindow.size() - 1; ++i) {
		ReadWindow[i] = ReadWindow[i] + ReadWindow[i + 1];
	}

	if (opts.dotPlot) {
		stringstream outNameStrm;
		outNameStrm << baseName + "." << r << ".ReadWindow.txt";
		ofstream baseDots(outNameStrm.str().c_str());
		for (int m=0; m < ReadWindow.size(); m++) {
			baseDots << ReadWindow[m] << "\t" << endl;
		}
		baseDots.close();
	}



	// ----------- sort fragments(starting points) first by second.pos, and then by first.pos ------------
	// store the number of starting points in every window in GenomeWindow
/*
	CartesianTargetSort<GenomeTuple>(matches.begin(), matches.end());     
	int GenomeStart = matches[0].second.pos;
	int GenomeEnd = matches[ch- 1].second.pos;
	if ((GenomeEnd - GenomeStart + 1)%50 != 0) {
		length = (GenomeEnd - GenomeStart + 1)/50 + 1;
	}
	else {
		length = (GenomeEnd - GenomeStart + 1)/50;
	}
	std::vector<int> GenomeWindow(length, 0);

	j = 0;
	for (int i = 0; i < GenomeWindow.size(); ++i) {
		while (j < matches.size() and matches[j].second.pos >= GenomeStart+ i*50 and matches[j].second.pos < GenomeStart + (i+1)*50) {
			++GenomeWindow[i];
			++j;
		}
	}
	for (int i = 0; i < GenomeWindow.size() - 1; ++i) {
		GenomeWindow[i] = GenomeWindow[i] + GenomeWindow[i + 1];
	}

	cerr << "GenomeWindow: " << GenomeWindow.size() << endl;
	cerr << "GenomeStart: " << GenomeStart << endl;

	if (opts.dotPlot) {
		stringstream outNameStrm;
		outNameStrm << baseName + "." << r << ".GenomeWindow.txt";
		ofstream baseDots(outNameStrm.str().c_str());
		for (int m=0; m < GenomeWindow.size(); m++) {
			baseDots << GenomeWindow[m] << "\t" << endl;
		}
		baseDots.close();
	}
*/

	//---------find the threshold in ReadWindow counts and pick ones that higher than the threshold ------------
	double mean = (std::accumulate(ReadWindow.begin(), ReadWindow.end(), 0))/ReadWindow.size();
	double sqtsum = 0;
	for (int i = 0; i < ReadWindow.size(); ++i) {
		sqtsum += std::pow((double)ReadWindow[i] - mean, 2); 
	}
	double stdev = std::sqrt(sqtsum/(ReadWindow.size()-1));	
	double threshold = mean + 3*stdev;

	std::vector<int> ReadIndex(ReadWindow.size(), 0);
	for (int i = 0; i < ReadWindow.size(); ++i) {
		if (ReadWindow[i] >= threshold) ReadIndex[i] = 1;
	}


/*
	int ReadBoundaryStart, ReadBoundaryEnd;
	if (!ReadIndex.empty()) {
		ReadBoundaryStart = (*ReadIndex.begin())*50 + ReadStart;
		ReadBoundaryEnd = (ReadIndex.back())*50 + 2*50 + ReadStart;		
	}
	else {
		ReadBoundaryStart = ReadEnd;
		ReadBoundaryEnd = ReadEnd;
	}
*/

	//cerr << "ReadIndex: " << ReadIndex << endl;
	//cerr << "ReadBoundaryStart: " << ReadBoundaryStart << " ReadBoundaryEnd: " << ReadBoundaryEnd << endl;


	//---------find the boundary in Genome ------------
/*
	mean = (std::accumulate(GenomeWindow.begin(), GenomeWindow.end(), 0))/GenomeWindow.size();
	sqtsum = 0;
	for (int i = 0; i < GenomeWindow.size(); ++i) {
		sqtsum += std::pow((double)GenomeWindow[i] - mean, 2); 
	}
	stdev = std::sqrt(sqtsum/(GenomeWindow.size()-1));	
	threshold = mean + 4.5*stdev;

	std::vector<int> GenomeIndex; 
	for (int i = 0; i < GenomeWindow.size(); ++i) {
		if (GenomeWindow[i] >= threshold) GenomeIndex.push_back(i);
	}
	
	int GenomeBoundaryStart, GenomeBoundaryEnd;
	if (!GenomeIndex.empty()) {
		GenomeBoundaryStart = (*GenomeIndex.begin())*50 + GenomeStart;
		GenomeBoundaryEnd = (GenomeIndex.back())*50 + 2*50 + GenomeStart;		
	}
	else {
		GenomeBoundaryStart = GenomeEnd;
		GenomeBoundaryEnd = GenomeEnd;	
	}

	//cerr << "GenomeIndex: " << GenomeIndex << endl;
	//cerr << "GenomeBoundaryStart: " << GenomeBoundaryStart << "  , GenomeBoundaryEnd: " << GenomeBoundaryEnd << endl;
*/


	// -------- Merge ---------
/*
	int start = 0;
	int end = 0;	
	ts = 0; te = 0;
	for (int j = 0; j < Boundary.size(); ++j) {

		start = Boundary[j].first*50 + ReadStart;
		end = (Boundary[j].second)*50 + ReadStart;

		while (ts < matches.size() and matches[ts].first.pos >= start and matches[ts].first.pos < end) {
			te = ts + 1;
			GenomePos qStart=matches[ts].first.pos, 
					  qEnd=matches[ts].first.pos + opts.globalK, 
					  tStart=matches[ts].second.pos, 
				      tEnd=matches[ts].second.pos + opts.globalK;		

			while (te < matches.size() and matches[te].first.pos >= start and matches[te].first.pos < end) {

				int rgap = matches[te].first.pos - matches[te - 1].first.pos - opts.globalK;
				int ggap = matches[te].second.pos - matches[te - 1].second.pos - opts.globalK;	

				if (gapdifference(te, matches) <= mergeOpts.maxDiag and std::max(std::abs(rgap), std::abs(ggap)) <= mergeOpts.maxGap) {
					qStart = min(qStart, matches[te].first.pos);
					qEnd   = max(qEnd, matches[te].first.pos + opts.globalK);
					tStart = min(tStart, matches[te].second.pos);
					tEnd   = max(tEnd, matches[te].second.pos + opts.globalK);				
					++te;
				}
				else {
					v.push_back(Cluster(ts, te, qStart, qEnd, tStart, tEnd, 0)); // I just put 0 here for the strand argument
					ts = te;
					break;
				}					
			}	

			if (te < matches.size() and matches[te].first.pos >= end) {
				v.push_back(Cluster(ts, te, qStart, qEnd, tStart, tEnd, 0));
				ts = te;				
			}

			if (te == matches.size()) {
				v.push_back(Cluster(ts, te, qStart, qEnd, tStart, tEnd, 0));
				ts = te;
			}
		}
		if (matches[ts].first.pos >= end) {
			continue;
		}
		if (ts == matches.size()) {
			continue;
		}

	}
*/




	// -----------------------Merge--------------------------
    Options mergeOpts = opts;
    mergeOpts.maxGap = 3000;
    mergeOpts.maxDiag = 40;
	int cs= 0;
	int ce = 0;
	std::vector<std::pair<int, int>> Boundary;
	std::vector<int> index; // TODO(Jingwen): this is only for debug. delete this later
	// set up the boundary according ReadIndex
	while (cs < ReadIndex.size()) {
		ce = cs + 1;

		while (ce < ReadIndex.size() and ReadIndex[ce - 1] == ReadIndex[ce]) {
			++ce;
		}
		if (ce < ReadIndex.size() and ReadIndex[ce - 1] != ReadIndex[ce]) {
			std::pair<int, int> s (cs, ce);
			Boundary.push_back(s);
			int t = ReadIndex[cs];
			index.push_back(t);
			cs = ce;
		}
		if (ce == ReadIndex.size()) {
			std::pair<int, int> s (cs, ce);
			Boundary.push_back(s);
			int t = ReadIndex[cs];
			index.push_back(t);
			cs = ce;
		}
	}

	// TODO(Jingwen): for debug 
	if (opts.dotPlot) {
		stringstream outNameStrm;
		outNameStrm << baseName + "." << r << ".Boundary.txt";
		ofstream baseDots(outNameStrm.str().c_str());
		for (int m=0; m < Boundary.size(); m++) {
			baseDots << ReadStart + Boundary[m].first*50 << "\t" << ReadStart + Boundary[m].second*50 << "\t" << index[m] << endl;
		}
		baseDots.close();
	}



	int ts = 0;
	int te = 0;
	// merge without any boundary
	if (Boundary.empty()) {
		while (ts < matches.size()) {
			te = ts + 1;
			GenomePos qStart  =matches[ts].first.pos, 
					  qEnd = matches[ts].first.pos + opts.globalK, 
					  tStart = matches[ts].second.pos, 
				      tEnd = matches[ts].second.pos + opts.globalK;	
			while (te < matches.size()) {
				int rgap = matches[te].first.pos - matches[te - 1].first.pos - opts.globalK;
				int ggap = matches[te].second.pos - matches[te - 1].second.pos - opts.globalK;	
				if (gapdifference(te, matches) <= mergeOpts.maxDiag and std::max(std::abs(rgap), std::abs(ggap)) <= mergeOpts.maxGap) {
					qStart = min(qStart, matches[te].first.pos);
					qEnd   = max(qEnd, matches[te].first.pos + opts.globalK);
					tStart = min(tStart, matches[te].second.pos);
					tEnd   = max(tEnd, matches[te].second.pos + opts.globalK);				
					++te;
				}
				else {
					v.push_back(Cluster(ts, te, qStart, qEnd, tStart, tEnd, 0)); // I just put 0 here for the strand argument
					ts = te;
					break;
				}	
			}	
			if (te == matches.size()) {
				v.push_back(Cluster(ts, te, qStart, qEnd, tStart, tEnd, 0));
				ts = te;
			}
		}
	}

	//Boundary is not empty
	int start = 0;
	int end = 0;	
	ts = 0; te = 0;
	for (int j = 0; j < Boundary.size(); ++j) {

		start = Boundary[j].first*50 + ReadStart;
		end = (Boundary[j].second)*50 + ReadStart;

		if (index[j] == 0) {
			// merge anchors in [start, end)

			while (ts < matches.size() and matches[ts].first.pos >= start and matches[ts].first.pos < end) {
				te = ts + 1;
				GenomePos qStart=matches[ts].first.pos, 
						  qEnd=matches[ts].first.pos + opts.globalK, 
						  tStart=matches[ts].second.pos, 
					      tEnd=matches[ts].second.pos + opts.globalK;		

				while (te < matches.size() and matches[te].first.pos >= start and matches[te].first.pos < end) {

					int rgap = matches[te].first.pos - matches[te - 1].first.pos - opts.globalK;
					int ggap = matches[te].second.pos - matches[te - 1].second.pos - opts.globalK;	

					if (gapdifference(te, matches) <= mergeOpts.maxDiag and std::max(std::abs(rgap), std::abs(ggap)) <= mergeOpts.maxGap) {
						qStart = min(qStart, matches[te].first.pos);
						qEnd   = max(qEnd, matches[te].first.pos + opts.globalK);
						tStart = min(tStart, matches[te].second.pos);
						tEnd   = max(tEnd, matches[te].second.pos + opts.globalK);				
						++te;
					}
					else {
						v.push_back(Cluster(ts, te, qStart, qEnd, tStart, tEnd, 0)); // I just put 0 here for the strand argument
						ts = te;
						break;
					}					
				}	

				if (te < matches.size() and matches[te].first.pos >= end) {
					v.push_back(Cluster(ts, te, qStart, qEnd, tStart, tEnd, 0));
					ts = te;				
				}

				if (te == matches.size()) {
					v.push_back(Cluster(ts, te, qStart, qEnd, tStart, tEnd, 0));
					ts = te;
				}
			}
			if (matches[ts].first.pos >= end) {
				continue;
			}
			if (ts == matches.size()) {
				continue;
			}
		}
		else {
			// Do not merge anchors in [start, end)
			while (ts < matches.size() and matches[ts].first.pos >= start and matches[ts].first.pos < end) {
				GenomePos qStart = matches[ts].first.pos, 
						  qEnd = matches[ts].first.pos + opts.globalK, 
						  tStart = matches[ts].second.pos, 
					      tEnd = matches[ts].second.pos + opts.globalK;	
				v.push_back(Cluster(ts, ts + 1, qStart, qEnd, tStart, tEnd, 0));				 
				++ts;
			}
		}

	}

/*
	stringstream outNameStrm;
	outNameStrm << "AnchorEfficiency.txt";
	ofstream baseDots;
	baseDots.open(outNameStrm.str().c_str(), std::ios::app);
	baseDots << v.size() << "\t" << matches.size() << endl;
	baseDots.close();
*/
	
}


# endif