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
#include <set>
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
	int length;
	if ((ReadEnd - ReadStart + 1)%10 != 0) {
		length = (ReadEnd - ReadStart + 1)/10 + 1;
	}
	else {
		length = (ReadEnd - ReadStart + 1)/10;
	}
	std::vector<int> ReadWindow(length, 0);
	
	int j = 0;
	for (int i = 0; i < ReadWindow.size(); ++i) { // The region that ReadWindow[i] represents is [ReadStart + i*10, ReadStart + (i+1)*10)
		while (j < matches.size() and matches[j].first.pos >= ReadStart + i*10 and matches[j].first.pos < ReadStart + (i+1)*10) {
			++ReadWindow[i];
			++j;
		}
	}

	for (int i = 0; i < ReadWindow.size() - 1; ++i) {
		ReadWindow[i] = ReadWindow[i] + ReadWindow[i + 1];
	}


	cerr << "ReadWindow: " << ReadWindow.size() << endl;

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
	CartesianTargetSort<GenomeTuple>(matches.begin(), matches.end());     
	int GenomeStart = matches[0].second.pos;
	int GenomeEnd = matches[ch- 1].second.pos;
	if ((GenomeEnd - GenomeStart + 1)%10 != 0) {
		length = (GenomeEnd - GenomeStart + 1)/10 + 1;
	}
	else {
		length = (GenomeEnd - GenomeStart + 1)/10;
	}
	std::vector<int> GenomeWindow(length, 0);

	j = 0;
	for (int i = 0; i < GenomeWindow.size(); ++i) {
		while (j < matches.size() and matches[j].second.pos >= GenomeStart+ i*10 and matches[j].second.pos < GenomeStart + (i+1)*10) {
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


	//---------find the boundary in Read ------------
	double mean = (std::accumulate(ReadWindow.begin(), ReadWindow.end(), 0))/ReadWindow.size();
	double sqtsum = 0;
	for (int i = 0; i < ReadWindow.size(); ++i) {
		sqtsum += std::pow((double)ReadWindow[i] - mean, 2); 
	}
	double stdev = std::sqrt(sqtsum/(ReadWindow.size()-1));	
	double threshold = mean + 4.5*stdev;

	std::vector<int> ReadIndex; 
	for (int i = 0; i < ReadWindow.size(); ++i) {
		if (ReadWindow[i] >= threshold) ReadIndex.push_back(i);
	}

	int ReadBoundaryStart, ReadBoundaryEnd;
	if (!ReadIndex.empty()) {
		ReadBoundaryStart = (*ReadIndex.begin())*10 + ReadStart;
		ReadBoundaryEnd = (ReadIndex.back())*10 + 2*10 + ReadStart;		
	}
	else {
		ReadBoundaryStart = ReadEnd;
		ReadBoundaryEnd = ReadEnd;
	}


	//cerr << "ReadIndex: " << ReadIndex << endl;
	//cerr << "ReadBoundaryStart: " << ReadBoundaryStart << " ReadBoundaryEnd: " << ReadBoundaryEnd << endl;


	//---------find the boundary in Genome ------------
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
		GenomeBoundaryStart = (*GenomeIndex.begin())*10 + GenomeStart;
		GenomeBoundaryEnd = (GenomeIndex.back())*10 + 2*10 + GenomeStart;		
	}
	else {
		GenomeBoundaryStart = GenomeEnd;
		GenomeBoundaryEnd = GenomeEnd;	
	}

	//cerr << "GenomeIndex: " << GenomeIndex << endl;
	//cerr << "GenomeBoundaryStart: " << GenomeBoundaryStart << "  , GenomeBoundaryEnd: " << GenomeBoundaryEnd << endl;


	// -------- Merge ---------

    Options mergeOpts = opts;
    mergeOpts.maxGap = 5000;
    mergeOpts.maxDiag = 20;
	int cs = 0;
	int ce = 0;

	while (cs < matches.size()) {
		ce = cs + 1;
		GenomePos qStart=matches[cs].first.pos, 
			qEnd=matches[cs].first.pos + opts.globalK, 
			tStart=matches[cs].second.pos, 
			tEnd=matches[cs].second.pos + opts.globalK;

		while (ce < matches.size() 
					and (matches[ce].first.pos < ReadBoundaryStart or matches[ce].first.pos >= ReadBoundaryEnd)
					and (matches[ce - 1].first.pos < ReadBoundaryStart or matches[ce - 1].first.pos >= ReadBoundaryEnd)
					and (matches[ce].second.pos < GenomeBoundaryStart or matches[ce].second.pos >= GenomeBoundaryEnd)
					and (matches[ce - 1].second.pos < GenomeBoundaryStart or matches[ce - 1].second.pos >= GenomeBoundaryEnd)) {

			int rgap = matches[ce].first.pos - matches[ce - 1].first.pos - opts.globalK;
			int ggap = matches[ce].second.pos - matches[ce - 1].second.pos - opts.globalK;

			if (gapdifference(ce, matches) <= mergeOpts.maxDiag and std::max(std::abs(rgap), std::abs(ggap)) <= mergeOpts.maxGap) {
				qStart = min(qStart, matches[ce].first.pos);
				qEnd   = max(qEnd, matches[ce].first.pos + opts.globalK);
				tStart = min(tStart, matches[ce].second.pos);
				tEnd   = max(tEnd, matches[ce].second.pos + opts.globalK);				
				++ce;
			}
			else {
				v.push_back(Cluster(cs, ce, qStart, qEnd, tStart, tEnd, 0)); // I just put 0 here for the strand argument
				cs = ce;
				break;
			}
		}
		if (ce < matches.size() and (!(matches[ce].first.pos < ReadBoundaryStart or matches[ce].first.pos >= ReadBoundaryEnd) 
				or !(matches[ce - 1].first.pos < ReadBoundaryStart or matches[ce - 1].first.pos >= ReadBoundaryEnd)
				or !(matches[ce].second.pos < GenomeBoundaryStart or matches[ce].second.pos >= GenomeBoundaryEnd)
				or !(matches[ce - 1].second.pos < GenomeBoundaryStart or matches[ce - 1].second.pos >= GenomeBoundaryEnd))) {
			v.push_back(Cluster(cs, ce, qStart, qEnd, tStart, tEnd, 0));
			cs = ce;
		}

		if (ce == matches.size()) {
			v.push_back(Cluster(cs, ce, qStart, qEnd, tStart, tEnd, 0));
			cs = ce;
		}

	}
}


# endif