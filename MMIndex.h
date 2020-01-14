#ifndef MMINDEX_H_
#define MMINDEX_H_
#include "TupleOps.h"
#include "Options.h"
#include "Genome.h"
#include "Sorting.h"
#include "MinCount.h"
#include <numeric>     /*iota */
#include <cmath>       /* ceil */
#include <utility>      // std::pair, std::make_pair

template<typename Tup>
class SortByPos {
 public:
	int operator() (const Tup &a, const Tup &b) {
		return (a.pos < b.pos);
	}
};


template<typename Tup>
void PrintPairs(vector<pair<Tup, Tup> > &mins, int k, int cluster=-1) {
	CartesianTargetSort<Tup>(mins);
	for(int i = 0; i < mins.size();i++) {
		string s;
		TupleToString(mins[i].first.t, k, s);
#ifdef _TESTING_
		if (cluster != -1) {
			cout << "clust\t" << cluster << "\t";
		}
		cout << i << "\t" << mins[i].first.pos << "\t" << mins[i].second.pos << "\t" << s << endl;
#endif
	}	
}
	

template<typename Tup>
void PrintIndex(vector<Tup> &minimizers, int k) {
	sort(minimizers.begin(), minimizers.end(), SortByPos<Tup>());
	for(int i = 0; i < minimizers.size();i++) {
		string s;
		TupleToString(minimizers[i].t, k, s);
		cout << i << "\t" << minimizers[i].pos << "\t" << s << endl;
	}	
}


template<typename Tup>
void RemoveFrequent(vector<Tup> &minimizers, int maxFreq) {
	int c=0,n=0;
	int before=minimizers.size();
	while(n < minimizers.size()) {
		int ne=n;
		while (ne < minimizers.size() && 
					 minimizers[ne].t == minimizers[n].t) { ne++;}
		if (ne - n < maxFreq) {
			int end = ne;
			for (ne = n; ne < end; ne++, c++) {
				minimizers[c] = minimizers[ne];
			}
		}
		n=ne;
	}
	minimizers.resize(c);
}


template<typename Tup>
void RemoveFrequent(vector<Tup> &minimizers, vector<bool> &remove) {
	int c = 0;

	for (int n = 0; n < minimizers.size(); n++) {
		if (remove[n] == 0) {
			minimizers[c] = minimizers[n];
			c++;
		}	
	}
	minimizers.resize(c);
}


template<typename Tup>
void CompressFrequent(vector<Tup> &minimizers) {
	int c=0,n=0;
	int before=minimizers.size();
	while(n < minimizers.size()) {		
		while (n < minimizers.size() && 
					 minimizers[n].t == minimizers[c].t) { n++;}
		
		minimizers[c] = minimizers[n-1];
		c++;
	}
	minimizers.resize(c);
}


class LocalIndex {
 public:
	int localIndexSize;
	int k;
	int w;
	int maxFreq;
	vector<LocalTuple>  minimizers;
	vector<uint64_t>    seqOffsets; // seqOffsets stores actual boundaries 
	vector<uint64_t>    tupleBoundaries; // tupleBoundaries stores the number of minimizers in the corresponding interval
	uint64_t offset;
	
	LocalIndex() { 
		k = 10;
		w=5; 
		offset=0;
		maxFreq=5;
		tupleBoundaries.push_back(0);
		seqOffsets.push_back(0);
		localIndexSize = 1 << (LOCAL_POS_BITS-1);
	}
	
	LocalIndex( LocalIndex &init) { 
		k=init.k;
		w=init.w;
		offset=0;
		maxFreq=init.maxFreq;
		localIndexSize = 1 << (LOCAL_POS_BITS-1);
		tupleBoundaries.push_back(0);
		seqOffsets.push_back(0);
	}		

	void Write(string filename) {
		ofstream fout(filename.c_str(), ios::out|ios::binary);
		fout.write((char*)&k, sizeof(int));
		fout.write((char*)&w, sizeof(int));
		int nRegions=seqOffsets.size();
		fout.write((char*)&nRegions, sizeof(int));
		fout.write((char*)&seqOffsets[0], sizeof(uint64_t)*seqOffsets.size());
		fout.write((char*)&tupleBoundaries[0], sizeof(uint64_t)*tupleBoundaries.size());
		uint64_t nMin = minimizers.size();
		fout.write((char*)&nMin, sizeof(uint64_t));
		fout.write((char*)&minimizers[0], sizeof(LocalTuple)*minimizers.size());
		fout.close();
	}

	int Read(string filename) {
		ifstream fin(filename.c_str(), ios::in|ios::binary);
		if (fin.good() == false or fin.eof() == true) {
			return 0;
		}
		fin.read((char*)&k, sizeof(int));
		fin.read((char*)&w, sizeof(int));
		int nRegions;
		fin.read((char*)&nRegions, sizeof(int));
		seqOffsets.resize(nRegions);		
		fin.read((char*)&seqOffsets[0], sizeof(uint64_t)*nRegions);
		tupleBoundaries.resize(nRegions);
		fin.read((char*)&tupleBoundaries[0], sizeof(uint64_t)*nRegions);
		uint64_t nMin;
		fin.read((char*) &nMin, sizeof(uint64_t));
		minimizers.resize(nMin);
		fin.read((char*)&minimizers[0], sizeof(LocalTuple)*nMin);
		fin.close();
		return 1;
	}
	
	int LookupIndex(uint64_t querySeqPos) {
		if (seqOffsets.size() == 0) {
			return 0;
		}
		assert(querySeqPos <= seqOffsets[seqOffsets.size()-1]);
		vector<uint64_t>::iterator it;
		it = lower_bound(seqOffsets.begin(), seqOffsets.end(), querySeqPos);
		//		while(it != seqOffsets.end() and *it == querySeqPos) { ++it;}
		int index = it - seqOffsets.begin();
		if (*it != querySeqPos) {
			return index - 1;
		}
		else {
			return index;
		}
	}

	void MinimizerBounds(uint64_t querySeqPos, uint64_t &lb, uint64_t &ub) {
		assert(querySeqPos < minimizers.size());
		int index = this->LookupIndex(querySeqPos);
		assert(index < tupleBoundaries.size());
		lb = tupleBoundaries[index];
		ub = tupleBoundaries[index+1];
	}

	void IndexSeq(char* seq, int seqLen) {
		int gi = 0;
		int nIndex = seqLen / localIndexSize;

		if (seqLen % localIndexSize != 0) {
			nIndex +=1;
		}
		GenomePos seqPos=0;

		vector<LocalTuple> locMinimizers;
		GenomePos netSize=0;
		
		for (int i = 0; i < nIndex; i++) {	
			locMinimizers.clear();
			StoreMinimizers<LocalTuple, SmallTuple>(&seq[seqPos], min((GenomePos) seqLen, (GenomePos) seqPos+localIndexSize) - seqPos, k, w, locMinimizers, false); // canonical = false
			//
			// Sort minimzers by tuple value.
			//
			sort(locMinimizers.begin(), locMinimizers.end());
			//
			// Remove frequenct tuples
			//
			RemoveFrequent(locMinimizers, maxFreq);

			//
			// Update local sequence pos (index in chrom).
			//
			seqPos+=(GenomePos)min((int)localIndexSize, (int) (seqLen - seqPos));

			//
			// Add boundaries representing the end of the current interval.
			//
			seqOffsets.push_back(offset+seqPos);

			//
			// Add minimizers and store where they end.
			//
			minimizers.insert(minimizers.end(), locMinimizers.begin(), locMinimizers.end());
			tupleBoundaries.push_back(minimizers.size());
			netSize+=minimizers.size();
		}
		//
		// Update offset for recently added sequence
		//
		offset+=seqLen;
	}

	void IndexFile(string &genome) {	
		gzFile f = gzopen(genome.c_str(), "r");
		kseq_t *ks = kseq_init(f);
		while (kseq_read(ks) >= 0) { 
			cerr << "Storing for "<< ks->name.s << endl;
			IndexSeq(ks->seq.s, ks->seq.l);
		}
	}

};


class Freminimizers {
public:
	vector<std::pair<uint32_t, uint32_t>> *freCount;
	vector<uint32_t> index;

	// Initialization
	Freminimizers (vector<std::pair<uint32_t, uint32_t>> *f) :freCount(f){
		index.resize(freCount->size());
		std::iota(index.begin(), index.end(), 0);
		cerr << "3" << endl;
		Sort();
		cerr << "4" << endl;
	}

	bool operator () (uint32_t i, uint32_t j) {
		assert(i < freCount->size());
		assert(j < freCount->size());
		return (*freCount)[i].first < (*freCount)[j].first; // TODO(Jingwen): "<" or "<="??
	}

	void Sort() {
		std::sort(index.begin(), index.end(), *this);
	}

	uint32_t & operator[] (const uint32_t &i) {
		return index[i];
	}

};

void CountSort(const vector<uint32_t> & Freq, const int & RANGE, const vector<bool> & Remove, vector<uint32_t> & Sortindex){
//void CountSort (const vector<std::pair<uint32_t, uint32_t>> & v, const int & RANGE, 
//						vector<std::pair<uint32_t, uint32_t>> & sortv) {

	// Create a count vector to store counts of each frequency
	vector<uint32_t> count(RANGE + 1, 0);
  
    // Store counts of each frequency in v
    for (uint32_t i = 0; i < Freq.size(); i++) {
    	if (Remove[i] == 0) {
    		++count[Freq[i]];
    	}
    }
  
    // Change count[i] so that count[i] now contains actual  
    // position of each frequency
    for (int i = 1; i <= RANGE; i++) {
        count[i] += count[i-1];     	
    }
 
    // Build the output sorted vector
    for (uint32_t i = 0; i < Freq.size() ; i++) { 
    	if (Remove[i] == 0) {
    		assert (Freq[i] <= RANGE);
    		Sortindex[count[Freq[i]] - 1] = i;
    		--count[Freq[i]];
    	}
    }  
}


void StoreIndex(string &genome, vector<GenomeTuple> &minimizers, Header &header, Options &opts) {	
	
	ifstream testGenome(genome.c_str());
	if (testGenome.good() == false or testGenome.eof()) {
		cerr << "Cannot open target " << genome << endl;
		exit(1);
	}
	gzFile f = gzopen(genome.c_str(), "r");

	kseq_t *ks = kseq_init(f);
	GenomePos offset=0;

	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		int prevMinCount = minimizers.size();
		cerr << "Storing for "<< ks->name.s << " " << prevMinCount << " " << offset << endl;
		StoreMinimizers<GenomeTuple, Tuple>(ks->seq.s, ks->seq.l, opts.globalK, opts.globalW, minimizers);
		
		for (GenomePos i=prevMinCount; i< minimizers.size(); i++) {
			minimizers[i].pos+=offset;
		}
		offset+= ks->seq.l;	
		header.Add(ks->name.s, offset);
	}
	kseq_destroy(ks);
	gzclose(f);
	cerr << "Sorting " << minimizers.size() << " minimizers" << endl;
	std::sort(minimizers.begin(), minimizers.end());

	//
	// Get the frequency for minimizers; Store the frequency in Freq;
	//
	vector<bool> Remove (minimizers.size(), 0);
	int RANGE = 0;
	vector<uint32_t> Freq(minimizers.size(), 0);

	int n = 0;
	uint32_t unremoved = 0;
	while (n < minimizers.size()) {
		int ne = n;
		while (ne < minimizers.size() && minimizers[ne].t == minimizers[n].t) {ne++;}
		int rz = 1;
		if (header.pos.back()/1000000000 > 1) {rz = header.pos.back()/1000000000;}
		if (ne - n > opts.minimizerFreq * rz) { // opts.minimizerFreq*rz is the rough threshold
			for (uint32_t i = n; i < ne; i++) {
				//Freq[i] = ne - n;
				Remove[i] = 1;
			}
		}
		else {
			//assert(ne - n != 0);
			for (uint32_t i = n; i < ne; i++) {
				Freq[i] = ne - n;
				++unremoved;
			}		
			if (ne -n > RANGE) {RANGE = ne - n;}
		}
		n = ne;
	}

	cerr << "RANGE: " << RANGE << endl;


	// Sort unremoved minimizers by frequency 
	// Use count sort
	//
	int winsize = 100;
	uint32_t sz = header.pos.back()/winsize;
	if (header.pos.back()/winsize % winsize > 0) sz += 1;
	vector<uint32_t> Sortindex(unremoved, 0);
	CountSort(Freq, RANGE, Remove, Sortindex);

	vector<uint32_t> winCount(sz, opts.NumOfminimizersPerWindow); // 50 is a parameter that can be changed 

	for (uint32_t s = 0; s < Sortindex.size(); s++) {
		uint32_t id = minimizers[Sortindex[s]].pos/winsize;
		if (winCount[id] > 0) {
			winCount[id] -= 1;
		}
		else {
			Remove[Sortindex[s]] = 1;
		}
	}

/*
	// TODO(Jingwen) Print out all the minimizer and delete this later
	if (opts.dotPlot) {
		stringstream outNameStrm;
		outNameStrm << "minimizers.txt";
		ofstream baseDots(outNameStrm.str().c_str());
		for (int m=0; m < minimizers.size(); m++) {
			baseDots << minimizers[m].t << "\t"
					 << minimizers[m].pos << "\t" 
					 << minimizers[m].pos + opts.globalK << "\t"
					 << Freq[m] << "\t" 
					 << Remove[m] << endl;				
		}
		baseDots.close();
	}
 */

	// TODO(Jingwen) Print out all the minimizer and delete this later
	if (opts.dotPlot) {
		stringstream outNameStrm;
		outNameStrm << "picked.unique.minimizers.txt";
		ofstream baseDots(outNameStrm.str().c_str());
		for (int m=0; m < minimizers.size(); m++) {
			if (Remove[m] == 0) {
				baseDots << minimizers[m].t << "\t"
						 << minimizers[m].pos << "\t" 
						 << minimizers[m].pos + opts.globalK << "\t"
						 << Freq[m] << "\t" 
						 << Remove[m] << endl;				
			}
		}
		baseDots.close();
	}
 
	//
	// Remove too frequent minimizers;
	//
	cerr << "Removing too frequent minimizers" << endl;
	RemoveFrequent (minimizers, Remove);
	cerr << "There are " << minimizers.size() << " minimizers left" << endl;

	//RemoveFrequent(minimizers, opts.globalMaxFreq);

}


int ReadIndex(string fn, vector<GenomeTuple> &index, Header &h, Options &opts) {
	ifstream fin(fn.c_str(), ios::in|ios::binary);
	if (fin.good() == false or fin.eof()) {
		return 0;
	}
	int64_t len;
	fin.read((char*) &len, sizeof(int64_t));	
	fin.read((char*) &opts.globalK, sizeof(int));
	h.Read(fin);
	index.resize(len);
	fin.read((char*) &index[0], sizeof(GenomeTuple)*len);
	return len;
}


void WriteIndex(string fn, vector<GenomeTuple> &index, Header &h, Options &opts) {
	ofstream fout(fn.c_str(), ios::out|ios::binary);
	int64_t minLength = index.size();
	fout.write((char*) &minLength, sizeof(int64_t));
	fout.write((char*) &opts.globalK, sizeof(int));
	h.Write(fout);
	fout.write((char*) &index[0], sizeof(GenomeTuple)* index.size());
	fout.close();
}

#endif
