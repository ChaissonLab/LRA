#ifndef MMINDEX_H_
#define MMINDEX_H_
#include "TupleOps.h"
#include "Options.h"
#include "Genome.h"
#include "Sorting.h"
#include "MinCount.h"

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
	vector<uint64_t>    seqOffsets;
	vector<uint64_t>    tupleBoundaries;
	uint64_t offset;

	void StoreLocalIndexSize(int s) {
		if (s > 0 and s < (1 << (LOCAL_POS_BITS-1))) {
			localIndexSize = s;
		}
	}

	LocalIndex() { 
		k = 10;
		w=5; 
		offset=0;
		maxFreq=5;
		tupleBoundaries.push_back(0);
		seqOffsets.push_back(0);
		localIndexSize = min(1 << (LOCAL_POS_BITS-1), 512);
	}
	
	LocalIndex( LocalIndex &init) { 
		k=init.k;
		w=init.w;
		offset=0;
		maxFreq=init.maxFreq;
		localIndexSize = init.localIndexSize;
		tupleBoundaries.push_back(0);
		seqOffsets.push_back(0);
	}		

	void Write(string filename) {
		ofstream fout(filename.c_str(), ios::out|ios::binary);
		fout.write((char*)&k, sizeof(int));
		fout.write((char*)&w, sizeof(int));
		fout.write((char*)&localIndexSize, sizeof(int));
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
		fin.read((char*)&localIndexSize, sizeof(int));
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
			StoreMinimizers<LocalTuple, SmallTuple>(&seq[seqPos], min((GenomePos) seqLen, (GenomePos) seqPos+localIndexSize) - seqPos, k, w, locMinimizers, false );
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



void StoreIndex(string &genome, 
								vector<GenomeTuple> &minimizers, 
								Header &header,
								Options &opts) {	
	
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
	RemoveFrequent(minimizers, opts.globalMaxFreq);

	//
	// Remove too frequent minimizers;
	//
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
