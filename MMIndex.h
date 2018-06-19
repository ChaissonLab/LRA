#ifndef MMINDEX_H_
#define MMINDEX_H_
#include "TupleOps.h"
#include "Options.h"



class Header {
 public:
	vector<string> names;
	vector<uint64_t> pos;
	
	int Find(uint64_t query) { 
		vector<uint64_t>::iterator it = lower_bound(pos.begin(), pos.end(), query);
		return it - pos.begin();
	}

	void Add(const char* name, uint64_t p) {
		names.push_back(string(name));
		pos.push_back(p);
	}
		
	void Write(ofstream &out) {
		int idxLen = names.size();
		out.write((char*) &idxLen, sizeof(int));
		for(int i=0; i < names.size();i++) {
			int nameLen=names[i].size();
			out.write((char*) &nameLen, sizeof(int));
			out.write((char*) names[i].c_str(), names[i].size());
		}
		out.write((char*) &pos[0], sizeof(int64_t)*pos.size());
	}

	void Read(ifstream &in) {
		int idxLen;
		in.read((char*) &idxLen, sizeof(int));
		names.resize(idxLen);
		pos.resize(idxLen);		
		for(int i=0; i < names.size(); i++) {			
			int nameLen;
			in.read((char*) &nameLen, sizeof(int));
			char *name = new char[nameLen+1];
			name[nameLen] = '\0';
			in.read((char*) name, nameLen);
			names[i] = name;
		}
		in.read((char*) &pos[0],sizeof(int64_t)*pos.size());
	}
};

class SortByPos {
 public:
	int operator() (const GenomeTuple &a, const GenomeTuple &b) {
		return (a.pos < b.pos);
	}
};

void PrintIndex(vector<GenomeTuple> &minimizers, int k) {
	sort(minimizers.begin(), minimizers.end(), SortByPos());
	for(int i = 0; i < minimizers.size();i++) {
		string s;
		TupleToString(minimizers[i].t, k, s);
		cout << i << "\t" << minimizers[i].pos << "\t" << s << endl;
	}	
}


class GenomeLocalIndex {
 public:
	static int localIndexSize;
	int k;
	int w;
	vector<LocalTuple>  minimizers;
	vector<uint64_t>    offsets;
	vector<uint64_t>    boundaries;
	GenomeLocalIndex() { k = 8; w=4;}
	void Write(string filename) {
		ofstream fout(filename.c_str(), ios::out|ios::binary);
		fout.write((char*)&k, sizeof(int));
		fout.write((char*)&w, sizeof(int));
		int nRegions=offsets.size();
		fout.write((char*)&nRegions, sizeof(int));
		fout.write((char*)&offsets[0], sizeof(uint64_t)*offsets.size());
		fout.write((char*)&boundaries[0], sizeof(uint64_t)*boundaries.size());
		uint64_t nMin = minimizers.size();
		fout.write((char*)&minimizers[0], sizeof(LocalTuple)*minimizers.size());
	}

	void Read(string &filename) {
		ifstream fin(filename.c_str(), ios::out|ios::binary);
		fin.read((char*)&k, sizeof(int));
		fin.read((char*)&w, sizeof(int));
		int nRegions;
		fin.read((char*)&nRegions, sizeof(int));
		offsets.resize(nRegions);		
		fin.read((char*)&offsets[0], sizeof(uint64_t)*nRegions);
		boundaries.resize(nRegions);
		fin.read((char*)&boundaries[0], sizeof(uint64_t)*nRegions);
		uint64_t nMin;
		fin.read((char*) &nMin, sizeof(uint64_t));
		minimizers.resize(nMin);
		fin.read((char*)&minimizers[0], sizeof(LocalTuple)*nMin);
	}
	
	int LookupBin(uint64_t query) {
		if (offsets.size() == 0) {
			return 0;
		}
		assert(query < offsets[offsets.size()-1]);
		vector<uint64_t>::iterator it;
		it = lower_bound(offsets.begin(), offsets.end(), query);
		int index = it - offsets.begin();
		return index;
	}
	int MinimizerBounds(uint64_t query, 
	
};

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

int GenomeLocalIndex::localIndexSize=0xFFFF;

void StoreLocalIndex(string &genome, 
										 GenomeLocalIndex &glIndex,
										 Options &opts) {	
	gzFile f = gzopen(genome.c_str(), "r");

	kseq_t *ks = kseq_init(f);
	int offset=0;
	uint64_t globalOffset=0;
	Options localOptions;
	int gi=0;	
	glIndex.boundaries.push_back(0);
	glIndex.offsets.push_back(0);
	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		cerr << "Storing for "<< ks->name.s << endl;

		int nIndex = ks->seq.l / GenomeLocalIndex::localIndexSize;
		if (ks->seq.l % GenomeLocalIndex::localIndexSize != 0) {
			nIndex +=1;
		}
		int curPos=0;

		vector<LocalTuple> minimizers;
		int netSize=0;
		
		for (int i = 0; i < nIndex; i++) {	
			minimizers.clear();
			StoreMinimizers<LocalTuple, SmallTuple>(&ks->seq.s[curPos], 
																							min((int) ks->seq.l, 
																									(int) curPos+GenomeLocalIndex::localIndexSize) - curPos,
																							glIndex.k, glIndex.w, minimizers );
			std::sort(minimizers.begin(), minimizers.end());
			RemoveFrequent(minimizers, 5);
			glIndex.offsets.push_back(globalOffset+curPos);
			curPos+=GenomeLocalIndex::localIndexSize;
			netSize+=minimizers.size();
			glIndex.minimizers.insert(glIndex.minimizers.end(), minimizers.begin(), minimizers.end());
			glIndex.boundaries.push_back(glIndex.minimizers.size());
			if (gi % 100 == 0) {
				cerr << "Storing " << ks->name.s << "\t" << gi << "\t" << curPos << "\t" << minimizers.size() << "\t" << sizeof(LocalTuple) * netSize << endl;
			}
			gi++;
		}
		globalOffset+=ks->seq.l;		
	}
}


void StoreIndex(string &genome, 
								vector<GenomeTuple> &minimizers, 
								Header &header,
								Options &opts) {	

	gzFile f = gzopen(genome.c_str(), "r");

	kseq_t *ks = kseq_init(f);
	int offset=0;

	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		cerr << "Storing for "<< ks->name.s << endl;
		StoreMinimizers<GenomeTuple, Tuple>(ks->seq.s, ks->seq.l, opts.k, opts.w, minimizers);
		
		for (int i=offset; i< minimizers.size(); i++) {
			minimizers[i].pos+=offset;
		}
		header.Add(ks->name.s, offset);
		offset+= ks->seq.l;	
	}
	kseq_destroy(ks);
	gzclose(f);
	cerr << "Sorting " << minimizers.size() << " minimizers" << endl;
	std::sort(minimizers.begin(), minimizers.end());
	RemoveFrequent(minimizers, opts.maxFreq);
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
	fin.read((char*) &opts.k, sizeof(int));
	h.Read(fin);
	index.resize(len);
	fin.read((char*) &index[0], sizeof(GenomeTuple)*len);
	return len;
}

void WriteIndex(string fn, vector<GenomeTuple> &index, Header &h, Options &opts) {
	ofstream fout(fn.c_str(), ios::out|ios::binary);
	int64_t minLength = index.size();
	fout.write((char*) &minLength, sizeof(int64_t));
	fout.write((char*) &opts.k, sizeof(int));
	h.Write(fout);
	fout.write((char*) &index[0], sizeof(GenomeTuple)* index.size());
	fout.close();
}

#endif
