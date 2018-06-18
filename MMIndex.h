#ifndef MMINDEX_H_
#define MMINDEX_H_
#include "TupleOps.h"
#include "Options.h"

void StoreIndex(string &genome, 
								vector<GenomeTuple> &minimizers, 
								Options &opts) {	

	gzFile f = gzopen(genome.c_str(), "r");

	kseq_t *ks = kseq_init(f);
	int offset=0;
	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		cerr << "Storing for "<< ks->name.s << endl;
		StoreMinimizers(ks, opts.k, opts.w, minimizers);
		for (int i=offset; i< minimizers.size(); i++) {
			minimizers[i].pos+=offset;
		}
		offset+= ks->seq.l;	
	}
	kseq_destroy(ks);
	gzclose(f);
	cerr << "Sorting " << minimizers.size() << " minimizers" << endl;
	std::sort(minimizers.begin(), minimizers.end());
	//
	// Remove too frequent minimizers;
	//
	int c=0,n=0;
	int before=minimizers.size();
	while(n < minimizers.size()) {
		int ne=n;
		while (ne < minimizers.size() && 
					 minimizers[ne].tuple == minimizers[n].tuple) { ne++;}
		if (ne - n < opts.maxFreq) {
			int end = ne;
			for (ne = n; ne < end; ne++, c++) {
				minimizers[c] = minimizers[ne];
			}
		}
		n=ne;
	}
  cerr << "Removed " << before - c << " redundant minimizers " << endl;
	minimizers.resize(c);
}

int ReadIndex(string fn, vector<GenomeTuple> &index, Options &opts) {
	ifstream fin(fn.c_str(), ios::in|ios::binary);
	if (fin.good() == false or fin.eof()) {
		return 0;
	}
	int64_t len;
	fin.read((char*) &len, sizeof(int64_t));	
	fin.read((char*) &opts.k, sizeof(int));
	InitMask(opts.k);
	index.resize(len);
	fin.read((char*) &index[0], sizeof(GenomeTuple)*len);
	return len;
}

void WriteIndex(string fn, vector<GenomeTuple> &index, Options &opts) {
	ofstream fout(fn.c_str(), ios::out|ios::binary);
	int64_t minLength = index.size();
	fout.write((char*) &minLength, sizeof(int64_t));
	fout.write((char*) &opts.k, sizeof(int));
	fout.write((char*) &index[0], sizeof(GenomeTuple)* index.size());
	fout.close();
}

#endif
