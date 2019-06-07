#include "htslib/kseq.h"
#include "bwa/kstring.h"
#include "bwa/bwa.h"
#include <iostream>

//KSEQ_DECLARE(gzFile)
KSEQ_INIT(gzFile, gzread)
#include "Genome.h"
#include "stdlib.h"


class ProcInfo {
public: 
	bwaidx_t *idx;
	vector<int> *pos;
	vector<int> *mult;
	int window;
	int start;
	int end;
	uint8_t *binseq;
	int k;
	int c;
	Genome *genome;
};

typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;


 smem_aux_t *smem_aux_init()
{
	smem_aux_t *a;
	a = (smem_aux_t*) calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = (bwtintv_v*) calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = (bwtintv_v*) calloc(1, sizeof(bwtintv_v));
	return a;
}

 void smem_aux_destroy(smem_aux_t *a)
{	
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

void CountKmers(ProcInfo* procInfo) {
	int window=procInfo->window;
	vector<int> mult(procInfo->window,0);
	vector<int> pos(procInfo->window,0);
	int k=procInfo->k;
	smem_aux_t *a = smem_aux_init();

	for (int p=procInfo->start; p < procInfo->end-k; p++) {
		bool foundN=false;
		for (int ki=0; ki < k; ki++) {
			if (procInfo->binseq[p+ki] > 3) { foundN = true; break;}
		}
		if (foundN == false) {
			bwt_smem1(procInfo->idx->bwt, k, (uint8_t*) &procInfo->binseq[p], 0, 1, &a->mem1, a->tmpv);
			assert(p%window < mult.size());
			mult[p%window] = a->mem1.a[0].x[2];
			pos[p%window]=p;
		}
		else {
			mult[p%window] = -1;
			pos[p%window]=0;
		}
		if ((p % window) == (window-1) or p == procInfo->end-k-1) {
			//
			// Reached the end of the window
			//
			int minMult = -1;
			int minMultPos = -1;
			int minMultIdx=-1;
			for (int i=0; i < min(window, procInfo->end - (p/window) * window);i++) {

				if ((minMult < 0 and mult[i] > 0) or (minMult > 0 and minMult > mult[i])) {
					minMult = mult[i];
					minMultPos = pos[i];
					minMultIdx=i;
				}
			}

			(*procInfo->pos)[p/window] = minMultPos;
			(*procInfo->mult)[p/window] = minMult;
		}
	}
}
int main(int argc, char* argv[]) {
	int optind=1;
	bwaidx_t *idx;
	Genome genome;
	if (argc < 5) {
		cerr <<" Usage: iik genome.fasta k window maxMult" << endl;
		exit(0);
	}
	string genomeName=argv[1];
	int k=atoi(argv[2]);
	int window=atoi(argv[3]);
	int maxMult=atoi(argv[4]);


	if ((idx = bwa_idx_load(argv[optind], BWA_IDX_BWT|BWA_IDX_BNS)) == 0) return 1;
	cerr << "Done reading index" << endl;
	genome.Read(genomeName);
	cerr << "Done reading genome" << endl;
	bwtintv_t mem;
	vector<int> mult(window,0), pos(window,0);
	
	ProcInfo procInfo;
	procInfo.idx = idx;

	int nproc=48;
	for (int c =0; c < genome.seqs.size(); c++) {
		std::map<int,int> hist;
		uint8_t *binseq = new uint8_t[genome.lengths[c]];
		for (int p=0; p < genome.lengths[c] ; p++) {
			binseq[p] =  genome.seqs[c][p] < 4? genome.seqs[c][p] : nst_nt4_table[(int)genome.seqs[c][p]];
		}
		vector<int> pos(genome.lengths[c]/window, -1);
		vector<int> mult(genome.lengths[c]/window, -1);

		vector<ProcInfo> procs(nproc);
		pthread_attr_t *threadAttr = new pthread_attr_t[nproc];
		pthread_t *threads = new pthread_t[nproc];
		for (int procIndex = 0; procIndex < nproc; procIndex++ ){
			pthread_attr_init(&threadAttr[procIndex]);
		}
		
		for (int i = 0; i < nproc; i++) {
			procs[i].idx=idx;
			procs[i].pos=&pos;
			procs[i].mult=&mult;
			procs[i].k=k; procs[i].window=window;
			procs[i].start = i * (genome.lengths[c]/nproc);
			procs[i].end = min((i+1)*(genome.lengths[c]/nproc), genome.lengths[c]);
			procs[i].genome = &genome;
			procs[i].c = c;
			procs[i].binseq = binseq;
			pthread_create(&threads[i], &threadAttr[i], (void* (*)(void*)) CountKmers, &procs[i]);
		}
		for (int procIndex = 0; procIndex < nproc; procIndex++) {
			pthread_join(threads[procIndex], NULL);
		}

		/*
		for (int i=0; i < mult.size(); i++) {
			if (mult[i] > 0) {
				if (hist.find(mult[i]) == hist.end()) {
					hist[mult[i]] = 0;
				}
				hist[mult[i]]++;
			}
		}

		for (map<int,int>::iterator it=hist.begin(); it != hist.end(); ++it) {
			if (it->first < 100) {
				cerr << it->first << "\t" << it->second << endl;
			}
		}
		*/

		int i;
		i=0;
		int npos=0;
		while (i < genome.lengths[c]/window) {
			int n=i;
			long totalMult=0;
			while (n < genome.lengths[c]/window and mult[n] <= maxMult and mult[n] > 0) { totalMult+=mult[n]; n++; npos+=1;}
			if (n > i){ 
				float avgMult = ((float)totalMult)/(n-i);
				if (n < genome.lengths[c]/window) {
					cout << genome.header.names[c] << "\t" << pos[i] << "\t" << pos[n-1] << "\t" << avgMult << endl;
				}
				else {
					cout << genome.header.names[c] << "\t" << pos[i] << "\t" << pos[pos.size()-1] << "\t" << avgMult << endl;
				}
			}
			i=n+1;
		}
		cerr << "done with " << genome.header.names[c] << "\t" << npos << endl;
	}
}
