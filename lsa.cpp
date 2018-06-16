#include "htslib/hts.h"
#include "htslib/kseq.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <zlib.h>
#include <vector>
#include "htslib/kseq.h"
#include <algorithm>
#include <queue>
#include "SeqUtils.h"
#include "TupleOps.h"
#include "MinCount.h"




int k=21;


void RunStore(int argc, char* argv[]) {
	// open query file for reading; you may use your favorite FASTA/Q parser
	gzFile f = gzopen(argv[0], "r");
	assert(f);
	kseq_t *ks = kseq_init(f);
	int w=10;
	vector<GenomeTuple> minimizers, simple;
	while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence
		StoreMinimizers(ks, k, w, minimizers);
	}

}


int main(int argc, char *argv[]) {

	if (argc < 2) {
		fprintf(stderr, "Usage: lsa options\n");
		return 1;
	}
	InitSeqMap();
	InitMask(k);

  int i;

}
