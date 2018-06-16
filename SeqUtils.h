#ifndef SEQ_UTILS_H_
#define SEQ_UTILS_H_
int seqMap[256];
void InitSeqMap() {
	fill(seqMap, seqMap+256,0);
	seqMap[int('a')] = 0;
	seqMap[int('A')] = 0;
	seqMap[int('c')] = 1;
	seqMap[int('C')] = 1;
	seqMap[int('g')] = 2;
	seqMap[int('G')] = 2;
	seqMap[int('t')] = 3;
	seqMap[int('T')] = 3;
}


#endif
